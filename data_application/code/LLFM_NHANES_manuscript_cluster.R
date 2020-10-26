rm(list=ls())
## get the number of cores and current bootstrap sample
scenario <- as.numeric(commandArgs(trailingOnly=TRUE))
b       <- scenario[1] # current bootstrap
n_b     <- scenario[2] # total number of bootstrap samples
n_cores <- scenario[3] # number of cores for parallelizing fit along the functional domain
outfile <- file.path("..","results", paste0("NHAMES_LLFM_bs", b, ".rds"))



## Check for packages needed to run analyses/install the rnhanesdata package.
pckgs <- c("devtools","mgcv","lme4","glmmTMB","tidyverse","parallel","doParallel","foreach")
sapply(pckgs, function(x) if(!require(x,character.only=TRUE,quietly=TRUE)) {
  install.packages(x)
  require(x, character.only=TRUE)
})
rm(list=c("pckgs"))
## Install the rnhanesdata package and dependencies.
## This may take a few minutes because of the size of the data package.
if(!require("rnhanesdata")){
  install_github("andrew-leroux/rnhanesdata",build_vignettes = FALSE)
  require("rnhanesdata")
}
## load activity count and wear/non-wear flag data
data("PAXINTEN_C"); data("PAXINTEN_D")
data("Flags_C"); data("Flags_D")
data("Covariate_C"); data("Covariate_D")

## re-code activity counts which are considered "non-wear" to be 0
## this doesn't impact much data, most estimated non-wear times correspond to 0 counts anyway
PAXINTEN_C[,paste0("MIN",1:1440)] <- PAXINTEN_C[,paste0("MIN",1:1440)]*Flags_C[,paste0("MIN",1:1440)]
PAXINTEN_D[,paste0("MIN",1:1440)] <- PAXINTEN_D[,paste0("MIN",1:1440)]*Flags_D[,paste0("MIN",1:1440)]


## mege 2003-2004 and 2005-2006 waves' data
PAXINTEN <- bind_rows(PAXINTEN_C, PAXINTEN_D)
Flags    <- bind_rows(Flags_C, Flags_D)
Covariate <- bind_rows(Covariate_C, Covariate_D) %>% 
  mutate(Age = RIDAGEEX/12) %>% 
  dplyr::select(SEQN, Gender, BMI_cat, Age)

## clear up the workspace
rm(PAXINTEN_C, PAXINTEN_D, Flags_C, Flags_D, Covariate_C, Covariate_D)

## 1) subset to good days of data (>= 10 hours estiamted wear time) and good quality data indicators (PAXCAL/PAXSTAT)
## 2) create day number of longitudinal information
## 3) re-code day of the week as factor variable (for ease of use in regressions)
PAXINTEN$nmins_wear <- rowSums(Flags[,paste0("MIN",1:1440)], na.rm=TRUE)
PAXINTEN <- 
  PAXINTEN %>% 
  mutate(good_day = as.numeric(nmins_wear >= 600 & PAXCAL %in% 1 & PAXSTAT %in% 1)) %>% 
  group_by(SEQN) %>% 
  mutate(day = 1:n(),
         n_good_days = sum(good_day)) %>% 
  ungroup() %>% 
  filter(good_day == 1, n_good_days >= 3) %>% 
  dplyr::select(-PAXCAL, -PAXSTAT) %>% 
  mutate(dow_fac = factor(WEEKDAY, levels=1:7, labels=c("Sun","Mon","Tue","Wed","Thu","Fri","Sat")))



## merge PA and covariate data
df_fit <- left_join(PAXINTEN, Covariate, by="SEQN")
## add PA to the data frame as a matrix
df_fit$X <- I(as.matrix(df_fit[,paste0("MIN",1:1440)]))
## replace any NAs with 0
df_fit$X[is.na(df_fit$X)] <- 0
## get rid of redundant columns
df_fit   <- 
  df_fit %>% 
  dplyr::select(-one_of(paste0("MIN",1:1440))) 
## create active/inactive indicators
df_fit$Y <- I(matrix(as.numeric(df_fit$X >= 100),ncol=1440,nrow=nrow(df_fit$X),byrow=FALSE))
## clear up the workspace
rm(PAXINTEN, Flags, Covariate)
gc()


### bin the active/inactive indicators
# binning window
tlen <- 1
nt   <- floor(1440/tlen)
tind <- seq(0,1,len=nt)
# create a list of indices for binning into tlen minute windows
inx_col_ls       <- split(1:1440, rep(1:nt,each=tlen))
df_fit$Y_bin     <- I(sapply(inx_col_ls, function(x) rowSums(df_fit$Y[,x,drop=FALSE])))
df_fit$X_bin     <- I(sapply(inx_col_ls, function(x) rowSums(df_fit$X[,x,drop=FALSE])))
# subset to just ages 18-30 for computational convenince and to avoid 
# issues associated with potential non-linear effects across age 
df_fit <- 
  df_fit %>% 
  filter(Age < 30 & Age > 18 & !is.na(BMI_cat)) 

## functional domain 
sind <- 1:ncol(df_fit$Y_bin)

## create empty matrices/arrays for holding results
uid <- unique(df_fit$SEQN)
nid <- length(uid)
coef_fixed_mat_BB  <- coef_fixed_mat_B  <- matrix(NA, ncol=10, nrow=length(sind))     # fixed effects
coef_random_mat_BB <- coef_random_mat_B <- array(NA, dim=c(length(sind), nid, 2))    # random effects (here only 2 per subject b/c we fit a random intercept + slope model)
eta_mat_BB <- eta_mat_B <-
  p_mat_BB <- p_mat_B <-matrix(NA, nrow= nrow(df_fit), ncol=length(sind)) # linear predictors

## scale and center continuous variables (helps with convergence issues)
df_fit$Age_sc <- scale(df_fit$Age)
df_fit$day_sc <- scale(df_fit$day)
## expit function
expit <- function(x) 1/(1+exp(-x))


## simulate common set of bootstrap indices for reproducibility
set.seed(91010)
# get subject ids
bs_id <-  matrix(c(uid,sample(uid, size=nid*(n_b-1), replace=TRUE)), nid, n_b, byrow=FALSE)
# create re-sampled data fram current subject bootstraps
bs_id_b <- bs_id[,b]

df_b        <- lapply(seq_along(bs_id_b), function(x) data.frame("SEQN"=bs_id_b[x],"SEQN_b"=x))
df_b        <- lapply(df_b, function(x) data.frame("SEQN_b"=x$SEQN_b, df_fit[which(df_fit$SEQN == x$SEQN),]))
df_b        <- bind_rows(df_b)



## fit the models
time_st <- Sys.time()
## set up cluster for parallel computing
cl <- makeCluster(n_cores,setup_strategy = "sequential")
registerDoParallel(cl)


time_start <- Sys.time()
results <- foreach(s = sind, .packages=c("glmmTMB", "tidyverse")) %dopar% {
  df_fit_s <- df_b %>% dplyr::select(-X,-Y,-Y_bin) %>% data.frame("Y" = df_b$Y_bin[,s])
  
  df_fit_s$nfail <- tlen - df_fit_s$Y
  
  time_st <- Sys.time()
  fit_B <- try(glmmTMB::glmmTMB(cbind(Y, nfail) ~ Gender + Age_sc + dow_fac + day_sc + (day_sc|SEQN_b),
                       data=df_fit_s, family=binomial))
  time_ed <- Sys.time()
  time_ed-time_st

  if(!inherits(fit_B, "try-error")){
    eta_hat_B <- predict(fit_B, type="link")
    p_hat_B   <- predict(fit_B, type="response")
    
    ret <- list(eta_hat = eta_hat_B,
                p_hat = p_hat_B,
                coef_fixed = fit_B$fit$par[1:10],
                coef_random = as.matrix(ranef(fit_B)$cond$SEQN))
    return(ret)
  } else {
    return(NULL)
  }
}
time_end <- Sys.time()
time_end-time_start
# stopCluster(cl)

write_rds(list("results" = results, "time_st" = time_start, time_end = time_end), file=outfile)


## save results
# 
# 
# 
# 
# beta_hat <- t(vapply(results, function(x) as.vector(x$coef_fixed), numeric(10)))
# 
# 
# 
# beta_hat_sc <- beta_hat
# beta_hat_sc[,3] <- beta_hat_sc[,3]/sd(df_fit$Age)
# beta_hat_sc[,10] <- beta_hat_sc[,10]/sd(df_fit$day)
# # tind <- seq(0,1,len=1440)
# beta_tilde <- apply(beta_hat_sc, 2, function(x) gam(x ~ s(tind, bs="cr", k = 15))$fitted.values)
# colnames(beta_hat) <- c("Intercept","Gender (Female)", "Age", 
#                         "DoW:Mon","DoW:Tue","DoW:Wed","DoW:Thu","DoW:Fri","DoW:Sat",
#                         "Day")
# 
# 
# ylims <- c(-2,3)
# textsize <- 2
# 
# xinx <- (c(1,6,12,18,23)*60+1)/1440
# xinx_lab <- c("01:00","06:00","12:00","18:00","23:00")
# jpeg("~/Desktop/coefs_NHANES_llfm.jpeg",height=1600,width=1800,quality=100)
# par(mfrow=c(4,5),oma=c(7,7,0,0))
# matplot(tind, cbind(beta_hat_sc[,1], beta_tilde[,1]),type='l',lty=c(1,1),lwd=2*textsize,xlab="",ylab="",
#         main=colnames(beta_hat)[1],xaxt='n',cex.axis=textsize,cex.lab=textsize,cex.main=textsize)
# axis(1,at=xinx, xinx_lab,cex.axis=textsize)
# legend("bottomright",c("Raw Estimate","Smoothed Estimate"), col=c("black","red"), lty=c(1,1), lwd=2*textsize, cex=textsize, bty='n')
# 
# 
# for(p in c(2,3,10,4:9,11:16)){
#   if(p %in% c(3,10)){
#     ylims <- c(-0.2,0.3)
#   } else {
#     ylims <- c(-2,3)
#   }
#   matplot(tind, cbind(beta_hat_sc[,p], beta_tilde[,p]),type='l',lty=c(1,1),lwd=2*textsize,xlab="",ylab="",
#           main=colnames(beta_hat)[p], ylim=ylims,xaxt='n',cex.axis=textsize,cex.lab=textsize,cex.main=textsize)
#   axis(1,at=xinx, xinx_lab,cex.axis=textsize)
#   abline(h=0,col='grey',lwd=textsize)
#   
# }
# # matplot(tind, cbind(beta_hat[,2], beta_tilde[,2]),type='l',lty=c(1,2),lwd=textsize,xlab="Time of Day",ylab="",
# #         main=colnames(beta_hat)[2], ylim=ylims)
# # abline(h=0,col='grey',lwd=textsize)
# # 
# # matplot(tind, cbind(beta_hat[,3], beta_tilde[,3]),type='l',lty=c(1,2),lwd=textsize,xlab="Time of Day",ylab="",
# #         main=colnames(beta_hat)[3], ylim=ylims)
# # abline(h=0,col='grey',lwd=textsize)
# # 
# # matplot(tind, cbind(beta_hat[,10], beta_tilde[,10]),type='l',lty=c(1,2),lwd=textsize,xlab="Time of Day",ylab="",
# #         main=colnames(beta_hat)[10], ylim=ylims)
# # abline(h=0,col='grey',lwd=textsize)
# 
# mtext(text="Time of Day (t)", side=1, cex=textsize,outer=T,line=1)
# mtext(text=expression(hat(beta)(t)), side=2, cex=textsize,outer=T,line=2)
# 
# dev.off()
# 
# 
# 
# 
# 
# 
# for(s in seq_along(sind)){
#   df_fit_s <- df_fit %>% dplyr::select(-X,-Y,-Y_bin) %>% data.frame(., "Y" = df_fit$Y_bin[,s])
#   # df_fit_s <- df_fit %>% dplyr::select(-X,-Y,-Y_bin) %>% data.frame(., "Y" = df_fit$X_bin[,s])
#   
#   df_fit_s$SEQN_fac <- as.factor(df_fit_s$SEQN)
#   df_fit_s$nfail <- tlen - df_fit_s$Y
#   
#   fit_B <- try(glmmTMB(cbind(Y, nfail) ~ Gender + Age_sc + dow_fac + day_sc + (day_sc|SEQN),
#                    data=df_fit_s, family=binomial))
# 
#   if(!inherits(fit_B, "try-error")){
#     eta_hat_B <- predict(fit_B, type="link")
#     p_hat_B   <- predict(fit_B, type="response")
#     
#     # eta_mat_B[,s] <- eta_hat_B
#     # p_mat_B[,s]   <- p_hat_B
#     # 
#     # coef_fixed_mat_B[s,]  <- fit_B$fit$par[1:10]
#     # coef_random_mat_B[s,,] <- as.matrix(ranef(fit_B)$cond$SEQN)
#     
#     ret <- list(eta_hat = eta_hat_B,
#                 p_hat = p_hat_b,
#                 coef_fixed = fit_B$fit$par[1:10],
#                 coef_random = as.matrix(ranef(fit_B)$cond$SEQN))
#     return(ret)
#   } 
#   
#   # fit_BB <- try(glmmTMB(cbind(Y, nfail) ~ Gender + Age_sc + dow_fac + day_sc + (day_sc|SEQN),
#   #                       data=df_fit_s, family=betabinomial))
#   # if(!inherits(fit_BB, "try-error")){
#   # 
#   #   eta_hat_BB <- predict(fit_BB, type="link")
#   #   p_hat_BB   <- predict(fit_BB, type="response")
#   # 
#   #   eta_mat_BB[,s] <- eta_hat_BB
#   #   p_mat_BB[,s]   <- p_hat_BB
#   # 
#   #   coef_fixed_mat_BB[s,]  <- fit_BB$fit$par[1:10]
#   #   coef_random_mat_BB[s,,] <- as.matrix(ranef(fit_BB)$cond$SEQN)
#   # }
#   # 
#   
#   
#   # print(mean((p_mat_BB[,s]- df_fit_s$Y/tlen)^2))
#   print(mean((p_mat_B[,s]- df_fit_s$Y/tlen)^2))
#   print(s)
# }
# time_end <- Sys.time()
# 
# time_end - time_st
# 
# 
# ylims <- range(c(coef_fixed_mat_B, coef_fixed_mat_BB))
# xind <- seq(0,1,len=144)
# xinx <- (c(1, 6, 12, 18, 23)*60 + 1)/1440
# xinx_lab <- c("01:00","06:00","12:00","18:00","23:00")
# textsize <- 2
# jpeg("~/Desktop/BetaBinomial_Coefs.jpeg",height=500,width=1200,quality=100)
# par(mfrow=c(1,2), las=1, mar=c(5,8,5,2))
# matplot(xind, coef_fixed_mat_B, xlab="Time of day",ylab=expression(hat(beta(t))),main="Binomial Mixed Model",ylim=ylims, type='l',xaxt='n',
#         cex.main=textsize,cex.lab=textsize,cex.axis=textsize)
# axis(1, at=xinx, labels=xinx_lab, cex.axis=textsize)
# matplot(xind, coef_fixed_mat_BB, xlab="Time of day",ylab=expression(hat(beta(t))),main="Beta-Binomial Mixed Model",ylim=ylims, type='l',xaxt='n',
#         cex.main=textsize,cex.lab=textsize,cex.axis=textsize)
# axis(1, at=xinx, labels=xinx_lab, cex.axis=textsize)
# legend("bottomright", c("Intercept","Gender (Female)","Age",
#                         paste0("DoW:", c("Mon","Tue","Wed","Thu","Fri","Sat")),
#                         "Day"), 
#        bty='n',cex=textsize*0.75,lwd=textsize,lty=1, col=1:10,ncol=2)
# dev.off()
# 
# textsize <- 3
# jpeg("~/Desktop/BetaBinomial_Subjs.jpeg",height=1200,width=1600,quality=100)
# par(mfrow=c(3,4))
# set.seed(100)
# inx <- sample(1:nrow(p_mat_B), size=12,replace=F)
# for(i in 1:12){
#   plot(xind, df_fit$Y_bin[inx[i],]/tlen, pch=16, col=rgb(0,0,0,1),ylim=c(0,1), cex=0.25*textsize, 
#        cex.main=textsize,cex.lab=textsize,cex.axis=textsize,ylab=expression(Pr(Active|X,b,t)))
#   lines(xind, p_mat_B[inx[i],], col=rgb(1,0,0,0.5),pch=16, cex=0.25*textsize)
#   lines(xind, p_mat_BB[inx[i],], col=rgb(0,0,1,0.5),pch=16, cex=0.25*textsize)
#   if(i == 1){
#     legend("topleft", c("Empirical","Binomial","Beta-Binomial"),  pch=16, col=c("black","red","blue"), cex=textsize*0.5, bty='n')
#   }
# }
# dev.off()
# 
# 
# 
# 
# ### plot (smoothed) \hat{\beta}
# ## note that age and day have been standardized, need to multiply by a constand (standard deviation across the data) to obtain 
# ## inference on the original scale (i.e. units of years/days, respectively)
# coef_fixed_mat_plt <- coef_fixed_mat
# rescale_age_day <- TRUE
# if(rescale_age_day){
#   coef_fixed_mat_plt[,3]  <- coef_fixed_mat_plt[,3]*sd(df_fit$Age)
#   coef_fixed_mat_plt[,10] <- coef_fixed_mat_plt[,10]*sd(df_fit$day)
# }
# ## Note that here I'm not plotting the intercept function just to make 
# ## comparisons across coefficients easier visually (range of y-axis increases substantially with the intercept)
# xlab_inx <- (c(1,6,12,18,23)*60 + 1)/(1440/length(sind))
# xlab     <- c("01:00","06:00","12:00","18:00","23:00")
# labels <- c("Intercept","Gender","Age","DOW: Monday","DOW: Tuesday","DOW:Wednesday","DOW:Thursday","DOW:Friday","DOW:Saturday","Day")
# ylims <- range(as.vector(coef_fixed_mat_plt[,c(2:10)]))
# par(mfrow=c(3,3))
# for(i in 2:10){
#   plot(sind, gam(coef_fixed_mat_plt[,i] ~ s(sind,k=10), method="REML")$fitted.values,type='l',main=labels[i],xlab="Time of Day", ylab=expression(hat(beta)),
#        xaxt='n',ylim=ylims,lwd=2)
#   points(sind, coef_fixed_mat_plt[,i], pch=16,col='red')
#   abline(h=0,col='grey',lty=2,lwd=2)
#   axis(1, at=xlab_inx, labels=xlab)
# }
# 
# 
# 
# ### plot (smoothed) \eta_ij
# set.seed(101)
# nplt <- 3
# nplt <- length(unique(uid))
# ## plot smoothed \eta_ij directly in the top row
# ## plot smoohted \eta_ij obtained from smoothing the \beta's in the bottom row
# par(mfrow=c(1,nplt))
# layout(matrix(1:(2*nplt),ncol=nplt,nrow=2,byrow=FALSE))
# uid_plt   <- sample(uid, size=nplt, replace=FALSE)
# inx_uid   <- which(df_fit$SEQN %in% uid_plt)
# ylims_eta <- range(eta_mat[inx_uid,])
# 
# ## get smoothed \betas
# for(k in 1:10){
#   coef_fixed_sm <- apply(coef_fixed_mat, 2, function(x) gam(x ~ s(sind, k=20), method="REML")$fitted.values)
# }
# ## get smoothed b_i(s)
# ## note that you can either smooth the intercept and slope separately, or smooth the sum of the two.
# ## here i smooth the sum
# ## Note that this is quite slow for all subjects (could be made much faster with some clever matrix algebra)
# ## so i only calculate Zmat for those who we are interested in plotting
# Xmat <- (model.matrix(fit) %*% t(coef_fixed_sm)) 
# Zmat <- matrix(NA, ncol=length(sind), nrow=nrow(df_fit))
# for(i in inx_uid){
#   inx_i <- which(uid == df_fit$SEQN[i])
#   Zmat[i,] <- gam(coef_random_mat[,inx_i,] %*% c(1, df_fit$day_sc[i])~ s(sind,k=10), metho="REML")$fitted.values
#   if(i %% 100 == 0) print(i)
# }
# eta_mat_sm_beta <- Xmat + Zmat
# for(i in 1:nplt){
#   Ji <- sum(df_fit$SEQN==uid_plt[i])
#   inx_i <- which(df_fit$SEQN == uid_plt[i])
#   ## get smoothed \eta_ij via direct smoothing
#   plot(sind, rep(-1000,length(sind)), xaxt='n',xlab="Time of Day", ylab=expression(eta[ij]),
#        xaxt='n',ylim=ylims_eta,lwd=2,main=paste0("SEQN",uid_plt[i]))
#   abline(h=0,col='grey',lty=2,lwd=2)
#   axis(1, at=xlab_inx, labels=xlab)
#   for(j in 1:Ji){
#     lines(sind, gam(eta_mat[inx_i[j],]~ s(sind, k=10), method="REML")$fitted.values,col=df_fit$day[inx_i[j]])
#   }
#   ## results from smoothing \beta first 
#   plot(sind, rep(-1000,length(sind)), xaxt='n',xlab="Time of Day", ylab=expression(eta[ij]),
#        xaxt='n',ylim=ylims_eta,lwd=2,main=paste0("SEQN",uid_plt[i]))
#   abline(h=0,col='grey',lty=2,lwd=2)
#   axis(1, at=xlab_inx, labels=xlab)
#   for(j in 1:Ji){
#     lines(sind, eta_mat_sm_beta[inx_i[j],],col=df_fit$day[inx_i[j]])
#   }
#   if(i == 1){
#     legend("bottom",paste0("Day ",1:7), col=1:7,lwd=2,bty='n')
#   }
# }
# 
# 
