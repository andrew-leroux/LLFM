rm(list=ls())
## get the number of cores and current bootstrap sample
scenario <- as.numeric(commandArgs(trailingOnly=TRUE))
b       <- scenario[1] # current bootstrap
n_b     <- scenario[2] # total number of bootstrap samples
s_b     <- scenario[3] # current functional domain

outfile <- file.path("..","results", paste0("NHANES_LLFM_bs", b, "_s", s_b, ".rds"))


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
df_fit$X <- I(as.matrix(df_fit[,paste0("MIN",s_b)]))
## replace any NAs with 0
df_fit$X[is.na(df_fit$X)] <- 0
## get rid of redundant columns
df_fit   <- 
  df_fit %>% 
  dplyr::select(-one_of(paste0("MIN",1:1440))) 
## create active/inactive indicators
# df_fit$Y <- I(matrix(as.numeric(df_fit$X >= 100),ncol=1440,nrow=nrow(df_fit$X),byrow=FALSE))
df_fit$Y <- as.numeric(df_fit$X >= 100)
## clear up the workspace
rm(PAXINTEN, Flags, Covariate)

# subset to just ages 18-30 for computational convenince and to avoid 
# issues associated with potential non-linear effects across age 
df_fit <- 
  df_fit %>% 
  filter(Age < 30 & Age > 18 & !is.na(BMI_cat))

## create empty matrices/arrays for holding results
uid <- unique(df_fit$SEQN)
nid <- length(uid)

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
time_start <- Sys.time()

df_b$nfail <- 1-df_b$Y
fit_B <- try(glmmTMB::glmmTMB(cbind(Y, nfail) ~ Gender + Age_sc + dow_fac + day_sc + (day_sc|SEQN_b),
                              data=df_b, family=binomial))

if(!inherits(fit_B, "try-error")){
  eta_hat_B <- predict(fit_B, type="link")
  p_hat_B   <- predict(fit_B, type="response")
  
  results <- list(eta_hat = eta_hat_B,
                  p_hat = p_hat_B,
                  coef_fixed = fit_B$fit$par[1:10],
                  coef_random = as.matrix(ranef(fit_B)$cond$SEQN_b))
} else {
  results <- NULL
}

time_end <- Sys.time()
time_end-time_start

write_rds(list("results" = results, "time_start" = time_start, "time_end" = time_end, "bootstrap"=b, "s"=s_b), file=outfile)

