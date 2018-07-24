library(AER)
library(readstata13)
library(tidyverse)

# R function adapted from Ian Gows' webpage:
# http://www.people.hbs.edu/igow/GOT/Code/cluster2.R.html.
clusterVCV <- function(data, fm, cluster1, cluster2=NULL) {
   
  require(sandwich)
  require(lmtest)
   
  # Calculation shared by covariance estimates
  est.fun <- estfun(fm)
  cols <- names(fs$model)[!startsWith(names(fs$model), '(')]
  inc.obs <- complete.cases(data[,cols])
   
  # Shared data for degrees-of-freedom corrections
  N  <- dim(fm$model)[1]
  NROW <- NROW(est.fun)
  K  <- fm$rank
   
  # Calculate the sandwich covariance estimate
  cov <- function(cluster) {
    cluster <- factor(cluster)
     
    # Calculate the "meat" of the sandwich estimators
    u <- apply(est.fun, 2, function(x) tapply(x, cluster, sum))
    meat <- crossprod(u)/N
     
    # Calculations for degrees-of-freedom corrections, followed 
    # by calculation of the variance-covariance estimate.
    # NOTE: NROW/N is a kluge to address the fact that sandwich uses the
    # wrong number of rows (includes rows omitted from the regression).
    M <- length(levels(cluster))
    # dfc <- M/(M-1) * (N-1)/(N-K)
    dfc <- N/(N-K)
    dfc * NROW/N * sandwich(fm, meat=meat)
  }
   
  # Calculate the covariance matrix estimate for the first cluster.
  cluster1 <- data[inc.obs,cluster1]
  cov1  <- cov(unlist(cluster1))
   
  if(is.null(cluster2)) {
    # If only one cluster supplied, return single cluster
    # results
    return(cov1)
  } else {
    # Otherwise do the calculations for the second cluster
    # and the "intersection" cluster.
    cluster2 <- data[inc.obs,cluster2]
    cluster12 <- paste(cluster1,cluster2, sep="")
     
    # Calculate the covariance matrices for cluster2, the "intersection"
    # cluster, then then put all the pieces together.
    cov2   <- cov(cluster2)
    cov12  <- cov(cluster12)
    covMCL <- (cov1 + cov2 - cov12)
     
    # Return the output of coeftest using two-way cluster-robust
    # standard errors.
    return(covMCL)
  }
}

workfile_data <- as.tibble(read.dta13('CCL_workfile_May29.dta'))
ols_data <- workfile_data %>%
  filter(year>=2013 & year<=2015) %>%
  mutate(prefect_dummy = as.factor(prefect),
         province_year_dummy = as.factor(province_year),
         year_dummy = as.factor(year),
         province_dummy = as.factor(province))

fs <- lm(d_export_cus_pw ~ d_export_btkrow_pw00 + prefect_dummy + year_dummy,
                  ols_data, weights=ols_data$pop1564_2010)
fn <- lm(d_export_cus_pw ~ prefect_dummy + year_dummy,
                  ols_data, weights=ols_data$pop1564_2010)

waldtest(fs, fn, vcov = clusterVCV(ols_data, fs, cluster1="province_dummy"))$F[2]
