############################################################
# Section 3: selective CIs in multiple testing
# helper scripts
############################################################
library(MASS)

gen_data <- function(tau) {
  ######################################
  # Generate Gaussian noise w specific
  # correlation structure using
  # Toeplitz matrix
  ######################################
  
  # Create data and fission split
  Sigma <- diag(rep(1, n))
  Y <- mvrnorm(1, mu, Sigma)
  Z <- mvrnorm(1, rep(0, n), Sigma)
  f_Y <- Y + tau*Z
  g_Y <- Y - (1/tau)*Z
  
  # Also return p-values for fissioned and non-fissioned data
  pvals <- 1 - pnorm(Y)
  pvals_fission <- 1 - pnorm(f_Y, sqrt((1+tau**2)))
  
  return(
    list(var = Sigma,
         Y = Y,
         Z = Z, 
         f_Y = f_Y,
         g_Y = g_Y,
         pvals = pvals,
         pvals_fission = pvals_fission)
  )
}


summary.BH <- function(pvals, H0, vals, mu, var, alpha.list,
                       debug=FALSE) {
  n <- length(pvals)
  nfrej <- sapply(alpha.list, function(alpha){
    khat <- max(c(0,which(sort(pvals)<=alpha*(1:n)/n)))
    alpha <- alpha * khat / n
    sum(pvals[!H0] < alpha, na.rm = TRUE)
  })
  ntrej <- sapply(alpha.list, function(alpha){
    khat <- max(c(0,which(sort(pvals)<=alpha*(1:n)/n)))
    alpha <- alpha * khat / n
    sum(pvals[H0] < alpha, na.rm = TRUE)
  })
  true_mean <- sapply(alpha.list, function(alpha){
    khat <- max(c(0,which(sort(pvals)<=alpha*(1:n)/n)))
    alpha <- alpha * khat / n
    mean(mu[pvals < alpha], na.rm = TRUE)
  })
  measured_mean <- sapply(alpha.list, function(alpha){
    khat <- max(c(0,which(sort(pvals)<=alpha*(1:n)/n)))
    alpha <- alpha * khat / n
    mean(vals[pvals < alpha], na.rm = TRUE)
  })
  nrej <- nfrej + ntrej
  FDP <- nfrej / pmax(nrej, 1)
  power <- ntrej / max(sum(H0), 1)
  sd <- sqrt(var/nrej)
  df <- data.frame(nrej = nrej, 
                   FDP = FDP, 
                   power = power, 
                   mean=measured_mean, 
                   true_mean=true_mean,
                   sd=sd)
  return(df)
}


summary.AdaPT <- function(adapt, H0, pvals, vals, mu, var){
  nfrej <- apply(adapt$s, 2, function(s){
    bool.mask <- (pvals <= s)
    sum(bool.mask[!H0], na.rm = TRUE)
  })
  ntrej <- apply(adapt$s, 2, function(s){
    bool.mask <- (pvals <= s)        
    sum(bool.mask[H0], na.rm = TRUE)
  })
  true_mean <- apply(adapt$s, 2, function(s){
    bool.mask <- (pvals <= s)        
    mean(mu[bool.mask], na.rm = TRUE)
  })
  measured_mean <- apply(adapt$s, 2, function(s){
    bool.mask <- (pvals <= s)        
    mean(vals[bool.mask], na.rm = TRUE)
  })
  
  nrej <- nfrej + ntrej
  FDP <- nfrej / pmax(nrej, 1)
  power <- ntrej / max(sum(H0),1)
  sd <- sqrt(var/nrej)
  df <- data.frame(nrej = nrej, FDP = FDP, power = power, mean=measured_mean, true_mean =true_mean,sd=sd)
  return(df)
}


summary.STAR <- function(STAR.obj, H0, vals, mu, var){
  nfrej <- apply(STAR.obj$mask, 2, function(x){
    sum(x[!H0], na.rm = TRUE)
  })
  ntrej <- apply(STAR.obj$mask, 2, function(x){
    sum(x[H0], na.rm = TRUE)
  })
  true_mean <- apply(STAR.obj$mask, 2, function(x){
    mean(mu[x])
  })
  measured_mean <- apply(STAR.obj$mask, 2, function(x){
    mean(vals[x])
  })
  nrej <- nfrej + ntrej
  FDP <- nfrej / pmax(nrej, 1)
  power <- ntrej / max(sum(H0),1)
  sd <- sqrt(var/nrej)
  df <- data.frame(nrej = nrej, FDP = FDP, power = power, mean=measured_mean, true_mean =true_mean, sd=sd)
  return(df)
}


interactive_testing_exp <- function(idx, debug=FALSE) {
  if (debug) print(paste('Exp run', idx))
  
  # Generate the data
  out <- gen_data(tau)
  for (param in names(out)) {assign(param, out[[param]])}
  
  # Reject with BH, AdaPT, STAR for full dataset
  STAR.obj1 <- STAR.convex(pvals, x,
                           alpha.list = alpha.list,
                           type = "model-assist",
                           update.mask.params = list(dir="min"),
                           num.steps.update.score = num.steps.update.score,
                           score.params = score.params)
  AdaPT.obj1 <- AdaPT(x, pvals, cov.formula = cov.formula,
                      q.list = alpha.list, plot.quiet = TRUE)
  
  # Compute the results for the three CI procedures for the full dataset
  BH.result_full <- summary.BH(pvals, H0, Y, mu, 1, alpha.list, debug=debug)
  STAR.result_full <- summary.STAR(STAR.obj1, H0, Y, mu, 1)
  adapt.result_full <- summary.AdaPT(AdaPT.obj1, H0, pvals, Y, mu, 1)
  
  #Reject with BH, AdaPT, STAR for fissioned data
  STAR.obj2 <- STAR.convex(pvals_fission, x, 
                           alpha.list = alpha.list,
                           type = "model-assist",
                           update.mask.params = list(dir = "min"),
                           num.steps.update.score = num.steps.update.score,
                           score.params = score.params)
  AdaPT.obj2 <- AdaPT(x, pvals_fission, cov.formula = cov.formula,
                      q.list = alpha.list, plot.quiet = TRUE)
  
  # Compute the results for the three CI procedures for the fissioned dataset
  var_fission <- 1 + (1/(tau**2))
  BH.result.fission <- summary.BH(pvals_fission, H0, g_Y, mu, var_fission, alpha.list)
  STAR.result.fission <- summary.STAR(STAR.obj2, H0, g_Y, mu, var_fission)
  adapt.result.fission <- summary.AdaPT(AdaPT.obj2, H0, pvals, g_Y, mu, var_fission)
  
  return(list(x=x,mu=mu,tau=tau,H0=H0,mu=mu,alt=alt,null=null,type=type,var=var,var_fission=var_fission,
              Y=Y, f_Y=f_Y,g_Y=g_Y,pvals=pvals,pvals_fission=pvals_fission,
              BH.result.fission=BH.result.fission, BH.result.full=BH.result_full,
              adapt.result.fission=adapt.result.fission,adapt.result.full=adapt.result_full,
              STAR.result.fission=STAR.result.fission, STAR.result.full = STAR.result_full,
              STAR.full.object=STAR.obj1,AdaPT.full.object=AdaPT.obj1,
              STAR.mask.object=STAR.obj2,AdaPT.mask.object=AdaPT.obj2))
}


