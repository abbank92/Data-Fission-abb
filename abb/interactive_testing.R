############################################################
# Section 3: selective CIs in multiple testing
############################################################

library(MASS)
library(mgcv)
library(RCurl)
library(parallel)

## AdaPT
source("STAR_code/AdaPT.R")
source("STAR_code/AdaPT_gam.R")
source('STAR_code/summarize_methods.R')
source("STAR_code/expr_func.R")
source('STAR_code/STAR_convex.R')

# Set variables (temp)
args <- commandArgs(trailingOnly=TRUE)
tau <- as.numeric(args[1])
seed <- as.numeric(args[2])
print(paste("tau:", tau, "seed:", seed))
set.seed(seed)

alt <- 2
null <- 0
grid_size <- 50
n.repeat <- 100

n <- grid_size**2
x1 <- seq(-100, 100, length.out = grid_size)
x2 <- seq(-100, 100, length.out = grid_size)
x <- expand.grid(x1, x2)
colnames(x) <- c("x1", "x2")
cov.formula <- "s(x1, x2)"
alpha.list <- seq(0.01, 0.3, 0.02)
num.steps.update.score <- 10
num.steps.gam <- 5
score.params <- list(cov.formula = cov.formula,
                     num.steps = num.steps.gam)

H0 <- apply(x, 1, function(coord){sum(coord^2) < 900})
mu <- ifelse(H0, alt, null)

# Now load the scripts to use instantiated variables
source('interactive_testing_helpers.R')

# Run the experiment n.repeat times
print_suppress <- function(idx){
  print(paste("Running trial", idx))
  sink('/dev/null')
  out <- interactive_testing_exp(idx)
  sink()
  return(out)
}
wrapper_fn <- function(idx) {
  res <- tryCatch(print_suppress(idx),
                  error=function(e){
                    print(e)
                    return(list())
                  })
}

# Run the experiment once
runtime <- system.time(out <- wrapper_fn(-1))
print(runtime)
save(file=paste0("abb_results/it_results/sample_exp.Rdata"), out)

# interactive_testing_res <- mclapply(1:n.repeat, wrapper_fn, mc.cores=detectCores())
# save(file=paste0("abb_results/it_results/interactive_tau_",tau,"_seed_",seed,".Rdata"), interactive_testing_res)
