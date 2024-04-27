############################################################
# Section 4: selective CIs in fixed-design linear regression
############################################################
library(dplyr)
library(glmnet)
library(parallel)
library(ggplot2)

# Function to generate the data
linear_data <- function(n, p, beta, influ, lev=TRUE) {
  X <- matrix(rnorm(n = n*p, sd = sqrt(sigma2)), ncol = p)
  Sigma <- diag(rep(1, p))
  if (lev){
    # Create the leverage datapoint
    x_lev <- apply(X, 2, max)
    X <- rbind(X, influ * x_lev)
    nt <- n+1
  } else {nt <- n}
  # Generate the observed Y data
  Y <- rnorm(nt, sqrt(sigma2)) + X %*% beta
  return(
    list(X=X, Y=Y)
  )
}

# Function that defines a single simulation run
run_sim_linear_regression <- function(n, p, beta, influ, sigma2, scale=1, lev=TRUE) {
  # Generate the data
  beta <- beta * scale
  dat <- linear_data(n, p, beta, influ, lev)
  X <- dat[['X']]
  Y <- dat[['Y']]
  if (lev) {
    n <- n + 1  # update n since we have added x_lev
  }
  Sigma <- sigma2 * diag(rep(1, p))
  
  # Create containers for results
  CIs <- list()
  selected <- list()
  projected <- list()
  
  # Conduct selection/inference with data fission
  Z <- rnorm(n, sd = sqrt(sigma2))
  f_Y <- Y + Z
  g_Y <- Y - Z
  select_model <- cv.glmnet(X, f_Y, family = 'gaussian')
  selected[['fission']] <- which(coef(select_model, s = 'lambda.1se') != 0)[-1] - 1
  if (length(selected[['fission']]) > 0) {
    inference <- lm(g_Y ~ X[,selected[['fission']]])
    CIs[['fission']] <- confint(inference)[-1,]
    projected[["fission"]] <- beta[selected[["fission"]]] +
      beta[-selected[["fission"]]] %*% 
      Sigma[-selected[["fission"]], selected[["fission"]]] %*%
      solve(Sigma[selected[["fission"]], selected[["fission"]]])
  } 
  else {
    selected[['fission']] <- NA
    CIs[['fission']] <- NA
    projected[['fission']] <- NA
  }
  
  # Conduct selection/inference with data splitting (50/50 split)
  n_select <- round(n / 2)
  split_ind <- sample(
    c(rep(0, n_select), rep(1, n - n_select))
  )
  select_model <- cv.glmnet(X[split_ind == 0,], Y[split_ind == 0], family = "gaussian")
  selected[["split"]] <- which(coef(select_model, s = 'lambda.1se') != 0)[-1] - 1
  if (length(selected[["split"]]) > 0) {
    inference <- lm(Y[split_ind == 1] ~ X[split_ind == 1, selected[["split"]]])
    CIs[['split']] <- confint(inference)[-1,]
    projected[["split"]] = beta[selected[["split"]]] +
      beta[-selected[["split"]]] %*% 
      Sigma[-selected[["split"]], selected[["split"]]] %*%
      solve(Sigma[selected[["split"]], selected[["split"]]])
  }
  else {
    selected[['split']] <- NA
    CIs[['split']] <- NA
    projected[['split']] <- NA
  }
  
  # Conduct selection/inference with full dataset
  select_model <- cv.glmnet(X, Y, family = "gaussian")
  selected[["full"]] = which(coef(select_model, s = 'lambda.1se') != 0)[-1] - 1
  if (length(selected[["full"]]) > 0) {
    inference <- lm(Y ~ X[,selected[["full"]]])
    CIs[['full']] <- confint(inference)[-1,]
    projected[["full"]] = beta[selected[["full"]]] +
      beta[-selected[["full"]]] %*% 
      Sigma[-selected[["full"]], selected[["full"]]] %*%
      solve(Sigma[selected[["full"]], selected[["full"]]])
  }
  else {
    selected[['full']] <- NA
    CIs[['full']] <- NA
    projected[['full']] <- NA
  }
  
  return(
    list(CIs = CIs, selected = selected, projected = projected)
  )
}

# Now we create the visualizations

# Figure 5 (ex plot, not using experiment results)
set.seed(42)
X1 <- runif(9, 0, 2)
X2 <- 5.5
Y1 <- 1 + rnorm(length(X1))
Y2 <- X2*2 + rnorm(1)
X <- c(X1,X2)
Y <- c(Y1,Y2)

Z <- rnorm(length(X))
f_Y <- Y - Z
g_Y <- Y + Z
split <- c(rep("Split 1",length(X)/2),rep("Split 2",length(X)/2))

df <- data.frame(X = X, Y = Y, f_Y= f_Y, g_Y= g_Y,split = split[c(sample(1:length(X)))])

ggplot(df,aes(x=X,y=Y,color=split,shape=split)) +
  geom_point() +
  geom_smooth(method = "lm", fill = NA, fullrange = TRUE, linetype = "dashed") +
  theme_minimal() +
  theme(legend.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = c(0.2,0.80),
        panel.background = element_rect(fill = "white", colour = "black"),
        legend.background = element_rect(fill = "white"))
ggsave("abb_results/fig5_splitplot.pdf", width = 3, height = 3, units = "in")

ggplot(df,aes(x=X,y=f_Y)) +
  geom_point(color="blue") +
  geom_smooth(method = "lm", color = "blue", linetype ="dashed", fill = NA, fullrange=TRUE) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black")) +
  ylab("f(Y)")
ggsave("abb_results/fig5_selection.pdf", width = 3, height = 3, units = "in")

ggplot(df,aes(x=X,y=g_Y)) +
  geom_point(color="darkgreen") +
  geom_smooth(method = "lm", color = "darkgreen", linetype ="dashed", fill = NA, fullrange=TRUE) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black")) +
  ylab("g(Y)")
ggsave("abb_results/fig5_inference.pdf", width = 3, height = 3, units = "in")

# Figure 6

# Set up parameters
trials <- 50
sigma2 <- 1
p <- 100
n <- 1000
beta <- c(1, 0, rep(1,20), rep(0, p - 31), rep(-1,9))
scale <- 0.2
fig6_wrapper <- function(i) {
  res <- run_sim_linear_regression(n, p, beta, 1, sigma2, scale, lev=FALSE)
  return(res)
}
set.seed(42)
fig6_res <- mclapply(1:trials, fig6_wrapper, mc.cores = detectCores())

save(fig6_res, file="abb_results/fig6_data.Rdata")

save_dir <- 'abb_results/ci_fig_'
sel_trial <- sample(1:trials, 1)
for (type in c('fission', 'split', 'full')) {
  res <- fig6_res[[sel_trial]]
  res_sel <- res$selected[[type]]
  res_cis <- res$CIs[[type]]
  df <- data.frame(x=1:p,
                   val=beta*scale,
                   sel=rep(0,p),
                   ci_l=rep(0,p),
                   ci_u=rep(0,p),
                   cover=rep(0,p))
  df[res_sel, 'sel'] = 1
  df$sel <- as.factor(df$sel)
  df[res_sel, "ci_l"] = res_cis[,1]
  df[res_sel, "ci_u"] = res_cis[,2]
  df$cover <- as.factor((df$val >= df$ci_l) & (df$val <= df$ci_u))
  df_filtered <- df[!(df$ci_l == df$ci_u), ]
  plt <- ggplot() + 
    geom_point(data=df, aes(x, val, color=sel)) +
    geom_errorbar(data=df_filtered, aes(x, ymin=ci_l, ymax=ci_u, color=cover)) +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          legend.position = "none") +
    ylab("Coefficients") + xlab("")
  ggsave(paste0(save_dir, type, ".png"), plt, width = 4, height = 3, units = "in")
}






# Figure 7 

# Set up parameters
trials <- 500
sigma2 <- 1
n <- 15  # not counting leverage point
p <- 20
beta <- c(1, rep(0, 15), 1, -1, 1, 0)
influ_seq <- seq(2, 6, length.out = 5)  # \gamma

set.seed(42)
results <- list()
for (influ in influ_seq) {
  # Create a wrapper function so we can compute efficiently
  exp_wrapper <- function(i) {
    res <- run_sim_linear_regression(n, p, beta, influ, sigma2)
  }
  exp_res <- mclapply(1:trials, exp_wrapper, mc.cores = detectCores())
  results[[as.character(influ)]] <- exp_res
}

# Save the results
save(results, file="abb_results/regression_linear_influential.Rdata")

# Read the results
load("abb_results/regression_linear_influential.Rdata")

# Create the plot
get_metrics <- function(res, x, beta) {
  if (is.na(res$selected[[x]][1])) {
    error_FCR = 0
    CI_length = NA
    FSR = NA
    power_sign = NA
    power_selected = 0 
    precision_selected = 0
  } else if (length(res$selected[[x]]) == 1) {
    error_FCR <- mean((res$CIs[[x]][1] > res$projected[[x]]) |  (res$CIs[[x]][2] < res$projected[[x]]))
    CI_length <- mean(res$CIs[[x]][2] - res$CIs[[x]][1])
    FSR <- mean((res$projected[[x]] > 0 & res$CIs[[x]][2] < 0) | (res$projected[[x]] <0 & res$CIs[[x]][1] > 0))
    power_sign <-mean((res$projected[[x]] > 0 & res$CIs[[x]][1] > 0) | (res$projected[[x]] <0 & res$CIs[[x]][2] <0))
    power_selected <- sum(abs(beta[res$selected[[x]]]) > 0) / (sum(abs(beta)>0))
    precision_selected <- sum(abs(beta[res$selected[[x]]]) > 0) / length(res$selected[[x]])
  } else {
    error_FCR <- mean((res$CIs[[x]][,1] > res$projected[[x]]) |  (res$CIs[[x]][,2] < res$projected[[x]]))
    CI_length <- mean(res$CIs[[x]][,2] - res$CIs[[x]][,1])
    FSR <- mean((res$projected[[x]] > 0 & res$CIs[[x]][,2] < 0) | (res$projected[[x]] <0 & res$CIs[[x]][,1] > 0))
    power_sign <-mean((res$projected[[x]] > 0 & res$CIs[[x]][,1] > 0) | (res$projected[[x]] <0 & res$CIs[[x]][,2] <0))
    power_selected <- sum(abs(beta[res$selected[[x]]]) > 0) / (sum(abs(beta)>0))
    precision_selected <- sum(abs(beta[res$selected[[x]]]) > 0) / length(res$selected[[x]])
  }
  return(c(error_FCR,
           CI_length,
           FSR,
           power_sign,
           power_selected,
           precision_selected))
}

for (influ in influ_seq) {
  agg_fission <- cbind("fission",
                       influ,
                       data.frame(
                         t(sapply(results[[as.character(influ)]], function(x) get_metrics(x,"fission", beta)))
                       ))
  agg_split   <- cbind("split",
                       influ,
                       data.frame(
                         t(sapply(results[[as.character(influ)]], function(x) get_metrics(x,"split", beta)))
                       ))
  agg_full    <- cbind("full",
                       influ,
                       data.frame(
                         t(sapply(results[[as.character(influ)]], function(x) get_metrics(x, "full", beta)))
                       ))
  colnames(agg_fission) <- c("type","scale","error_FCR","CI_length","FSR","power_sign","power_selected","precision_selected")
  colnames(agg_split) <- colnames(agg_fission)
  colnames(agg_full) <- colnames(agg_fission)
  agg <- rbind(agg_fission,agg_split,agg_full)
  if(influ == influ_seq[1]){
    res_agg = agg
  }
  else{
    res_agg = rbind(res_agg,agg)
  }
}

res_agg$type = as.character(res_agg$type)
df = aggregate(res_agg$error_FCR ~ res_agg$scale + res_agg$type, FUN = mean)
colnames(df) <-c("signal","type","FCR")
FCR_plot <- df %>%
  ggplot( aes(x=signal, y=FCR, group=type, color=type)) +
  geom_line(aes(linetype = type, color = type), size = 1.5) +
  geom_point(aes(shape = type, color = type), size = 3) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 15),
        legend.position = "none", legend.text = element_text(size = 15)) +
  xlab("Leverage Parameter") +
  ylab("FCR") +
  geom_hline(yintercept=0.1)+
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0,1))

df = aggregate(res_agg$CI_length ~ res_agg$scale + res_agg$type, FUN = median)
colnames(df) <-c("signal","type","CI Length")
CI_length_plot <- df %>%
  ggplot( aes(x=signal, y=`CI Length`, group=type, color=type)) +
  geom_line(aes(linetype = type, color = type), size = 1.5) +
  geom_point(aes(shape = type, color = type), size = 3) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 15),
        legend.position = "none", legend.text = element_text(size = 15)) +
  xlab("Leverage Parameter") +
  ylab("CI Length (given selection)") 

df = aggregate(res_agg$power_selected ~ res_agg$scale + res_agg$type, FUN = mean)
colnames(df) <-c("signal","type","power_selected")
power_selected_plot <- df %>%
  ggplot( aes(x=signal, y=power_selected, group=type, color=type)) +
  geom_line(aes(linetype = type, color = type), size = 1.5) +
  geom_point(aes(shape = type, color = type), size = 3) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 15),
        legend.position = "none", legend.text = element_text(size = 15)) +
  xlab("Leverage Parameter") +
  ylab("Power Selected") +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0,1)) 

df = aggregate(res_agg$precision_selected ~ res_agg$scale + res_agg$type, FUN = mean)
colnames(df) <-c("signal","type","precision_selected")
precision_selected_plot <- df %>%
  ggplot( aes(x=signal, y=precision_selected, group=type, color=type)) +
  geom_line(aes(linetype = type, color = type), size = 1.5) +
  geom_point(aes(shape = type, color = type), size = 3) +
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 15),
        legend.position = "none", legend.text = element_text(size = 15)) +
  xlab("Leverage Parameter") +
  ylab("Precision Selected") +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0,1))

save_dir <- 'abb_results/regression_'
ggsave(paste0(save_dir, 'fcr.pdf'), FCR_plot, width = 3, height = 3, units = "in")
ggsave(paste0(save_dir, 'ci_len.pdf'), CI_length_plot, width = 3, height = 3, units = "in")
ggsave(paste0(save_dir, 'power_sel.pdf'), power_selected_plot, width = 3, height = 3, units = "in")
ggsave(paste0(save_dir, 'prec_sel.pdf'), precision_selected_plot, width = 3, height = 3, units = "in")




