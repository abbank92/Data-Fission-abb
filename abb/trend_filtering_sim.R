############################################################
# Section 6: trend filtering for spectroscopy
############################################################

library(cplm)
library(trendfiltering)

# Define our function to solve for the knots given a fitted sequence
find_the_knots <- function(fit_Y, X, deg=1, tol=0.01) {
  # Compute difference matrices
  first_deriv <- diff(fit_Y)/diff(X)
  second_deriv <- diff(first_deriv) / diff(X)[2:(length(X)-1)]
  if (deg == 1) {
    knots <- which(abs(second_deriv) > tol) + 2
  } else if (deg == 2) {
    third_deriv <- diff(second_deriv)/diff(X)[3:(length(X)-1)]
    knots <- which(abs(third_deriv) > tol) + 3
  }
  return(knots)
}

# Define our function to compute uniform CIs fro trend filtering
# Per Koenker, R (2011) "Additive models for quantile regression..."
unif_CI <- function(n, bs_X, X, Y, alpha){
  E <- eigen(solve(t(bs_X) %*% bs_X))
  B <- E$vectors %*% diag(sqrt(pmax(0, E$values))) %*% t(E$vectors)
  BX1 <- B %*% t(bs_X[-1,])
  BX1 <- BX1 / sqrt(apply(BX1^2, 2, sum))
  BX0 <- B %*% t(bs_X[-n, ])
  BX0 <- BX0 / sqrt(apply(BX0^2, 2, sum))
  kappa <- sum(sqrt(apply((BX1 - BX0)^2, 2, sum)))
  v <- n - k - 1
  cvu <- uniroot(
    function(x) {kappa*(1 + x^2/v)^(-v/2)/pi + 2*(1 - pt(x, df = v)) - alpha},
    c(0, max(2*kappa/alpha/pi, 10)))$root
  
  true_trend <- lm(X ~ bs_X)$fitted.values
  pred_trend <- predict(
    lm(Y ~ bs_X), se.fit=TRUE, level = 1 - alpha
  )
  width <- cvu * pred_trend$se.fit
  
  return(list(predicted = cbind(pred_trend$fit, pred_trend$fit - width, pred_trend$fit + width, true_trend),
              c = cvu, 
              mean_se = mean(pred_trend$se.fit)))
}


# Define our function to compute the trend filtering
filt_the_trend <- function(astro_obj, type='SURE', onesd_rule=TRUE, deg=1, alpha=0.1){
  # Split out features and other required vars
  X <- astro_obj$wavelength
  Y <- astro_obj$flux
  std <- astro_obj$flux_std_err
  n <- length(X)
  noise <- rnorm(n, sd = std)
  f_Y <- Y + noise
  g_Y <- Y - noise
  
  # Perform the trendfiltering
  if (type=='SURE'){
    tf <- sure_trendfilter(X, f_Y, weights=1/std, k=deg)
    lambda_min <- tf$lambda_min
    lambda_1se <- tf$lambda_1se
  } else {  # use cv_trendfilter
    tf <- cv_trendfilter(X, f_Y, weights=1/std, k=deg)
    lambda_min <- tf$lambda_min[3]
    lambda_1se <- tf$lambda_1se[3]
  }
  
  # Select knots based on one SD rule or min lambda 
  if (onesd_rule) {
    fit_Y <- tf$fitted_values[, which(tf$lambda == lambda_1se)]
  } else {
    fit_Y <- tf$fitted_values[, which(tf$lambda == lambda_min)]
  }
  knots <- find_the_knots(fit_Y, X, deg = deg)
  
  # Compute the basis
  k <- length(knots) + deg
  basis <- tp(1:n, knots = knots, degree =  deg, k=k)
  bs_X <- cbind(rep(1, n), basis$X, basis$Z)
  n_x <- ncol(bs_X)
  P_CI <- predict(
    lm(g_Y ~ bs_X), 
    interval="confidence", 
    se.fit = TRUE, 
    weights = 1/std,
    level = 1 - alpha)
  U_CI <- unif_CI(
    n = n, 
    bs_X = bs_X,
    X = 1:n, 
    Y = g_Y,
    alpha = alpha)$predicted[,2:3]
  
  # Create table of results
  df <- data.frame(cbind(X, Y, P_CI$fit, U_CI))
  colnames(df) <- c("Position",
                    "Actual",
                    "Predicted",
                    "CI_Low",
                    "CI_High",
                    "UCI_Low",
                    "UCI_High")
  return(df)
}

quasar <- readRDS("trendfiltering-astro-data/fig3top_quasar_spectrum.rds")
galaxy <- readRDS("trendfiltering-astro-data/fig3middle_galaxy_spectrum.rds")
stellar <- readRDS("trendfiltering-astro-data/fig3bottom_stellar_spectrum.rds")

# Filter to the wavelengths we care about
max_wvln <- 5500
min_wvln <- 4000

quasar <- quasar[quasar$wavelength <= max_wvln ,]
quasar <- quasar[quasar$wavelength >= min_wvln ,]

galaxy <- galaxy[galaxy$wavelength <= max_wvln ,]
galaxy <- galaxy[galaxy$wavelength >= min_wvln ,]

stellar <- stellar[stellar$wavelength <= max_wvln ,]
stellar <- stellar[stellar$wavelength >= min_wvln ,]

# Fit the model
quasar_fit <- filt_the_trend(quasar)
galaxy_fit <- filt_the_trend(galaxy)
stellar_fit <- filt_the_trend(stellar)


# Now let's reproduce the plots!
make_space_plot <- function(space_fit, plot_title='hi there', xlim = c(4000, 5000), ylim = c(0, 50)) {
  p <- ggplot(space_fit) +
    geom_line(aes(Position, Actual, linetype='Actual')) +
    geom_line(aes(Position, Predicted, linetype='Fitted')) + 
    geom_line(aes(Position, UCI_Low), color = 'orange', linetype='longdash') +
    geom_line(aes(Position, UCI_High, color='Uniform CI'), linetype='longdash') +
    geom_line(aes(Position, CI_Low), color = 'blue', linetype='twodash') +
    geom_line(aes(Position, CI_High, color='Pointwise CI'), linetype='twodash') +
    # theme_minimal() +
    theme(
      legend.position = 'bottom',
      panel.background = element_rect(fill = "white", colour = "black"),
      panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
      panel.grid.minor = element_line(colour = "grey"),
    ) +
    xlab(TeX("$\\lambda (Å)")) +
    ylab(TeX("$f(\\lambda) ( 10^{-7} ergs/s/cm^{2}/Å ) $")) +
    ggtitle(plot_title) + 
    scale_color_manual(NULL, values = c('blue','orange')) +
    scale_linetype_manual(NULL, values = c('solid','dashed')) +
    facet_zoom(
      xlim = xlim,
      ylim = ylim,
      zoom.size = 0.5,
      show.area = TRUE,
      horizontal = TRUE,
      split = FALSE
    )
  return(p)
}
quasar_p <- make_space_plot(quasar_fit,
                            plot_title = parse(text= paste("Quasar ( ~", expression(z %~~% 2.4),")",sep="")),
                            xlim = c(4700, 4800),
                            ylim = c(3, 13))
ggsave('abb/abb_results/quasar_fit.png', quasar_p)

galaxy_p <- make_space_plot(galaxy_fit,
                            plot_title = parse(text= paste("Galaxy ( ~", expression(z %~~% 2.4),")",sep="")),
                            xlim = c(4700, 4800),
                            ylim = c(5, 9))
ggsave('abb/abb_results/galaxy_fit.png', galaxy_p)

stellar_p <- make_space_plot(stellar_fit,
                             plot_title = parse(text= paste("Star ( ~", expression(z %~~% 0.0),")",sep="")),
                             xlim = c(4000, 4200),
                             ylim = c(35, 100))
ggsave('abb/abb_results/stellar_fit.png', stellar_p)

