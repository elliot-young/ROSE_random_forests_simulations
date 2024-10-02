# Produces the relevant plots for the estimators in Simulation 3

# Load data from simulations in simulation_3.R
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
unw_col = 'black'
rose1_col = brewer.pal(n = 8, name = "Dark2")[2]
rose2_col = brewer.pal(n = 8, name = "Dark2")[6]
eff_col = brewer.pal(n = 8, name = "Set1")[2]

SS = list()
for (n in 10000*c(1,2,4,8)) SS[[paste0(n)]] = n
sqBias_div_Var <- Coverage <- data.frame(unw=rep(0,4), semieff=rep(0,4), roseJ2=rep(0,4), roseJ1=rep(0,4))
sqBias_div_Var_var <- data.frame(unw=rep(0,4), semieff=rep(0,4), roseJ2=rep(0,4), roseJ1=rep(0,4))

sqBias_div_Var$unw <- unlist(lapply(BETAS, function(x) mean(x$unw-1)^2/mean((x$unw-1)^2)*100))
sqBias_div_Var$semieff <- unlist(lapply(BETAS, function(x) mean(x$semieff-1)^2/mean((x$semieff-1)^2)*100))
sqBias_div_Var$roseJ2 <- unlist(lapply(BETAS, function(x) mean(x$roseJ2-1)^2/mean((x$roseJ2-1)^2)*100))
sqBias_div_Var$roseJ1 <- unlist(lapply(BETAS, function(x) mean(x$roseJ1-1)^2/mean((x$roseJ1-1)^2)*100))
sqBias_div_Var_var$unw <- unlist(lapply(BETAS, function(x) var(x$unw)/4000*(200*mean(x$unw-1)/mean((x$unw-1)^2))^2 )) # Delta method
sqBias_div_Var_var$semieff <- unlist(lapply(BETAS, function(x) var(x$semieff)/4000*(200*mean(x$semieff-1)/mean((x$semieff-1)^2))^2 ))
sqBias_div_Var_var$roseJ2 <- unlist(lapply(BETAS, function(x) var(x$roseJ2)/4000*(200*mean(x$roseJ2-1)/mean((x$roseJ2-1)^2))^2 ))
sqBias_div_Var_var$roseJ1 <- unlist(lapply(BETAS, function(x) var(x$roseJ1)/4000*(200*mean(x$roseJ1-1)/mean((x$roseJ1-1)^2))^2 ))

Coverage$unw <- unlist(Map(function(x,y,z) 100*mean(as.numeric((x$unw-1-qnorm(0.975)*sqrt(y$unw/z))*(x$unw-1+qnorm(0.975)*sqrt(y$unw/z))<0)), BETAS, VARS, SS))
Coverage$semieff <- unlist(Map(function(x,y,z) 100*mean(as.numeric((x$semieff-1-qnorm(0.975)*sqrt(y$semieff/z))*(x$semieff-1+qnorm(0.975)*sqrt(y$semieff/z))<0)), BETAS, VARS, SS))
Coverage$roseJ2 <- unlist(Map(function(x,y,z) 100*mean(as.numeric((x$roseJ2-1-qnorm(0.975)*sqrt(y$roseJ2/z))*(x$roseJ2-1+qnorm(0.975)*sqrt(y$roseJ2/z))<0)), BETAS, VARS, SS))
Coverage$roseJ1 <- unlist(Map(function(x,y,z) 100*mean(as.numeric((x$roseJ1-1-qnorm(0.975)*sqrt(y$roseJ1/z))*(x$roseJ1-1+qnorm(0.975)*sqrt(y$roseJ1/z))<0)), BETAS, VARS, SS))
Coverage_var <- 100^2*(Coverage/100)*(1-Coverage/100)/4000

plot1data <- data.frame(
  sample_size = rep(10000*2^(0:3), times = 4),
  mse = c(sqBias_div_Var$unw,
          sqBias_div_Var$roseJ1,
          sqBias_div_Var$roseJ2,
          sqBias_div_Var$semieff),
  sd = c(qnorm(0.975)*sqrt(sqBias_div_Var_var$unw),
         qnorm(0.975)*sqrt(sqBias_div_Var_var$roseJ1),
         qnorm(0.975)*sqrt(sqBias_div_Var_var$roseJ2),
         qnorm(0.975)*sqrt(sqBias_div_Var_var$semieff)),
  estimator = rep(c("Unweighted", "Rose random forest (J=1)", "Rose random forest (J=2)", "Semiparametric efficient"), each = 4)
)
plot1data$estimator <- factor(plot1data$estimator, levels = c("Unweighted", "Rose random forest (J=1)", "Rose random forest (J=2)", "Semiparametric efficient"))

plot1 <- ggplot(plot1data, aes(x = sample_size, y = mse, color = estimator)) +
  geom_line(position = position_dodge(width = 0.04), size = 1) +
  geom_point(position = position_dodge(width = 0.04), size = 2) +
  geom_errorbar(aes(ymin = mse - sd, ymax = mse + sd), position = position_dodge(width = 0.04), width = 0.1) +
  labs(title = "Contribution of squared bias to MSE",
       x = "Sample size",
       y = "Contribution of squared bias to MSE",
       color = "Estimator") +
  theme_minimal() +
  scale_x_log10(breaks = 10000*2^(0:3), labels = c("10 000","20 000","40 000","80 000")) +
  scale_y_continuous(breaks = c(0,10,20,30), labels = c("0%","10%","20%","30%")) +
  expand_limits(y = c(0,30)) +
  scale_color_manual(values = c("Unweighted" = unw_col,
                                "Rose random forest (J=1)" = rose1_col,
                                "Rose random forest (J=2)" = rose2_col,
                                "Semiparametric efficient" = eff_col))


plot2data <- data.frame(
  sample_size = rep(10000*2^(0:3), times = 4),
  cov = c(Coverage$unw,
          Coverage$roseJ1,
          Coverage$roseJ2,
          Coverage$semieff),
  sd = c(qnorm(0.975)*sqrt(Coverage_var$unw),
         qnorm(0.975)*sqrt(Coverage_var$roseJ1),
         qnorm(0.975)*sqrt(Coverage_var$roseJ2),
         qnorm(0.975)*sqrt(Coverage_var$semieff)),
  estimator = rep(c("Unweighted", "Rose random forest (J=1)", "Rose random forest (J=2)", "Semiparametric efficient"), each = 4)
)
plot2data$estimator <- factor(plot1data$estimator, levels = c("Unweighted", "Rose random forest (J=1)", "Rose random forest (J=2)", "Semiparametric efficient"))

plot2 <- ggplot(plot2data, aes(x = sample_size, y = cov, color = estimator)) +
  geom_line(position = position_dodge(width = 0.04), size = 1) +
  geom_point(position = position_dodge(width = 0.04), size = 2) +
  geom_errorbar(aes(ymin = cov - sd, ymax = cov + sd), position = position_dodge(width = 0.04), width = 0.1) +
  labs(title = "Coverage (of 95% confidence intervals)",
       x = "Sample size",
       y = "Coverage (of 95% confidence intervals)",
       color = "Estimator") +
  theme_minimal() +
  scale_x_log10(breaks = 10000*2^(0:3), labels = c("10 000","20 000","40 000","80 000")) +
  scale_y_continuous(breaks = c(90,92,94,96), labels = c("90%","92%","94%","96%")) +
  expand_limits(y = c(90,96)) +
  scale_color_manual(values = c("Unweighted" = unw_col,
                                "Rose random forest (J=1)" = rose1_col,
                                "Rose random forest (J=2)" = rose2_col,
                                "Semiparametric efficient" = eff_col)) +
  geom_abline(intercept=95, slope=0, linetype = "dotted", color = "grey", size=1)


get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
cairo_pdf("~/Downloads/sim3plot.pdf",11,4)
grid.arrange(plot1+theme(legend.position="none"), plot2+theme(legend.position="none"), get_legend(plot1), ncol=3, widths=c(3, 3, 1.2))
dev.off()

