# Code for producing plots for Simulation 2

library(ggplot2)
library(RColorBrewer)
coul <- brewer.pal(4, "Dark2")
coul_rgb <- col2rgb(coul)
library(gridExtra)
library(ggpubr)
get_legend<-function(myggplot){
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}

ours_rel_eff <- ourRsq_rel_eff <- numeric()
trim_val <- 0.0
T <- 1000
# Load data
BETAS_ALL_ns = readRDS("~/Documents/Simulation 2/BETAS.RData")
CVS_ALL_ns = readRDS("~/Documents/Simulation 2/CVS.RData")
VARS_ALL_ns = readRDS("~/Documents/Simulation 2/VAR.RData")
ns <- 2^(1:8)*10^2
for (n in ns) {
    BETAS_ALL = BETAS_ALL_ns[[paste0("n_",n)]]
    CVS_ALL = CVS_ALL_ns[[paste0("n_",n)]]
    VARS_ALL = VARS_ALL_ns[[paste0("n_",n)]]
    {
        betas_unw <- BETAS_ALL[[1]]$UNW
        betas_ora <- BETAS_ALL[[1]]$ORACLE
        betas_ours <- rep(0,T)
        betas_ourRsq <- rep(0,T)
        MD = length(BETAS_ALL)
        for (tt in 1:T) {
            cvs_here = rep(0,MD+1)
            cvs_here[1] = CVS_ALL[[1]]$ourRsq0[tt]
            for (pp in 1:MD) {
                cvs_here[1+pp] = CVS_ALL[[pp]]$ourRsq1[tt]
            }
            placement = which.min(cvs_here)
            if (placement==1) betas_ourRsq[tt] = betas_unw[tt]
            if (placement!=1) betas_ourRsq[tt] = BETAS_ALL[[placement-1]]$ourRsq1[tt]
        }
        for (tt in 1:T) {
            cvs_here = rep(0,MD+1)
            cvs_here[1] = VARS_ALL[[1]]$ours3[tt]
            for (pp in 1:MD) {
                cvs_here[1+pp] = VARS_ALL[[pp]]$ours3[tt]
            }
            placement = which.min(cvs_here)
            if (placement==1) betas_ours[tt] = betas_unw[tt]
            if (placement!=1) betas_ours[tt] = BETAS_ALL[[placement-1]]$ours3[tt]
        }
        
        mse_unw = mean((betas_unw-1)^2,trim=trim_val)
        mse_ora = mean((betas_ora-1)^2,trim=trim_val)
        mse_ours = mean((betas_ours-1)^2,trim=trim_val)
        mse_ourRsq = mean((betas_ourRsq-1)^2,trim=trim_val)
        
    }
    ours_rel_eff <- append(ours_rel_eff, ((mse_unw-mse_ours)/(mse_unw-mse_ora)*100))
    ourRsq_rel_eff <- append(ourRsq_rel_eff, ((mse_unw-mse_ourRsq)/(mse_unw-mse_ora)*100))
}

# Sole purpose a quick workaround for a bespoke legend
{
    MAXDEPTH = 1:16
    y_ours = y_rsq = rep(1,length(MAXDEPTH))
    MSE.DF <- data.frame(Estimators=c(rep("ROSE random forests",length(MAXDEPTH)),rep("CART random forests",length(MAXDEPTH)),rep("Unweighted",length(MAXDEPTH))), SampleSize=rep(MAXDEPTH,3), MSE=c(y_ours,y_rsq,rep(1,length(y_rsq))))
    legend_metadata = rbind(MSE.DF,MSE.DF[1:16,]); legend_metadata[49:64,1] = rep("Oracle",16)
    legend_metadata[65:80,] = legend_metadata[49:64,]
    legend_plot = ggplot(legend_metadata, aes(x=SampleSize, y=MSE, group=Estimators, color=Estimators)) +
        geom_line(position=position_dodge(width = 0.001), size=0.8) +
        geom_point(position=position_dodge(width = 0.001), size=1.4) +
        theme_bw() +
        labs(x = expression(lambda), y = "Ratio of MSE over Oracle") +
        scale_color_manual(
            values=c(coul[2],coul[1],"black",coul[4]),
            labels=c("ROSE random forests","CART random forests","Unweighted","Oracle"),
            name="Estimator",
            guide = guide_legend(override.aes = list(
                shape = c(16, 16, NA, NA),
                linetype = c("solid","solid","solid","solid")
            ))
        )
    legend = get_legend(legend_plot)
}
# Plots
RellEff.DF <- data.frame(Estimators=c(rep("ROSE random forests",length(ns)),rep("CART random forests",length(ns))), NumberOfObservations=rep(ns,2), RelEff=c(ours_rel_eff,ourRsq_rel_eff))
plot_RELEFF <- ggplot(RellEff.DF, aes(x=NumberOfObservations, y=RelEff, group=Estimators, color=Estimators, linetype=Estimators)) +
    geom_line(position=position_dodge(width = 0.0), size=0.8) + #0.05
    geom_point(position=position_dodge(width = 0.0), size=1.4) + #0.05
    theme_minimal() +
    geom_abline(intercept = 100, slope = 0, col=coul[4], size=0.8) +
    geom_abline(intercept = 0, slope = 0, col="black", linetype="solid", size=0.8) +
    labs(x = "Sample size", y = "Relative efficiency") +
    scale_color_manual(values=c(coul[1], coul[2])) +
    scale_linetype_manual(values=c(1,1)) +
    scale_x_continuous(breaks = ns, minor_breaks = NULL, trans='log10') +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("Simulation 2: Comparing splitting rules") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylim(0, 100)
plot_RELEFF



#cairo_pdf("~/Desktop/ROSE Paper/sim2.pdf", 5, 4)
#cairo_pdf("~/Downloads/sim2.pdf", 6.7, 4)
cairo_pdf("~/Downloads/sim2-ROSE.pdf", 6.4, 4)
ggarrange(plot_RELEFF, ncol=1, nrow=1, common.legend = TRUE, legend="right", legend.grob = legend)
dev.off()


# Plots for Simulation 2 in appendix

# For the case: n = 6400
BETAS_ALL = BETAS_ALL_ns$n_6400
CVS_ALL = CVS_ALL_ns$n_6400
VARS_ALL = VARS_ALL_ns$n_6400

trim_val=0.0
{library(ggplot2)
  library(RColorBrewer)
  coul <- brewer.pal(4, "Dark2")
  coul_rgb <- col2rgb(coul)
  library(gridExtra)
  library(ggpubr)
  get_legend<-function(myggplot){
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }}
mse_ours_eachdepth = mse_ourRsq_eachdepth = rep(0,16)
mse_ours_eachdepth[1] = mean((BETAS_ALL[[1]]$UNW-1)^2,trim=trim_val)
mse_ourRsq_eachdepth[1] = mean((BETAS_ALL[[1]]$UNW-1)^2,trim=trim_val)
for (maxdepth in 1:15) {
  mse_ours_eachdepth[maxdepth+1] = mean((BETAS_ALL[[maxdepth]]$ours3-1)^2,trim=trim_val)
  mse_ourRsq_eachdepth[maxdepth+1] = mean((BETAS_ALL[[maxdepth]]$ourRsq1-1)^2,trim=trim_val)
}
mse_ours_eachdepth = mse_ours_eachdepth/mean((BETAS_ALL[[1]]$UNW-1)^2,trim=trim_val)
mse_ourRsq_eachdepth = mse_ourRsq_eachdepth/mean((BETAS_ALL[[1]]$UNW-1)^2,trim=trim_val)
{
  MAXDEPTH = 1:16
  y_ours = y_rsq = rep(1,length(MAXDEPTH))
  MSE.DF <- data.frame(Estimators=c(rep("ROSE random forests",length(MAXDEPTH)),rep("CART random forests",length(MAXDEPTH)),rep("Unweighted",length(MAXDEPTH))), SampleSize=rep(MAXDEPTH,3), MSE=c(y_ours,y_rsq,rep(1,length(y_rsq))))
  legend_metadata = rbind(MSE.DF,MSE.DF[1:16,]); legend_metadata[49:64,1] = rep("Oracle",16)
  legend_metadata[65:80,] = legend_metadata[49:64,]
  legend_metadata[65:80,1] = rep("Cross-validated over depth",16)
  legend_plot = ggplot(legend_metadata, aes(x=SampleSize, y=MSE, group=Estimators, color=Estimators)) +
    geom_line(position=position_dodge(width = 0.001), size=0.8) +
    geom_point(position=position_dodge(width = 0.001), size=1.4) +
    theme_bw() +
    labs(x = expression(lambda), y = "Ratio of MSE over Oracle") +
    scale_color_manual(
      values=c(coul[2],coul[1],"black",coul[4],"grey"),
      labels=c("ROSE random forests","CART random forests","Unweighted","Oracle","Cross-validated over depth"),
      name="Estimator",
      guide = guide_legend(override.aes = list(
        shape = c(16, 16, NA, NA, NA),
        linetype = c("solid","solid","solid","solid","dashed")
      ))
    )
  legend = get_legend(legend_plot)
}
# ONLY PURPOSE TO FIND mse_ora, mse_unw, mse_ours and mse_ourRsq
{
  betas_unw <- BETAS_ALL[[1]]$UNW
  betas_ora <- BETAS_ALL[[1]]$ORACLE
  betas_ours <- rep(0,T)
  betas_ourRsq <- rep(0,T)
  MD = length(BETAS_ALL)
  for (tt in 1:T) {
    cvs_here = rep(0,MD+1)
    cvs_here[1] = CVS_ALL[[1]]$ourRsq0[tt]
    for (pp in 1:MD) {
      cvs_here[1+pp] = CVS_ALL[[pp]]$ourRsq1[tt]
    }
    placement = which.min(cvs_here)
    if (placement==1) betas_ourRsq[tt] = betas_unw[tt]
    if (placement!=1) betas_ourRsq[tt] = BETAS_ALL[[placement-1]]$ourRsq1[tt]
  }
  for (tt in 1:T) {
    cvs_here = rep(0,MD+1)
    cvs_here[1] = VARS_ALL[[1]]$ours3[tt]
    for (pp in 1:MD) {
      cvs_here[1+pp] = VARS_ALL[[pp]]$ours3[tt]
    }
    placement = which.min(cvs_here)
    if (placement==1) betas_ours[tt] = betas_unw[tt]
    if (placement!=1) betas_ours[tt] = BETAS_ALL[[placement-1]]$ours3[tt]
  }

  mse_unw = mean((betas_unw-1)^2,trim=trim_val)
  mse_ora = mean((betas_ora-1)^2,trim=trim_val)
  mse_ours = mean((betas_ours-1)^2,trim=trim_val)
  mse_ourRsq = mean((betas_ourRsq-1)^2,trim=trim_val)

}
MSE.DF <- data.frame(Estimators=c(rep("ROSE random forests",16),rep("CART random forests",16)), SampleSize=rep(0:15,2), MSE=c(mse_ours_eachdepth,mse_ourRsq_eachdepth))
plot_SIM2 <- ggplot(MSE.DF, aes(x=SampleSize, y=MSE, group=Estimators, color=Estimators, linetype=Estimators)) +
  geom_line(position=position_dodge(width = 0.0), size=0.8) + #0.05
  geom_point(position=position_dodge(width = 0.0), , size=1.4) + #0.05
  #theme_bw() +
  theme_minimal() +
  geom_abline(intercept = mse_ora/mse_unw, slope = 0, col=coul[4], size=0.8) +
  geom_abline(intercept = mse_ours/mse_unw, slope = 0, col=coul[2], linetype="dashed", size=0.8) +
  geom_abline(intercept = mse_ourRsq/mse_unw, slope = 0, col=coul[1], linetype="dashed", size=0.8) +
  geom_abline(intercept = 1, slope = 0, col="black", linetype="solid", size=0.8) +
  labs(x = "Tree depth (hyperparameter)", y = "Mean Squared Error") +
  scale_color_manual(values=c(coul[1], coul[2])) +
  scale_linetype_manual(values=c(1,1)) +
  scale_x_continuous(breaks = c(0,MAXDEPTH), minor_breaks = NULL) +
  #scale_y_continuous(trans='log10') +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Simulation 2: Comparing splitting rules") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(0.3, 1.0)
plot_SIM2

cairo_pdf("~/Downloads/sim2_6400_maxdepth_tuning.pdf", 6.5, 4)
ggarrange(plot_SIM2, ncol=1, nrow=1, common.legend = TRUE, legend="right", legend.grob = legend)
dev.off()
