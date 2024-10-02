# Produces the relevant tables and histogram plots for the estimators in Simulation 1

# --------------------- Simulation 1a ---------------------
betas = read.csv("~/Documents/Simulation 1a/Betas.csv")
vars = read.csv("~/Documents/Simulation 1a/Vars.csv")
n = 20000
Output <- data.frame(
  sq.bias = c(
    mean(betas$unw-1)^2,
    mean(betas$ours2-1)^2,
    mean(betas$ourRsq2-1)^2,
    mean(betas$wrong-1)^2
  ),
  variance = c(
    mean((betas$unw-mean(betas$unw))^2),
    mean((betas$ours2-mean(betas$ours2))^2),
    mean((betas$ourRsq2-mean(betas$ourRsq2))^2),
    mean((betas$wrong-mean(betas$wrong))^2)
  ),
  mse = c(
    mean((betas$unw-1)^2),
    mean((betas$ours2-1)^2),
    mean((betas$ourRsq2-1)^2),
    mean((betas$wrong-1)^2)
  ),
  coverage = c(
    mean( ( betas$unw-1-qnorm(0.975)*sqrt(vars$unw/n) ) * ( betas$unw-1+qnorm(0.975)*sqrt(vars$unw/n) ) < 0 ),
    mean( ( betas$ours2-1-qnorm(0.975)*sqrt(vars$ours2/n) ) * ( betas$ours2-1+qnorm(0.975)*sqrt(vars$ours2/n) ) < 0 ),
    mean( ( betas$ourRsq2-1-qnorm(0.975)*sqrt(vars$ourRsq2/n) ) * ( betas$ourRsq2-1+qnorm(0.975)*sqrt(vars$ourRsq2/n) ) < 0 ),
    mean( ( betas$wrong-1-qnorm(0.975)*sqrt(vars$wrong/n) ) * ( betas$wrong-1+qnorm(0.975)*sqrt(vars$wrong/n) ) < 0 )
  )
)
rownames(Output) <- c("Unweighted", "ROSE Random Forest", "Loceff (RSS) Random Forest", "Semiparametric Efficient Estimator")
print("Table of estimators in Simulation 1a:")
print(Output)

# Histograms for estimators in Simulation 1a (version in appendix)
cairo_pdf("apphist_1a.pdf", 8, 6)
par(mfrow=c(2,2), mar=c(2, 4, 2, 1)+0.1 )
hist(betas$unw,breaks=seq(-10,10,by=0.002),  xlim=c(0.95,1.05), ylim=c(0,45), probability=TRUE, main="Unweighted Estimator", xlab="")
axis(2, at = seq(0, 40, by = 10))
robvar = mean((betas$unw-median(betas$unw))^2,trim=0.5)*2.1899
lines(seq(0.9,1.1,by=0.0001),1/sqrt(2*pi*robvar)*exp(-(seq(0.9,1.1,by=0.0001)-1)^2/(2*robvar)),col="maroon",lwd=3)
hist(betas$wrong,breaks=seq(-10,10,by=0.002),  xlim=c(0.95,1.05), ylim=c(0,45), probability=TRUE, main="'Semiparametric Efficient' Estimator", xlab="")
axis(2, at = seq(0, 40, by = 10))
robvar = mean((betas$wrong-median(betas$wrong))^2,trim=0.5)*2.1899
lines(seq(0.9,1.1,by=0.0001),1/sqrt(2*pi*robvar)*exp(-(seq(0.9,1.1,by=0.0001)-1)^2/(2*robvar)),col="maroon",lwd=3)
hist(betas$ours3,breaks=seq(-10,10,by=0.002),  xlim=c(0.95,1.05), ylim=c(0,45), probability=TRUE, main="ROSE Estimator", xlab="")
axis(2, at = seq(0, 40, by = 10))
robvar = mean((betas$ours3-median(betas$ours3))^2,trim=0.5)*2.1899
lines(seq(0.9,1.1,by=0.0001),1/sqrt(2*pi*robvar)*exp(-(seq(0.9,1.1,by=0.0001)-1)^2/(2*robvar)),col="maroon",lwd=3)
hist(betas$ourRsq3,breaks=seq(-10,10,by=0.002),  xlim=c(0.95,1.05), ylim=c(0,45), probability=TRUE, main="Robust 'Locally Efficient' Estimator", xlab="")
axis(2, at = seq(0, 40, by = 10))
robvar = mean((betas$ourRsq3-median(betas$ourRsq3))^2,trim=0.5)*2.1899
lines(seq(0.9,1.1,by=0.0001),1/sqrt(2*pi*robvar)*exp(-(seq(0.9,1.1,by=0.0001)-1)^2/(2*robvar)),col="maroon",lwd=3)
dev.off()


# Histograms of estimators of Simulation 1 (version in main text)
cairo_pdf("maintext_hist_1a.pdf", 11.5, 3.4)
par(mfrow=c(1,3), mar=c(2, 4, 2, 1)+0.1 )
hist(betas$unw,breaks=seq(-10,10,by=0.002),  xlim=c(0.95,1.05), probability=TRUE, main="Unweighted Estimator", xlab="", ylim=c(0,45))
axis(2, at = seq(5, 45, by = 5))
robvar = mean((betas$unw-median(betas$unw))^2,trim=0.5)*2.1899
lines(seq(0.9,1.1,by=0.0001),1/sqrt(2*pi*robvar)*exp(-(seq(0.9,1.1,by=0.0001)-1)^2/(2*robvar)),col="maroon",lwd=3)
hist(betas$wrong,breaks=seq(-10,10,by=0.002),  xlim=c(0.95,1.05), probability=TRUE, main="'Semiparametric Efficient' Estimator", xlab="", ylim=c(0,45))
axis(2, at = seq(5, 45, by = 5))
robvar = mean((betas$wrong-median(betas$wrong))^2,trim=0.5)*2.1899
lines(seq(0.9,1.1,by=0.0001),1/sqrt(2*pi*robvar)*exp(-(seq(0.9,1.1,by=0.0001)-1)^2/(2*robvar)),col="maroon",lwd=3)
hist(betas$ours2,breaks=seq(-10,10,by=0.002),  xlim=c(0.95,1.05), probability=TRUE, main="Robust Semiparametric Efficient (ROSE) Estimator", xlab="", ylim=c(0,45))
axis(2, at = seq(5, 45, by = 5))
robvar = mean((betas$ours2-median(betas$ours2))^2,trim=0.5)*2.1899
lines(seq(0.9,1.1,by=0.0001),1/sqrt(2*pi*robvar)*exp(-(seq(0.9,1.1,by=0.0001)-1)^2/(2*robvar)),col="maroon",lwd=3)
dev.off()

# --------------------- Simulation 1b ---------------------

betas = read.csv("~/Documents/Simulation 1b/Betas.csv")
vars = read.csv("~/Documents/Simulation 1b/Vars.csv")
n = 20000

Output <- data.frame(
  sq.bias = c(
    mean(betas$unw-1,trim=0.01)^2,
    mean(betas$ours3-1,trim=0.01)^2,
    mean(betas$ourRsq3-1,trim=0.01)^2,
    mean(betas$wrong-1,trim=0.01)^2
  ),
  variance = c(
    mean((betas$unw-mean(betas$unw,trim=0.01))^2,trim=0.01),
    mean((betas$ours3-mean(betas$ours3,trim=0.01))^2,trim=0.01),
    mean((betas$ourRsq3-mean(betas$ourRsq3,trim=0.01))^2,trim=0.01),
    mean((betas$wrong-mean(betas$wrong,trim=0.01))^2,trim=0.01)
  ),
  mse = c(
    mean((betas$unw-1)^2,trim=0.01),
    mean((betas$ours3-1)^2,trim=0.01),
    mean((betas$ourRsq3-1)^2,trim=0.01),
    mean((betas$wrong-1)^2,trim=0.01)
  ),
  coverage = c(
    mean( ( betas$unw-1-qnorm(0.975)*sqrt(vars$unw/n) ) * ( betas$unw-1+qnorm(0.975)*sqrt(vars$unw/n) ) < 0 ),
    mean( ( betas$ours3-1-qnorm(0.975)*sqrt(vars$ours3/n) ) * ( betas$ours3-1+qnorm(0.975)*sqrt(vars$ours3/n) ) < 0 ),
    mean( ( betas$ourRsq3-1-qnorm(0.975)*sqrt(vars$ourRsq3/n) ) * ( betas$ourRsq3-1+qnorm(0.975)*sqrt(vars$ourRsq3/n) ) < 0 ),
    mean( ( betas$wrong-1-qnorm(0.975)*sqrt(vars$wrong/n) ) * ( betas$wrong-1+qnorm(0.975)*sqrt(vars$wrong/n) ) < 0 )
  )
)
rownames(Output) <- c("Unweighted", "ROSE Random Forest", "Loceff (RSS) Random Forest", "Semiparametric Efficient Estimator")
print("Table of estimators in Simulation 1b:")
print(Output)


cairo_pdf("apphist_1b.pdf", 8, 6)
par(mfrow=c(2,2), mar=c(2, 4, 2, 1)+0.1 )
hist(betas$unw,breaks=seq(-10,10,by=0.002),  xlim=c(0.95,1.05), ylim=c(0,120), probability=TRUE, main="Unweighted Estimator", xlab="")
axis(2, at = seq(0, 120, by = 20))
robvar = mean((betas$unw-median(betas$unw))^2,trim=0.5)*2.1899
lines(seq(0.9,1.1,by=0.0001),1/sqrt(2*pi*robvar)*exp(-(seq(0.9,1.1,by=0.0001)-1)^2/(2*robvar)),col="maroon",lwd=3)
hist(betas$wrong,breaks=seq(-10,10,by=0.002),  xlim=c(0.95,1.05), ylim=c(0,120), probability=TRUE, main="'Semiparametric Efficient' Estimator", xlab="")
axis(2, at = seq(0, 120, by = 20))
robvar = mean((betas$wrong-median(betas$wrong))^2,trim=0.5)*2.1899
lines(seq(0.9,1.1,by=0.0001),1/sqrt(2*pi*robvar)*exp(-(seq(0.9,1.1,by=0.0001)-1)^2/(2*robvar)),col="maroon",lwd=3)
hist(betas$ours3,breaks=seq(-10,10,by=0.002),  xlim=c(0.95,1.05), ylim=c(0,120), probability=TRUE, main="ROSE Estimator", xlab="")
axis(2, at = seq(0, 120, by = 20))
robvar = mean((betas$ours3-median(betas$ours3))^2,trim=0.5)*2.1899
lines(seq(0.9,1.1,by=0.0001),1/sqrt(2*pi*robvar)*exp(-(seq(0.9,1.1,by=0.0001)-1)^2/(2*robvar)),col="maroon",lwd=3)
hist(betas$ourRsq3,breaks=seq(-10,10,by=0.002),  xlim=c(0.95,1.05), ylim=c(0,120), probability=TRUE, main="Robust 'Locally Efficient' Estimator", xlab="")
axis(2, at = seq(0, 120, by = 20))
robvar = mean((betas$ourRsq3-median(betas$ourRsq3))^2,trim=0.5)*2.1899
lines(seq(0.9,1.1,by=0.0001),1/sqrt(2*pi*robvar)*exp(-(seq(0.9,1.1,by=0.0001)-1)^2/(2*robvar)),col="maroon",lwd=3)
dev.off()
