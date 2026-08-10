# Reproducible code for Simulation 5
# Note that one can fit a rose random forest estimator directly using [https://github.com/elliot-young/rose].
#   In this simulation, the four estimators share a number of nuisance functions,
#   and so for computational speed these estimators are estimated alongside
#   one another (as below).

# Load relevant rose random forest functions
{
  library(rpart)
  library(ranger)
  library(mgcv)
  library(glmnet)

  itemp_rf <- function(y, offset, parms, wt) {
    if (ncol(y) != 2) {
      stop("Matrix of response must be a 2 column matrix")
    }
    if (!missing(parms) && length(parms) > 0){
      warning("parameter argument ignored")
    }
    if (length(offset)) y <- y - offset
    sfun <- function(weigh, avar, ylevel, digits) {
      paste(" xisq=", format(signif(weigh[,1], digits)), " xiepsq=", format(signif(weigh[,2], digits)), ", AsymVar=", format(signif(avar, digits)), sep='')
    }
    environment(sfun) <- .GlobalEnv
    list(y = y, parms = NULL, numresp = 2, numy = 2, summary=sfun)
  }
  etemp_rf <- function(y, wt, parms) {
    sum_xisq = sum(y[,1])
    sum_xiepsq = sum(y[,2])
    avar = 10000*length(y[,1]) - sum_xisq*sum_xisq/sum_xiepsq
    list(label = cbind(sum_xisq, sum_xiepsq), deviance = avar)
  }
  stemp_rf <- function(y, wt, x, parms, continuous) {
    n <- dim(y)[1]
    if (continuous) {
      # Continuous x variable
      total_temp_xisq <- sum(y[,1])
      total_temp_xiepsq <- sum(y[,2])
      left_temp_xisq <- cumsum(y[,1])[-n]
      left_temp_xiepsq <- cumsum(y[,2])[-n]
      right_temp_xisq <- total_temp_xisq - left_temp_xisq
      right_temp_xiepsq <- total_temp_xiepsq - left_temp_xiepsq
      lavar <- left_temp_xisq * left_temp_xisq / left_temp_xiepsq
      ravar <- right_temp_xisq * right_temp_xisq / right_temp_xiepsq
      goodness <- lavar + ravar - total_temp_xisq*total_temp_xisq/total_temp_xiepsq
      list(goodness = goodness, direction = rep(1,length(x)-1))
    } else {
      stop("Not currently supported for discrete X variables")
    }
  }
  ulist_rf <- list(eval = etemp_rf, split = stemp_rf, init = itemp_rf)

}


generate_data <- function(n) {
  alpha <- 2
  Z <- rnorm(n)
  m <- function(z) 0
  f <- function(z) 0
  X <- EnvStats::rpareto(n, ((alpha-1)/alpha)*m(Z), alpha)
  Y <- X + f(Z) + rnorm(n,0,sd=abs(sqrt(X)))
  rdf <- data.frame(Z=Z, X=X, Y=Y, logX=log(X))
  return(rdf)
}

maxdepth <- maxdepth_only1 <- 1
minbucket <- 100
subsizeprop <- 0.5
numoftrees <- 500
K <- 2

all_betas <- all_vars <- list()
for (n in 250*2^(c(0:8))) {

  betas <- vars <- data.frame()
  estimators <- c("semieff", "unw", "rose_X", "rose_log")

  while (nrow(betas) < 100) {print(nrow(betas))# 4000
    beta_hat_num <- beta_hat_den <- V_num <- list()
    for (estimator in estimators) {
      beta_hat_num[[estimator]] <- numeric(K)
      beta_hat_den[[estimator]] <- numeric(K)
      V_num[[estimator]] <- numeric(K)
    }

    set.seed(nrow(betas)+1)
    rdf <- generate_data(n=n)
    K <- 2
    index <- rbinom(nrow(rdf), 1, 0.5)

    for (k in seq_len(K)) {
      rdf_train <- rdf[index==k-1,]
      rdf_test <- rdf[index!=k-1,]

      # E(Y|X,Z) = X\beta + f(Z)
      YXZ_gam_model <- mgcv::gam(Y ~ X + s(Z,bs="cr"), data = rdf_train)
      YXZ_predicted_terms <- predict(YXZ_gam_model, newdata = rdf_test, type = "terms")
      f_Z_test <- YXZ_predicted_terms[, "s(Z)"]
      model_intercept <- attr(YXZ_predicted_terms, "constant")
      f_Z_test <- f_Z_test + model_intercept

      #Y_m_f_train <- rdf_train$Y - predict(YXZ_gam_model, newdata=rdf_train)
      Y_m_f_test <- rdf_test$Y - f_Z_test


      # E(X|Z)
      X_model <- mgcv::gam(X ~ s(Z,bs="cr"), data = rdf_train)
      RX_train <- rdf_train$X - predict(X_model, newdata=rdf_train)
      RX_test <- rdf_test$X -predict(X_model, newdata=rdf_test)

      # E(log(X)|Z)
      logX_model <- mgcv::gam(logX ~ s(Z,bs="cr"), data=rdf_train)
      logRX_train <- rdf_train$logX - predict(logX_model, newdata=rdf_train)
      logRX_test <- rdf_test$logX - predict(logX_model, newdata=rdf_test)

      RX1_test <- logRX_test
      RX2_test <- RX_test

      epsq_test <- rdf_test$Y - predict(YXZ_gam_model, newdata=rdf_test)
      epsq_test <- epsq_test^2




      ep <- rdf_train$Y - predict(YXZ_gam_model, newdata=rdf_train)
      epsq <- ep^2
      xisq1 <- logRX_train*rdf_train$X #logRX_train^2
      xiepsq1 <- logRX_train^2 * epsq
      xisq2 <- RX_train*rdf_train$X #RX_train * logRX_train
      xiepsq2 <- RX_train^2 * epsq
      xi1 <- logRX_train
      xi2 <- RX_train
      DTdata <- cbind(xisq1=xisq1,xiepsq1=xiepsq1, xisq2=xisq2,xiepsq2=xiepsq2, rdf_train, xi1=xi1, xi2=xi2, epsq=epsq)

      #ep <- rdf_train$Y - predict(YXZ_gam_model, newdata=rdf_train)
      #epsq <- ep^2
      #xisq1 <- logRX_train^2
      #xiepsq1 <- logRX_train^2 * epsq
      #xisq2 <- RX_train * logRX_train
      #xiepsq2 <- RX_train^2 * epsq
      #xi1 <- logRX_train
      #xi2 <- RX_train
      #DTdata <- cbind(xisq1=xisq1,xiepsq1=xiepsq1, xisq2=xisq2,xiepsq2=xiepsq2, rdf_train, xi1=xi1, xi2=xi2, epsq=epsq)

      # Fit rose forests (J=1 and J=2 simultaneously)
      w_test_rf_1 <- matrix(0,dim(rdf_test)[1],numoftrees)
      w_test_rf_1_ONLY1_NUM <- matrix(0,dim(rdf_test)[1],numoftrees)
      w_test_rf_1_ONLY1_DEN <- matrix(0,dim(rdf_test)[1],numoftrees)
      for (fsts in seq_len(numoftrees)) {
        set.seed(fsts) # Reproducibility of rose forest
        # Trees
        selection <- sample(seq_len(dim(rdf_train)[1]), floor(subsizeprop*dim(rdf_train)[1]), replace=T)
        DTdata_sub <- DTdata[selection,]#bootstrap

        w_fit_1_ONLY1 <- rpart(paste("cbind(xisq1, xiepsq1) ~ ",paste(c("Z"),collapse="+"), collapse=""), data = DTdata_sub, method = ulist_rf, cp=0, minbucket=minbucket, maxdepth=max(c(maxdepth)))
        #w_fit_1_CORE <- rpart(paste("cbind(xisq1, xiepsq1) ~ ",paste(c("Z"),collapse="+"), collapse=""), data = DTdata_sub, method = ulist_rf, cp=0, minbucket=minbucket, maxdepth=max(c(maxdepth)))
        #w_fit_1 <- snip.rpart(w_fit_1_CORE, toss=node_labels[[maxdepth]])
        #w_fit_1_ONLY1 <- snip.rpart(w_fit_1_CORE, toss=node_labels[[maxdepth_only1]])

        pred <- predict(w_fit_1_ONLY1, rdf_test, type="matrix")
        w_test_rf_1_ONLY1_NUM[,fsts] <- pred[,1]
        w_test_rf_1_ONLY1_DEN[,fsts] <- pred[,2]
      }

      #w_test1 <- rowMeans(w_test_rf_1)
      #w_test2 <- rowMeans(w_test_rf_2)
      #beta_hat_num[["rose_log"]][k] <- sum((w_test1*RX1_test+w_test2*RX2_test)*  (Y_m_f_test)  )
      #beta_hat_den[["rose_log"]][k] <- sum((w_test1*RX1_test+w_test2*RX2_test)*  (rdf_test$X)  )
      #V_num[["rose_log"]][k] <- sum((w_test1*RX1_test+w_test2*RX2_test)^2 * epsq_test )

      w_test1_ONLY1 <- rowMeans(w_test_rf_1_ONLY1_NUM)/rowMeans(w_test_rf_1_ONLY1_DEN)
      beta_hat_num[["rose_X"]][k] <- sum((w_test1_ONLY1*RX1_test)*  (Y_m_f_test)  )
      beta_hat_den[["rose_X"]][k] <- sum((w_test1_ONLY1*RX1_test)*  (rdf_test$X)  )
      V_num[["rose_X"]][k] <- sum((w_test1_ONLY1*RX1_test)^2 * epsq_test )

      # Unweighted estimator
      beta_hat_num[["unw"]][k] <- sum(Y_m_f_test*RX_test)
      beta_hat_den[["unw"]][k] <- sum(RX_test*rdf_test$X)
      V_num[["unw"]][k] <- sum(RX_test^2*epsq_test)


      # FIT weighted RF with M(X)=X
      w_test_rf_1 <- matrix(0,dim(rdf_test)[1],numoftrees)
      w_test_rf_1_ONLY1_NUM <- matrix(0,dim(rdf_test)[1],numoftrees)
      w_test_rf_1_ONLY1_DEN <- matrix(0,dim(rdf_test)[1],numoftrees)
      for (fsts in seq_len(numoftrees)) {
        set.seed(fsts) # Reproducibility of rose forest
        # Trees
        selection <- sample(seq_len(dim(rdf_train)[1]), floor(subsizeprop*dim(rdf_train)[1]), replace=T)
        DTdata_sub <- DTdata[selection,]#bootstrap

        w_fit_1_ONLY1 <- rpart(paste("cbind(xisq2, xiepsq2) ~ ",paste(c("Z"),collapse="+"), collapse=""), data = DTdata_sub, method = ulist_rf, cp=0, minbucket=minbucket, maxdepth=max(c(maxdepth)))
        pred <- predict(w_fit_1_ONLY1, rdf_test, type="matrix")
        w_test_rf_1_ONLY1_NUM[,fsts] <- pred[,1]
        w_test_rf_1_ONLY1_DEN[,fsts] <- pred[,2]
      }
      w_test1_ONLY1 <- rowMeans(w_test_rf_1_ONLY1_NUM)/rowMeans(w_test_rf_1_ONLY1_DEN)
      beta_hat_num[["rose_log"]][k] <- sum((w_test1_ONLY1*RX2_test)*  (Y_m_f_test)  )
      beta_hat_den[["rose_log"]][k] <- sum((w_test1_ONLY1*RX2_test)*  (rdf_test$X)  )
      V_num[["rose_log"]][k] <- sum((w_test1_ONLY1*RX2_test)^2 * epsq_test )


      # Semiparametric efficient estimator
      sigmasq_model <- mgcv::gam(epsq ~ s(X,bs="cr")+s(Z,bs="cr"), data = cbind(epsq,rdf_train))
      W_train <- 1/predict(sigmasq_model, newdata=rdf_train)
      W_test <- 1/predict(sigmasq_model, newdata=rdf_test)
      W_train <- pmax(W_train, 0)
      W_test <- pmax(W_test, 0)

      h_mod <- mgcv::gam(X ~ s(Z,bs="cr"), data = rdf_train, weights = W_train)
      h_test <- predict(h_mod, newdata=rdf_test)
      RtildeX1 <- rdf_test$X - h_test

      beta_hat_num[["semieff"]][k] <- sum(W_test*Y_m_f_test*RtildeX1)
      beta_hat_den[["semieff"]][k] <- sum(W_test*rdf_test$X*RtildeX1)
      V_num[["semieff"]][k] <- sum(W_test^2*epsq_test*RtildeX1^2)

    }

    beta_single <- Map(function(x1,x2) sum(x1)/sum(x2), beta_hat_num, beta_hat_den)
    betas_single <- as.data.frame(t(unlist(beta_single)))
    var_single <- Map(function(x1,x2) n*sum(x1)/(sum(x2))^2, V_num, beta_hat_den)
    vars_single <- as.data.frame(t(unlist(var_single)))
    betas <- rbind(betas, betas_single)
    vars <- rbind(vars, vars_single)

  }
  all_betas[[paste0(n)]] <- betas
  all_vars[[paste0(n)]] <- vars

}




sizes <- 250*2^(c(0:8))
methods <- c("semieff", "unw", "rose_X", "rose_log")

mse_table <- sapply(sizes, function(n) {
  sapply(methods, function(m) {
    mean((all_betas[[as.character(n)]][[m]] - 1)^2)
  })
})

mse_table <- t(mse_table)
colnames(mse_table) <- methods
rownames(mse_table) <- sizes

print(mse_table)

