# Reproducible code for Simulation 3
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
expit = function(x) exp(x)/(1+exp(x))
generate_data <- function() {
  Z <- matrix(rnorm(n*3),n,3)
  X <- rep(0,n)
  p <- 0.1+0.8*expit(1*Z[,1])
  B <- rbinom(n,1,p)
  for (ii in 1:n) {
    mu <- expit(Z[ii,1])+expit(Z[ii,2])+expit(Z[ii,3])
    sigma_sq <- 0.1
    if (B[ii]==0) {
      if (sigma_sq!=0) X[ii] <- rgamma(1,shape=mu^2/sigma_sq, scale=sigma_sq/mu)
      if (sigma_sq==0) X[ii] <- mu
    }
  }
  Y <- rep(0,n)
  for (ii in 1:n) {
    mu <- X[ii] + 1*(expit(Z[ii,1])+expit(Z[ii,2])+expit(Z[ii,3]))
    if (X[ii]==0) sigma_sq <- 0.9
    if (X[ii]!=0 & X[ii]<1.5) sigma_sq <- 0.4
    if (X[ii]>=1.5) sigma_sq <- 0.1
    if (sigma_sq!=0) Y[ii] <- rgamma(1,shape=mu^2/sigma_sq, scale=sigma_sq/mu)
    if (sigma_sq==0) Y[ii] <- mu
  }
  Xzero <- as.numeric(X==0)
  rdf <- data.frame(Z=Z, X=X, Xzero=Xzero, Y=Y)
  return(rdf)
}
maxdepth <- maxdepth_only1 <- 5
minbucket <- 100
subsizeprop <- 1
numoftrees <- 500
K <- 2
node_labels <- list()
for (maxdepth in 1:9) node_labels[[maxdepth]] <- (2^(maxdepth)) : (2^(maxdepth+1)-1)

# Hyperparameters for each sample size (see file Simulation 3 hyperparameters cross validation)
n_hyperparameters <- list(
  `10000` = list(
    X_model_MNS = 100,
    Y_model_MNS = 30,
    sigsq_model_MNS = 1000,
    h_model_MNS = 20,
    X_model_SP = 0.20,
    Y_model_SP = 0.10,
    sigsq_model_SP = 0.80,
    h_model_SP = 0.05,
    Xzero_model_MNS = 100,
    Xzero_model_SP = 0.20
  ),
  `20000` = list(
    X_model_MNS = 100,
    Y_model_MNS = 50,
    sigsq_model_MNS = 1000,
    h_model_MNS = 30,
    X_model_SP = 0.10,
    Y_model_SP = 0.10,
    sigsq_model_SP = 0.50,
    h_model_SP = 0.05,
    Xzero_model_MNS = 500,
    Xzero_model_SP = 0.50
  ),
  `40000` = list(
    X_model_MNS = 200,
    Y_model_MNS = 30,
    sigsq_model_MNS = 2000,
    h_model_MNS = 50,
    X_model_SP = 0.20,
    Y_model_SP = 0.05,
    sigsq_model_SP = 0.50,
    h_model_SP = 0.05,
    Xzero_model_MNS = 200,
    Xzero_model_SP = 0.20
  ),
  `80000` = list(
    X_model_MNS = 200,
    Y_model_MNS = 50,
    sigsq_model_MNS = 2000,
    h_model_MNS = 50,
    X_model_SP = 0.10,
    Y_model_SP = 0.05,
    sigsq_model_SP = 0.50,
    h_model_SP = 0.05,
    Xzero_model_MNS = 500,
    Xzero_model_SP = 0.20
  ),
  `160000` = list(
    X_model_MNS = 200,
    Y_model_MNS = 50,
    sigsq_model_MNS = 2000,
    h_model_MNS = 100,
    X_model_SP = 0.10,
    Y_model_SP = 0.05,
    sigsq_model_SP = 0.20,
    h_model_SP = 0.05,
    Xzero_model_MNS = 500,
    Xzero_model_SP = 0.20
  )
)

all_betas <- all_vars <- list()
for (n in 2^(0:4)*10^4) {

  betas <- vars <- data.frame()
  estimators <- c("semieff", "unw", "roseJ1", "roseJ2")

  while (nrow(betas) < 4000) {
    beta_hat_num <- beta_hat_den <- V_num <- list()
    for (estimator in estimators) {
      beta_hat_num[[estimator]] <- numeric(K)
      beta_hat_den[[estimator]] <- numeric(K)
      V_num[[estimator]] <- numeric(K)
    }

    set.seed(nrow(betas)+1)
    rdf <- generate_data()
    K <- 2
    index <- rbinom(nrow(rdf), 1, 0.5)

    for (k in seq_len(K)) {
      rdf_train <- rdf[index==k-1,]
      rdf_test <- rdf[index!=k-1,]

      # P(X=0|Z)
      X2_model <- grf::probability_forest(X=rdf_train[,1:3], Y=as.factor(rdf_train[,"Xzero"]), num.trees=500, sample.fraction = n_hyperparameters[[paste0(n)]]$Xzero_model_SP, honesty = FALSE, min.node.size = n_hyperparameters[[paste0(n)]]$Xzero_model_MNS)
      RX2_train <- rdf_train$Xzero - predict(X2_model, rdf_train[,1:3])$predictions[,2]
      RX2_test <- rdf_test$Xzero - predict(X2_model, rdf_test[,1:3])$predictions[,2]
      # E(Y|Z)
      Y_model <- ranger(Y ~ Z.1+Z.2+Z.3, data=rdf_train, min.node.size = n_hyperparameters[[paste0(n)]]$Y_model_MNS, sample.fraction = n_hyperparameters[[paste0(n)]]$Y_model_SP, num.trees <- 500)
      RY_train <- rdf_train$Y - predict(Y_model, rdf_train)$predictions
      RY_test <- rdf_test$Y - predict(Y_model, rdf_test)$predictions
      # E(X|Z)
      X_model <- ranger(X ~ Z.1+Z.2+Z.3, data=rdf_train, min.node.size = n_hyperparameters[[paste0(n)]]$X_model_MNS, sample.fraction = n_hyperparameters[[paste0(n)]]$X_model_SP, num.trees <- 500)
      RX_train <- rdf_train$X - predict(X_model, rdf_train)$predictions
      RX_test <- rdf_test$X - predict(X_model, rdf_test)$predictions

      RX1_train <- RX_train
      RX1_test <- RX_test
      beta_train <- coefficients(lm(RY_train ~ RX1_train - 1))

      ep <- RY_train-beta_train[1]*RX1_train
      epsq <- ep^2
      xisq1 <- RX1_train^2
      xiepsq1 <- RX1_train^2 * epsq
      xisq2 <- RX2_train * RX1_train
      xiepsq2 <- RX2_train^2 * epsq
      xi1 <- RX1_train
      xi2 <- RX2_train
      epsq_test <- (RY_test-beta_train[1]*RX1_test)^2

      DTdata <- cbind(xisq1=xisq1,xiepsq1=xiepsq1, xisq2=xisq2,xiepsq2=xiepsq2, rdf_train, xi1=xi1, xi2=xi2, epsq=epsq)

      # Semiparametric efficient estimator
      sigmasq_model <- ranger(epsq ~ X+Z.1+Z.2+Z.3, data = cbind(epsq,rdf_train), min.node.size = n_hyperparameters[[paste0(n)]]$sigsq_model_MNS, sample.fraction = n_hyperparameters[[paste0(n)]]$sigsq_model_SP, num.trees = 500)
      W_train <- 1/predict(sigmasq_model, rdf_train)$predictions
      W_test <- 1/predict(sigmasq_model, rdf_test)$predictions
      h_mod <- ranger(X ~ Z.1+Z.2+Z.3, data = rdf_train, min.node.size = n_hyperparameters[[paste0(n)]]$h_model_MNS, sample.fraction = n_hyperparameters[[paste0(n)]]$h_model_SP, num.trees = 500, case.weights = W_train)
      h_test <- predict(h_mod, rdf_test)$predictions
      RtildeX1 <- rdf_test$X - h_test

      beta_hat_num[["semieff"]][k] <- sum(W_test*RY_test*RtildeX1)
      beta_hat_den[["semieff"]][k] <- sum(W_test*RX1_test*RtildeX1)
      V_num[["semieff"]][k] <- sum(W_test^2*(RY_test-beta_train*RX1_test)^2*RtildeX1^2)

      # Fit rose forests (J=1 and J=2 simultaneously)
      w_test_rf_1 <- matrix(0,dim(rdf_test)[1],numoftrees)
      w_test_rf_2 <- matrix(0,dim(rdf_test)[1],numoftrees)
      w_test_rf_1_ONLY1_NUM <- matrix(0,dim(rdf_test)[1],numoftrees)
      w_test_rf_1_ONLY1_DEN <- matrix(0,dim(rdf_test)[1],numoftrees)
      for (fsts in seq_len(numoftrees)) {
        set.seed(fsts) # Reproducibility of rose forest
        # Trees
        selection <- sample(seq_len(dim(rdf_train)[1]), floor(subsizeprop*dim(rdf_train)[1]), replace=T)
        DTdata_sub <- DTdata[selection,]#bootstrap
        w_fit_1_CORE <- rpart(paste("cbind(xisq1, xiepsq1) ~ ",paste(c("Z.1","Z.2","Z.3"),collapse="+"), collapse=""), data = DTdata_sub, method = ulist_rf, cp=0, minbucket=minbucket, maxdepth=max(c(maxdepth,maxdepth_only1)))
        w_fit_1 <- snip.rpart(w_fit_1_CORE, toss=node_labels[[maxdepth]])
        w_fit_1_ONLY1 <- snip.rpart(w_fit_1_CORE, toss=node_labels[[maxdepth_only1]])
        w_fit_2 <- rpart(paste("cbind(xisq2, xiepsq2) ~ ",paste(c("Z.1","Z.2","Z.3"),collapse="+"), collapse=""), data = DTdata_sub, method = ulist_rf, cp=0, minbucket=minbucket, maxdepth=maxdepth)

        leaves1 <- which(w_fit_1$frame$var=="<leaf>")
        num_leaves1 <- length(leaves1)
        leaves2 <- which(w_fit_2$frame$var=="<leaf>")
        num_leaves2 <- length(leaves2)

        # Calculate F11, F22
        F11 <- diag(w_fit_1$frame$yval2[leaves1,2])
        F22 <- diag(w_fit_2$frame$yval2[leaves2,2])
        # Calculate F12
        w_fit_1_where <- w_fit_1$where
        w_fit_2_where <- w_fit_2$where
        w_fit_1_where_1ton_indexed <- match(w_fit_1_where, leaves1)
        w_fit_2_where_1ton_indexed <- match(w_fit_2_where, leaves2)
        xi1xi2epsq <- DTdata_sub$xi1 * DTdata_sub$xi2 * DTdata_sub$epsq
        F12 <- matrix(0,num_leaves1,num_leaves2)
        for (i in 1:length(xi1)) {
          F12[w_fit_1_where_1ton_indexed[i], w_fit_2_where_1ton_indexed[i]] <-
            F12[w_fit_1_where_1ton_indexed[i], w_fit_2_where_1ton_indexed[i]] + xi1xi2epsq[i]
        }

        Fall <- rbind(cbind(F11, F12),cbind(t(F12), F22))

        phi1 <- w_fit_1$frame$yval2[leaves1,1]
        phi2 <- w_fit_2$frame$yval2[leaves2,1]
        phiall <- c(phi1, phi2)

        opweights <- solve(Fall, phiall)
        opweights_1 <- opweights[1:num_leaves1]
        opweights_2 <- opweights[(1+num_leaves1):(num_leaves1+num_leaves2)]

        # Overwrite original decision trees
        fresh_w_fit_1 <- w_fit_1
        fresh_w_fit_2 <- w_fit_2
        fresh_w_fit_1$frame$yval[leaves1] <- opweights_1
        fresh_w_fit_2$frame$yval[leaves2] <- opweights_2

        w_fit_test_1 <- predict(fresh_w_fit_1, newdata=rdf_test, type="vector")
        w_fit_test_2 <- predict(fresh_w_fit_2, newdata=rdf_test, type="vector")

        w_test_rf_1[,fsts] <- w_fit_test_1
        w_test_rf_2[,fsts] <- w_fit_test_2

        w_test_rf_1_ONLY1_NUM[,fsts] <- predict(w_fit_1_ONLY1, rdf_test, type="matrix")[,1]
        w_test_rf_1_ONLY1_DEN[,fsts] <- predict(w_fit_1_ONLY1, rdf_test, type="matrix")[,2]
      }

      w_test1 <- rowMeans(w_test_rf_1)
      w_test2 <- rowMeans(w_test_rf_2)
      beta_hat_num[["roseJ2"]][k] <- sum((w_test1*RX1_test+w_test2*RX2_test)*  (RY_test)  )
      beta_hat_den[["roseJ2"]][k] <- sum((w_test1*RX1_test+w_test2*RX2_test)*  (RX1_test)  )
      V_num[["roseJ2"]][k] <- sum((w_test1*RX1_test+w_test2*RX2_test)^2 * epsq_test )

      w_test1_ONLY1 <- rowMeans(w_test_rf_1_ONLY1_NUM)/rowMeans(w_test_rf_1_ONLY1_DEN)
      beta_hat_num[["roseJ1"]][k] <- sum((w_test1_ONLY1*RX1_test)*  (RY_test)  )
      beta_hat_den[["roseJ1"]][k] <- sum((w_test1_ONLY1*RX1_test)*  (RX1_test)  )
      V_num[["roseJ1"]][k] <- sum((w_test1_ONLY1*RX1_test)^2 * epsq_test )

      # Unweighted estimator
      beta_hat_num[["unw"]][k] <- sum(RY_test*RX1_test)
      beta_hat_den[["unw"]][k] <- sum(RX1_test^2)
      V_num[["unw"]][k] <- sum(RX1_test^2*epsq_test)

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

