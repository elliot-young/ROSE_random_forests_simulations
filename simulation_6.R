# Reproducible code for Simulation 6
# Note that one can fit a rose random forest estimator directly using [https://github.com/elliot-young/rose].
#   In this simulation, the four estimators share a number of nuisance functions,
#   and so for computational speed these estimators are estimated alongside
#   one another (as below).

# Load relevant rose random forest functions - J=1
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
#Load relevant rose random forest functions - J=2
{
  library(rpart)
  library(ranger)
  library(mgcv)
  library(glmnet)
  J=2

  itemp_rf_J2 <- function(y, offset, parms, wt) {
    if (ncol(y) != 6) {
      stop("Matrix of response must be a 6-column matrix")
    }
    if (!missing(parms) && length(parms) > 0) {
      warning("parameter argument ignored")
    }
    if (length(offset)) y <- y - offset

    sfun <- function(weigh, avar, ylevel, digits) {
      paste(" Score=", format(signif(-avar, digits)), sep = '')
    }
    environment(sfun) <- .GlobalEnv

    list(y = y, parms = NULL, numresp = 6, numy = 6, summary = sfun)
  }

  etemp_rf_J2 <- function(y, wt, parms) {
    sums <- colSums(y)
    detV <- sums[3] * sums[6] - sums[4] * sums[5]

    score <- if (detV > 1e-10) {
      (sums[1]^2 * sums[6] - sums[1] * sums[2] * (sums[4] + sums[5]) + sums[2]^2 * sums[3]) / detV
    } else {
      0
    }

    # Deviancy is minimized by rpart, so we negate the score to maximize it
    list(label = sums, deviance = -score)
  }
  etemp_rf_J2 <- function(y, wt, parms) {
    n <- nrow(y)
    sums <- colSums(y)
    detV <- sums[3] * sums[6] - sums[4] * sums[5]

    score <- if (detV > 1e-10) {
      (sums[1]^2 * sums[6] - sums[1] * sums[2] * (sums[4] + sums[5]) + sums[2]^2 * sums[3]) / detV
    } else {
      0
    }

    # Ensure deviance is strictly positive so cp checks pass
    deviance_val <- 10000 * n - score

    list(label = sums, deviance = deviance_val)
  }

  stemp_rf_J2 <- function(y, wt, x, parms, continuous) {
    n <- dim(y)[1]
    if (continuous) {
      tot <- colSums(y)

      # Cumulative sums across candidate splits
      l1 <- cumsum(y[, 1])[-n]; r1 <- tot[1] - l1
      l2 <- cumsum(y[, 2])[-n]; r2 <- tot[2] - l2
      l3 <- cumsum(y[, 3])[-n]; r3 <- tot[3] - l3
      l4 <- cumsum(y[, 4])[-n]; r4 <- tot[4] - l4
      l5 <- cumsum(y[, 5])[-n]; r5 <- tot[5] - l5
      l6 <- cumsum(y[, 6])[-n]; r6 <- tot[6] - l6

      # Vectorized score evaluator: C^T V^{-1} C
      eval_score <- function(y1, y2, y3, y4, y5, y6) {
        detV <- y3 * y6 - y4 * y5
        ifelse(detV > 1e-10,
               (y1^2 * y6 - y1 * y2 * (y4 + y5) + y2^2 * y3) / detV,
               0)
      }

      l_score <- eval_score(l1, l2, l3, l4, l5, l6)
      r_score <- eval_score(r1, r2, r3, r4, r5, r6)
      tot_score <- eval_score(tot[1], tot[2], tot[3], tot[4], tot[5], tot[6])

      goodness <- l_score + r_score - tot_score
      list(goodness = goodness, direction = rep(1, length(x) - 1))
    } else {
      stop("Not currently supported for discrete X variables")
    }
  }

  ulist_rf_J2 <- list(eval = etemp_rf_J2, split = stemp_rf_J2, init = itemp_rf_J2)
  #ulist_rf <- list(eval = etemp_rf, split = stemp_rf, init = itemp_rf)

}

n=1000
generate_data = function(n) {
  m <- function(z) 4*exp(2*z)/(1+exp(2*z))
  f <- function(z) 4*z^2/(1+z^2)
  X <- 1 + m(Z) + rgamma(n, 1, 20)
  Y <- X + f(Z) + 1.0*rnorm(n,0,sd = sqrt(X/10))
  rdf <- data.frame(Z=Z, X=X, Y=Y, logX=((X-min(X))/(max(X)-min(X)))^2)
  return(rdf)
}




numoftrees = 100
minbucket = 10
maxdepth = 5

betas <- vars <- vars2 <- data.frame()
estimators <- c("semieff", "unw", "roseJ1", "roseJ2")

# Helper functions
compute_optimal_weights <- function(A, normalize = c("none", "sum1", "unit_score")) {
  normalize <- match.arg(normalize)

  Y1 <- A[, 1]
  Y2 <- A[, 2]
  Y3 <- A[, 3]
  Y4 <- (A[, 4] + A[, 5]) / 2
  Y6 <- A[, 6]

  detV <- Y3 * Y6 - Y4^2

  w1_raw <- (Y6 * Y1 - Y4 * Y2) / detV
  w2_raw <- (Y3 * Y2 - Y4 * Y1) / detV

  invalid <- detV <= 1e-10
  w1_raw[invalid] <- 0
  w2_raw[invalid] <- 0

  if (normalize == "sum1") {
    total <- w1_raw + w2_raw
    w1 <- w1_raw / total
    w2 <- w2_raw / total
  } else if (normalize == "unit_score") {
    score <- w1_raw * Y1 + w2_raw * Y2
    w1 <- w1_raw / score
    w2 <- w2_raw / score
  } else {
    w1 <- w1_raw
    w2 <- w2_raw
  }

  data.frame(w1 = w1, w2 = w2)
}


while (nrow(betas) < 1000) {print(nrow(betas))# 4000
  beta_hat_num <- beta_hat_den <- beta_hat_den2 <- V_num <- list()
  for (estimator in estimators) {
    beta_hat_num[[estimator]] <- numeric(K)
    beta_hat_den[[estimator]] <- numeric(K)
    beta_hat_den2[[estimator]] <- numeric(K)
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
    YXZ_gam_model <- mgcv::gam(Y ~ X + s(Z,bs="cr", k=20), data = rdf_train)
    YXZ_predicted_terms <- predict(YXZ_gam_model, newdata = rdf_test, type = "terms")
    f_Z_test <- YXZ_predicted_terms[, "s(Z)"]
    model_intercept <- attr(YXZ_predicted_terms, "constant")
    f_Z_test <- f_Z_test + model_intercept

    #Y_m_f_train <- rdf_train$Y - predict(YXZ_gam_model, newdata=rdf_train)
    Y_m_f_test <- rdf_test$Y - f_Z_test


    # E(X|Z)
    X_model <- mgcv::gam(X ~ s(Z,bs="cr", k=20), data = rdf_train)
    RX_train <- rdf_train$X - predict(X_model, newdata=rdf_train)
    RX_test <- rdf_test$X -predict(X_model, newdata=rdf_test)

    # E(log(X)|Z)
    sqX_model <- mgcv::gam(sqX ~ s(Z,bs="cr", k=20), data=rdf_train)
    logRX_train <- rdf_train$sqX - predict(sqX_model, newdata=rdf_train)
    logRX_test <- rdf_test$sqX - predict(sqX_model, newdata=rdf_test)

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
    xi1xi2epsq <- RX_train * logRX_train * epsq

    DTdata <- cbind(xisq1=xisq1,xiepsq1=xiepsq1, xisq2=xisq2,xiepsq2=xiepsq2, rdf_train, xi1=xi1, xi2=xi2, epsq=epsq, xi1xi2epsq=xi1xi2epsq)

    DT_sub <- data.frame(resp = 1/DTdata$epsq, wts = DTdata$xisq1*DTdata$xiepsq1)

    # library(grf)

    #ROSErf <- regression_forest(
    #  DTdata[,5:7],
    #  DT_sub$resp,
    #  sample.weights = DT_sub$wts,
    #  num.trees = 500,
    #  honesty=FALSE,
    #  min.node.size=100
    #)
    #
    # w_test1 <- predict(ROSErf, rdf_test[,1:3])
    #
    ##w_test1_ONLY1 <- rowMeans(w_test_rf_1_ONLY1_NUM)/rowMeans(w_test_rf_1_ONLY1_DEN)
    #beta_hat_num[["roseJ1"]][k] <- sum((w_test1*RX1_test)*  (RY_test)  )
    #beta_hat_den[["roseJ1"]][k] <- sum((w_test1*RX1_test)*  (RX1_test)  )
    #V_num[["roseJ1"]][k] <- sum((w_test1*RX1_test)^2 * epsq_test )





    # Fit rose forests (J=1 and J=2 simultaneously)
    forest_fits <- vector("list", numoftrees)
    ###w_test_rf_1 <- w_test_rf_2 <- w_test_rf_3 <- w_test_rf_4 <- w_test_rf_5 <- w_test_rf_6 <- matrix(0,dim(rdf_test)[1],numoftrees)
    #w_test_rf_1_ONLY1_NUM <- matrix(0,dim(rdf_test)[1],numoftrees)
    #w_test_rf_1_ONLY1_DEN <- matrix(0,dim(rdf_test)[1],numoftrees)
    for (fsts in seq_len(numoftrees)) {
      set.seed(fsts) # Reproducibility of rose forest
      # Trees
      selection <- sample(seq_len(dim(rdf_train)[1]), floor(subsizeprop*dim(rdf_train)[1]), replace=T)
      DTdata_sub <- DTdata[selection,]#bootstrap
      w_fit_1 <- w_fit_1_CORE <- rpart(paste("cbind(xisq1, xisq2, xiepsq1, xi1xi2epsq, xi1xi2epsq, xiepsq2) ~ ",paste(c("Z"),collapse="+"), collapse=""), data = DTdata_sub, method = ulist_rf_J2, cp=0, minbucket=minbucket, maxdepth=max(c(maxdepth)))
      #w_fit_1 <- snip.rpart(w_fit_1_CORE, toss=node_labels[[maxdepth]])
      #w_fit_1_ONLY1 <- snip.rpart(w_fit_1_CORE, toss=node_labels[[maxdepth_only1]])

      #leaves1 <- which(w_fit_1$frame$var=="<leaf>")
      #num_leaves1 <- length(leaves1)

      forest_fits[[fsts]] <- w_fit_1

      #predict(w_fit_1, rdf_test, type="matrix")[,1]
      #A = predict(w_fit_1, rdf_test, type="matrix")

      # Usage:
      #A <- predict(w_fit_1, rdf_test, type = "matrix")
      #weights_df <- compute_optimal_weights(A, normalize = "sum1")
    }

    preds_list <- lapply(forest_fits, function(tree) {
      predict(tree, rdf_test, type = "matrix")
    })

    # 4. Aggregate component sums across all trees
    A_forest <- Reduce("+", preds_list) / length(preds_list)

    # 5. Extract optimal aggregated weights
    optimal_weights <- compute_optimal_weights(A_forest, normalize = "none")

    #w_test1 <- rowMeans(w_test_rf)
    beta_hat_num[["roseJ2"]][k] <- sum((optimal_weights[,1]*RX1_test+optimal_weights[,2]*RX2_test)*  (Y_m_f_test)  )
    beta_hat_den[["roseJ2"]][k] <- sum((optimal_weights[,1]*RX1_test+optimal_weights[,2]*RX2_test)*  (rdf_test$X)  )
    V_num[["roseJ2"]][k] <- sum((optimal_weights[,1]*RX1_test+optimal_weights[,2]*RX2_test)^2 * epsq_test )

    ##beta_hat_num[["roseJ1"]][k] <- sum((w_test1_ONLY1*RX1_test)*  (RY_test)  )
    ##beta_hat_den[["roseJ1"]][k] <- sum((w_test1_ONLY1*RX1_test)*  (RX1_test)  )
    ##V_num[["roseJ1"]][k] <- sum((w_test1_ONLY1*RX1_test)^2 * epsq_test )










    w_test_rf_1_ONLY1_NUM <- matrix(0,dim(rdf_test)[1],numoftrees)
    w_test_rf_1_ONLY1_DEN <- matrix(0,dim(rdf_test)[1],numoftrees)
    for (fsts in seq_len(numoftrees)) {
      set.seed(fsts) # Reproducibility of rose forest
      # Trees
      selection <- sample(seq_len(dim(rdf_train)[1]), floor(subsizeprop*dim(rdf_train)[1]), replace=T)
      DTdata_sub <- DTdata[selection,]#bootstrap
      w_fit_1_ONLY1 <- rpart(paste("cbind(xisq1, xiepsq1) ~ ",paste(c("Z"),collapse="+"), collapse=""), data = DTdata_sub, method = ulist_rf, cp=0, minbucket=minbucket, maxdepth=max(c(maxdepth)))

      w_test_rf_1_ONLY1_NUM[,fsts] <- predict(w_fit_1_ONLY1, rdf_test, type="matrix")[,1]
      w_test_rf_1_ONLY1_DEN[,fsts] <- predict(w_fit_1_ONLY1, rdf_test, type="matrix")[,2]
    }

    w_test1_ONLY1 <- rowMeans(w_test_rf_1_ONLY1_NUM)/rowMeans(w_test_rf_1_ONLY1_DEN)
    beta_hat_num[["roseJ1"]][k] <- sum((w_test1_ONLY1*RX1_test)*  (Y_m_f_test)  )
    beta_hat_den[["roseJ1"]][k] <- sum((w_test1_ONLY1*RX1_test)*  (rdf_test$X)  )
    V_num[["roseJ1"]][k] <- sum((w_test1_ONLY1*RX1_test)^2 * epsq_test )



    plot(rdf_test$Z, w_test1_ONLY1)



    # Semiparametric efficient estimator
    sigmasq_model <- mgcv::gam(epsq ~ s(X,bs="cr", k=20)+s(Z,bs="cr", k=20), data = cbind(epsq,rdf_train))
    W_train <- 1/predict(sigmasq_model, newdata=rdf_train)
    W_test <- 1/predict(sigmasq_model, newdata=rdf_test)
    W_train <- pmax(W_train, 0)
    W_test <- pmax(W_test, 0)

    h_mod <- mgcv::gam(X ~ s(Z,bs="cr", k=20), data = rdf_train, weights = W_train)
    h_test <- predict(h_mod, newdata=rdf_test)
    RtildeX1 <- rdf_test$X - h_test

    beta_hat_num[["semieff"]][k] <- sum(W_test*Y_m_f_test*RtildeX1)
    beta_hat_den[["semieff"]][k] <- sum(W_test*rdf_test$X*RtildeX1)
    beta_hat_den2[["semieff"]][k] <- sum(W_test*RX_test*RtildeX1)
    V_num[["semieff"]][k] <- sum(W_test^2*epsq_test*RtildeX1^2)

    # Unweighted estimator
    beta_hat_num[["unw"]][k] <- sum(Y_m_f_test*RX_test)
    beta_hat_den[["unw"]][k] <- sum(RX_test*rdf_test$X)
    beta_hat_den2[["unw"]][k] <- sum(RX_test^2)
    V_num[["unw"]][k] <- sum(RX_test^2*epsq_test)

  }
  print(beta_hat_num$unw/beta_hat_den$unw)
  print(beta_hat_num$unw)
  print(beta_hat_den$unw)
  print(sqrt(V_num$unw/beta_hat_den$unw^2/n))
  beta_single <- Map(function(x1,x2) mean(x1/x2), beta_hat_num, beta_hat_den)
  betas_single <- as.data.frame(t(unlist(beta_single)))
  var_single <- Map(function(x1,x2) n*mean(x1/x2^2), V_num, beta_hat_den)
  vars_single <- as.data.frame(t(unlist(var_single)))
  var_single2 <- Map(function(x1,x2) n*mean(x1/x2^2), V_num, beta_hat_den2)
  vars_single2 <- as.data.frame(t(unlist(var_single2)))
  #betas <- rbind(betas, betas_single)
  #vars <- rbind(vars, vars_single)



  betas <- rbind(betas, cbind(betas_single))
  vars <- rbind(vars, cbind(vars_single))
  vars2 <- rbind(vars2, cbind(vars_single2))
  if (dim(betas)[1]>=4) {
    plot(density(betas$semieff), ylim=c(0,9), xlim=c(0.5,1.5), col="blue")
    lines(density(betas$unw),col="black")
    lines(density(betas$roseunw5),col="grey")
    lines(density(betas$roseJ22),col="red")
    print(colMeans((betas-1)^2))
    print(colMedians((betas-1)^2))
    #print(colMeans((betas-1))^2/colMeans((betas-1)^2)); print(colMeans((betas-1))^2/colMeans((betas-1)^2));
    print(mean((betas$unw-1-qnorm(0.8)*sqrt(vars$unw/n))*(betas$unw-1+qnorm(0.8)*sqrt(vars$unw/n))<0)); print(mean((betas$roseunw5-1-qnorm(0.8)*sqrt(vars$roseunw5/n))*(betas$roseunw5-1+qnorm(0.8)*sqrt(vars$roseunw5/n))<0)); print(mean((betas$roseJ22-1-qnorm(0.8)*sqrt(vars$roseJ22/n))*(betas$roseJ22-1+qnorm(0.8)*sqrt(vars$roseJ22/n))<0))
  }
}




