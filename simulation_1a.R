# Reproducible code for Simulation 1a
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

  # For loceff CART random forest:
  rsqitemp <- function(y, offset, parms, wt) {
    if (is.matrix(y) && ncol(y) > 1)
      stop("Matrix response not allowed")
    if (!missing(parms) && length(parms) > 0)
      warning("parameter argument ignored")
    if (length(offset)) y <- y - offset
    sfun <- function(yval, dev, wt, ylevel, digits ) {
      paste(" mean=", format(signif(yval, digits)), ", MSE=" , format(signif(dev/wt, digits)), sep = '')
    }
    environment(sfun) <- .GlobalEnv
    list(y = c(y), parms = NULL, numresp = 1, numy = 1, summary = sfun)
  }
  rsqetemp <- function(y, wt, parms) {
    wmean <- sum(y*wt)/sum(wt)
    cart <- sum(wt*(y-wmean)^2)
    list(label = wmean, deviance = cart)
  }
  rsqstemp <- function(y, wt, x, parms, continuous) {
    # Center y
    n <- length(y)
    y <- y- sum(y*wt)/sum(wt)
    if (continuous) {
      # continuous x variable
      temp <- cumsum(y*wt)[-n]
      left.wt  <- cumsum(wt)[-n]
      right.wt <- sum(wt) - left.wt
      lmean <- temp/left.wt
      rmean <- -temp/right.wt
      goodness <- (left.wt*lmean^2 + right.wt*rmean^2)/sum(wt*y^2)
      list(goodness = goodness, direction = sign(lmean))
    } else {
      # Categorical X variable
      ux <- sort(unique(x))
      wtsum <- tapply(wt, x, sum)
      ysum  <- tapply(y*wt, x, sum)
      means <- ysum/wtsum

      # For anova splits, we can order the categories by their means
      #  then use the same code as for a non-categorical
      ord <- order(means)
      n <- length(ord)
      temp <- cumsum(ysum[ord])[-n]
      left.wt  <- cumsum(wtsum[ord])[-n]
      right.wt <- sum(wt) - left.wt
      lmean <- temp/left.wt
      rmean <- -temp/right.wt
      list(goodness= (left.wt*lmean^2 + right.wt*rmean^2)/sum(wt*y^2),
           direction = ux[ord])
    }

  }
  rsqulist <- list(eval = rsqetemp, split = rsqstemp, init = rsqitemp)

  # For rose random forest:
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
    avar = 10e6*length(y[,1]) - sum_xisq*sum_xisq/sum_xiepsq
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

# Hyperparameters for nuisance regressions (see simulation 1 hyperparameters cross validation)
params_D = list(m.try=ceiling(sqrt(10)), min.node.size=20, sample.fraction=0.05)
params_Y = list(m.try=ceiling(sqrt(10)), min.node.size=20, sample.fraction=0.05)
params_W = list(m.try=ceiling(sqrt(11)), min.node.size=100, sample.fraction=1.0)
params_WD = list(m.try=ceiling(sqrt(10)), min.node.size=20, sample.fraction=0.8)

# Node labels (for rose random forests)
node_labels <- list()
for (maxdepth in 1:9) {
  node_labels[[maxdepth]] <- (2^(maxdepth)) : (2^(maxdepth+1)-1)
}
candidate.variables <- c("X.1","X.2","X.3","X.4","X.5","X.6","X.7","X.8","X.9","X.10")
m.try <- ceiling(sqrt(length(candidate.variables)))

# Data generating distribution (Simulation 1a)
expit = function(x) exp(x)/(1+exp(x))
generate_data <- function(n) {
  covmat <- toeplitz(ARMAacf(ar = c(0.9), lag.max = 10-1))
  g_0 = function(x) expit(x[,1]) + expit(x[,2]) + expit(x[,3]) + expit(x[,4]) + expit(x[,5])
  m_0 = function(x) expit(x[,1]) + expit(x[,2]) + expit(x[,3]) + expit(x[,4]) + expit(x[,5])
  beta = 1
  X = as.matrix(MASS::mvrnorm(n = n, rep(0,10), covmat))
  p <- (exp(3*X[,1])/(1+exp(3*X[,1])))
  p[p<0.01] <- 0.01
  B <- rbinom(n, 1, p)
  zeta <- B/p - 1
  xi <- sqrt(1+zeta)*rnorm(n)
  epsilon <- sqrt(1+zeta)*(p^(1/4))*rnorm(n)
  D <- m_0(X) + xi
  Y <- beta*D + g_0(X) + epsilon
  rdf <- data.frame(X=X, D=D, Y=Y)
  return(rdf)
}

# Run Simulation 1a
betas = vars = cvs_cart = data.frame()
while ((dim(betas)[1])<4000) {
  estimators <- c("semieff", "unw", paste0("rose", 1:10), paste0("loceff", 1:10))
  CV_cart <- CV_rose <- list()
  for (estimator in paste0("loceff", 0:10)) CV_cart[[estimator]] <- 0
  for (estimator in paste0("rose", 0:10)) CV_rose[[estimator]] <- 0

  set.seed(1 + dim(betas)[1])
  rdf <- generate_data(n=20000)
  K <- 2
  index <- caret::createFolds(rdf[,1], k=K)

  beta_hat_num <- beta_hat_den <- V_num <- list()
  estimators <- c("semieff", "unw", paste0("rose", 1:10), paste0("loceff", 1:10))
  for (estimator in estimators) {
    beta_hat_num[[estimator]] <- numeric(K)
    beta_hat_den[[estimator]] <- numeric(K)
    V_num[[estimator]] <- numeric(K)
  }

  for (k in seq_len(K)) {
    if (k==1) {
      rdf_train <- rdf[index$Fold1,] # for K=2 only
      rdf_test <- rdf[index$Fold2,] # for K=2 only
    } else if (k==2) {
      rdf_train <- rdf[index$Fold2,] # for K=2 only
      rdf_test <- rdf[index$Fold1,] # for K=2 only
    }


    set.seed(1) # reproducibility of E(D|X) randomforest
    m_model <- ranger(D ~ X.1+X.2+X.3+X.4+X.5+X.6+X.7+X.8+X.9+X.10, data=rdf_train, min.node.size=params_D$min.node.size, sample.fraction=params_D$sample.fraction, mtry=params_D$m.try, num.trees = 500)
    set.seed(1) # reproducibility of E(Y|X) randomforest
    l_model <- ranger(Y ~ X.1+X.2+X.3+X.4+X.5+X.6+X.7+X.8+X.9+X.10, data=rdf_train, min.node.size=params_Y$min.node.size, sample.fraction=params_Y$sample.fraction, mtry=params_Y$m.try, num.trees = 500)

    RD_train <- rdf_train$D - predict(m_model, rdf_train)$predictions
    RD_test <- rdf_test$D - predict(m_model, rdf_test)$predictions
    RY_train <- rdf_train$Y - predict(l_model, rdf_train)$predictions
    RY_test <- rdf_test$Y - predict(l_model, rdf_test)$predictions
    beta_train <- sum(RY_train*RD_train)/sum(RD_train^2)
    eps_hat <- RY_train - beta_train * RD_train
    xi_hat <- RD_train

    epsq <- eps_hat^2
    xisq <- xi_hat^2
    xiepsq <- xi_hat^2*eps_hat^2

    # For the efficient influence function estimator
    set.seed(1) # reproducibility of E(\sigma^2(D,X)|X) randomforest
    w_fit <- ranger(epsq ~ D + X.1+X.2+X.3+X.4+X.5+X.6+X.7+X.8+X.9+X.10, data=cbind(epsq,rdf_train), min.node.size=params_W$min.node.size, sample.fraction=params_W$sample.fraction, mtry=params_W$m.try, num.trees = 500)
    w_test <- 1/predict(w_fit, rdf_test)$predictions
    w_train <- 1/predict(w_fit, rdf_train)$predictions
    set.seed(1)  # reproducibility of h(X) randomforest
    h <- predict(ranger(D ~ X.1+X.2+X.3+X.4+X.5+X.6+X.7+X.8+X.9+X.10, data=cbind(WD=w_train*rdf_train$D, rdf_train), min.node.size=params_WD$min.node.size, sample.fraction=params_WD$sample.fraction, mtry=params_WD$m.try, num.trees = 500, case.weights=w_train), rdf_test)$predictions
    RtildeD <- rdf_test$D - h
    beta_hat_num[["semieff"]][k] <- sum(w_test*RY_test*RtildeD)
    beta_hat_den[["semieff"]][k] <- sum(w_test*RD_test*RtildeD)
    V_num[["semieff"]][k] <- sum(w_test^2*(RY_test-beta_train*RD_test)^2*RtildeD^2)

    # Fit rose forest
    w_test_rf_NUM <- w_test_rf_DEN <- list()
    for (depths in 1:10) {
      w_test_rf_NUM[[depths]] <- matrix(0,dim(rdf_test)[1],500)
      w_test_rf_DEN[[depths]] <- matrix(0,dim(rdf_test)[1],500)
    }
    for (fsts in seq_len(500)) {
      set.seed(fsts) # Reproducibility of rose forest
      variables <- sample(candidate.variables,m.try,replace=FALSE)
      DTdata <- cbind(xisq=xisq,xiepsq=xiepsq,rdf_train)
      DTdata <- DTdata[sample(seq_len(dim(rdf_train)[1]), floor(dim(rdf_train)[1]), replace=T),]#bootstrap
      w_fit <- rpart(paste("cbind(xisq, xiepsq) ~ ",paste(variables,collapse="+"), collapse=""), data = DTdata, method = ulist_rf, cp=0, minbucket=20, maxdepth=10)
      for (depths in 1:10) {
        if (depths < 10) {
          w_fit_snipped <- snip.rpart(w_fit, toss=node_labels[[depths]])
        } else {
          w_fit_snipped <- w_fit
        }
        w_test_rf_NUM[[depths]][,fsts] <- predict(w_fit_snipped, rdf_test, type="matrix")[,1]
        w_test_rf_DEN[[depths]][,fsts] <- predict(w_fit_snipped, rdf_test, type="matrix")[,2]
      }
    }
    w_test <- Map("/", lapply(w_test_rf_NUM, rowMeans), lapply(w_test_rf_DEN, rowMeans))
    for (depths in 1:10) {
      beta_hat_num[[paste0("rose", depths)]][k] <- sum(w_test[[depths]]*RY_test*RD_test)
      beta_hat_den[[paste0("rose", depths)]][k] <- sum(w_test[[depths]]*RD_test*RD_test)
      V_num[[paste0("rose", depths)]][k] <- sum(w_test[[depths]]^2*(RY_test-beta_train*RD_test)^2*RD_test^2)
    }


    # Fit loceff CART random forest
    preds0 = mean(epsq) # for a CART CV criterion
    CV_cart[[paste0("loceff",0)]] = CV_cart[[paste0("loceff",0)]] + sum((preds0-(RY_test-beta_train*RD_test)^2)^2)
    # Fit loceff CART random forests
    w_test_cart_rf <- list()
    for (depths in 1:10) w_test_cart_rf[[depths]] <- matrix(0,dim(rdf_test)[1],500)
    for (fsts in seq_len(500)) {
      set.seed(fsts) # Reproducibility of CART forest
      variables <- sample(candidate.variables,m.try,replace=FALSE)
      DTdata <- cbind(epsq=epsq,rdf_train)
      DTdata <- DTdata[sample(seq_len(dim(rdf_train)[1]), floor(dim(rdf_train)[1]), replace=T),]#bootstrap
      w_cart_fit <- rpart(paste("cbind(epsq) ~ ",paste(variables,collapse="+"), collapse=""), data = DTdata, method = rsqulist, cp=0, minbucket=20, maxdepth=10)
      for (depths in 1:10) {
        if (depths < 10) {
          w_fit_snipped <- snip.rpart(w_cart_fit, toss=node_labels[[depths]])
        } else {
          w_fit_snipped <- w_cart_fit
        }
        w_test_cart_rf[[depths]][,fsts] <- predict(w_fit_snipped, rdf_test, type="matrix")
      }
    }
    preds <- lapply(w_test_cart_rf, rowMeans)
    w_test <- lapply(preds, function(x) 1/x)
    for (depths in 1:10) {
      CV_cart[[paste0("loceff",depths)]] = CV_cart[[paste0("loceff",depths)]] + sum((preds[[depths]]-(RY_test-beta_train*RD_test)^2)^2)
      beta_hat_num[[paste0("loceff", depths)]][k] <- sum(w_test[[depths]]*RY_test*RD_test)
      beta_hat_den[[paste0("loceff", depths)]][k] <- sum(w_test[[depths]]*RD_test*RD_test)
      V_num[[paste0("loceff", depths)]][k] <- sum(w_test[[depths]]^2*(RY_test-beta_train*RD_test)^2*RD_test^2)
    }

    # Unweighted estimator
    beta_hat_num[["unw"]][k] <- sum(RY_test*RD_test)
    beta_hat_den[["unw"]][k] <- sum(RD_test*RD_test)
    V_num[["unw"]][k] <- sum((RY_test-beta_train*RD_test)^2*RD_test^2)
  }

  beta_single <- Map(function(x1,x2) sum(x1)/sum(x2), beta_hat_num, beta_hat_den)
  betas_single <- as.data.frame(t(unlist(beta_single)))
  var_single <- Map(function(x1,x2) n*sum(x1)/(sum(x2))^2, V_num, beta_hat_den)
  vars_single <- as.data.frame(t(unlist(var_single)))
  betas = rbind(betas, betas_single)
  vars = rbind(vars, vars_single)
  cvs_cart = rbind(cvs_cart, as.data.frame(t(unlist(CV_cart))))

  write.csv(betas,paste0("~/Documents/Simulation 1a/betas.csv"))
  write.csv(vars,paste0("~/Documents/Simulation 1a/vars.csv"))
}
