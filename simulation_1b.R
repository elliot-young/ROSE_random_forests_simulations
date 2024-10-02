# Reproducible code for Simulation 1b
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
params_D = list(m.try=ceiling(sqrt(10)), min.node.size=10, sample.fraction=0.05)
params_Y = list(m.try=ceiling(sqrt(11)), min.node.size=10, sample.fraction=1.0)
params_GY = list(m.try=ceiling(sqrt(10)), min.node.size=10, sample.fraction=0.1)
params_W = list(m.try=ceiling(sqrt(11)), min.node.size=100, sample.fraction=0.05)
params_WD = list(m.try=ceiling(sqrt(10)), min.node.size=10, sample.fraction=0.05)
sqrt_link_inverse = function(x) x^2
sqrt_link_deriv = function (x) 1/2 * x^(-1/2)
link_itself = function(x) sqrt(x)

# Node labels (for rose random forests)
node_labels <- list()
for (maxdepth in 1:9) {
  node_labels[[maxdepth]] <- (2^(maxdepth)) : (2^(maxdepth+1)-1)
}
candidate.variables <- c("X.1","X.2","X.3","X.4","X.5","X.6","X.7","X.8","X.9","X.10")
m.try <- ceiling(sqrt(length(candidate.variables)))

# Data generating distribution (Simulation 1b)
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
  mu <- m_0(X)
  sigma_sq <- B/p
  shape <- mu^2/sigma_sq
  scale <- sigma_sq/mu
  D <- rep(Inf,n)
  for (i in 1:n) {
    if (sigma_sq[i]==0) D[i] = mu[i]
    if (sigma_sq[i]>0) D[i] <- rgamma(1, shape = shape[i], scale = scale[i]) #statmod::rinvgauss(1,mean=mu[i], shape=lambda[i]) #
  }
  eta <- beta*D + 0 +1*g_0(X)
  mu <- (eta)^2
  sigma_sq <- 0.01 * mu^2 * B / sqrt(p)
  shape <- mu^2/sigma_sq
  scale <- sigma_sq/mu
  Y <- rep(Inf,n)
  for (i in 1:n) {
    if (sigma_sq[i]==0) Y[i] = mu[i]
    if (sigma_sq[i]>0) Y[i] <- rgamma(1, shape = shape[i], scale = scale[i])
  }
  rdf <- data.frame(X=X, D=D, Y=Y)
  return(rdf)
}
betas = vars = cvs_cart = data.frame()

while ((dim(betas)[1])<4000) {
  estimators <- c("semieff", "unw", paste0("rose", 1:10), paste0("loceff", 1:10))
  CV_cart <- CV_rose <- list()
  for (estimator in paste0("loceff", 0:10)) CV_cart[[estimator]] <- 0
  for (estimator in paste0("rose", 0:10)) CV_rose[[estimator]] <- 0

  set.seed(1 + dim(betas)[1])
  rdf <- generate_data(n=20000)
  K <- 2
  index <- createFolds(rdf[,1], k=K)

  f_w_test <- list()
  estimators <- c("semieff", "unw", paste0("rose", 1:10), paste0("loceff", 1:10))
  for (estimator in estimators) {
    f_w_test[[estimator]] <- list()
  }
  res_x_on_z_theta_k <- y_theta_k <- fitted_Gy_on_z_theta_k <- list()

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
    set.seed(1) # reproducibility of E(Y|X,Z) randomforest
    l_first_model <- ranger(Y ~ D+X.1+X.2+X.3+X.4+X.5+X.6+X.7+X.8+X.9+X.10, data=rdf_train, min.node.size=params_Y$min.node.size, sample.fraction=params_Y$sample.fraction, mtry=params_Y$m.try, num.trees = 500)
    G_dataset <- cbind(rdf_train, G_vals=link_itself(predict(l_first_model, rdf_train)$predictions))
    set.seed(1) # reproducibility of E(g(E(Y|X,Z))|Z) randomforest
    l_second_model <- ranger(G_vals ~ X.1+X.2+X.3+X.4+X.5+X.6+X.7+X.8+X.9+X.10, data=G_dataset, min.node.size=params_GY$min.node.size, sample.fraction=params_GY$sample.fraction, mtry=params_GY$m.try, num.trees = 500)
    fitted_G_on_z_train <- predict(l_second_model, rdf_train)$predictions
    fitted_G_on_z_test <- predict(l_second_model, rdf_test)$predictions

    RD_train <- rdf_train$D - predict(m_model, rdf_train)$predictions
    RD_test <- rdf_test$D - predict(m_model, rdf_test)$predictions

    eval_theta_k <- function(theta) {
      mu <- sqrt_link_inverse( RD_train*theta + fitted_G_on_z_train )
      eps <- sqrt_link_deriv(mu) * ( rdf_train$Y - mu )
      score <- sum( RD_train * eps )/length(RD_train)
      return(score)
    }
    beta_train <- uniroot(eval_theta_k,interval=c(0,2))$root

    mu_nuisance <- sqrt_link_inverse( RD_train*beta_train + fitted_G_on_z_train )
    d_theta_eps_nuisance <- -RD_train
    eps_nuisance <- sqrt_link_deriv(mu_nuisance) * ( rdf_train$Y - mu_nuisance )
    eps_sq_nuisance <- eps_nuisance^2

    res_x_on_z_theta_k[[k]] <- RD_test
    y_theta_k[[k]] <- rdf_test$Y
    fitted_Gy_on_z_theta_k[[k]] <- fitted_G_on_z_test

    epsq <- eps_sq_nuisance
    xisq <- RD_train^2
    xiepsq <- RD_train^2 * eps_sq_nuisance

    # semieff: For the efficient influence function (semieffly biased) estimator
    set.seed(1) # reproducibility of E(\sigma^2(D,X)|X) randomforest
    w_fit <- ranger(epsq ~ D + X.1+X.2+X.3+X.4+X.5+X.6+X.7+X.8+X.9+X.10, data=cbind(epsq,rdf_train), min.node.size=params_W$min.node.size, sample.fraction=params_W$sample.fraction, mtry=params_W$m.try, num.trees = 500)
    w_test <- 1/predict(w_fit, rdf_test)$predictions
    w_train <- 1/predict(w_fit, rdf_train)$predictions
    h <- predict(ranger(D ~ X.1+X.2+X.3+X.4+X.5+X.6+X.7+X.8+X.9+X.10, data=cbind(WD=w_train*rdf_train$D, rdf_train), min.node.size=params_WD$min.node.size, sample.fraction=params_WD$sample.fraction, mtry=params_WD$m.try, num.trees = 500, case.weights=w_train), rdf_test)$predictions
    RtildeD <- rdf_test$D - h
    f_w_test[["semieff"]][[k]] <- w_test*RtildeD

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
      DTdata <- DTdata[sample(seq_len(dim(rdf_train)[1]), floor(m.trydim(rdf_train)[1]), replace=T),]#bootstrap
      w_fit <- rpart(paste("cbind(xisq, xiepsq) ~ ",paste(variables,collapse="+"), collapse=""), data = DTdata, method = ulist_rf, cp=0, 20=20, maxdepth=10)
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
      f_w_test[[paste0("rose", depths)]][[k]] <- w_test[[depths]]*RD_test
    }


    # Fit CART random forest
    eval_theta_k <- function(theta) {
      mu <- sqrt_link_inverse( RD_test*theta + fitted_G_on_z_test )
      eps <- sqrt_link_deriv(mu) * ( rdf_test$Y - mu )
      score <- sum( RD_test * eps )
      return(score)
    }
    beta_test <- uniroot(eval_theta_k,interval=c(0,2))$root
    ep_test <- sqrt_link_deriv(sqrt_link_inverse( RD_test*beta_test + fitted_G_on_z_test )) * ( rdf_test$Y - sqrt_link_inverse( RD_test*beta_test + fitted_G_on_z_test ) )

    preds0 = mean(epsq)
    CV_cart[[paste0("loceff",0)]] = CV_cart[[paste0("loceff",0)]] + sum((preds0-(ep_test)^2)^2)
    # Fit CART forests
    w_test_cart_rf <- list()
    for (depths in 1:10) w_test_cart_rf[[depths]] <- matrix(0,dim(rdf_test)[1],500)
    for (fsts in seq_len(500)) {
      set.seed(fsts) # Reproducibility of rose forest
      variables <- sample(candidate.variables,m.try,replace=FALSE)
      DTdata <- cbind(epsq=epsq,rdf_train)
      DTdata <- DTdata[sample(seq_len(dim(rdf_train)[1]), floor(m.trydim(rdf_train)[1]), replace=T),]#bootstrap
      w_cart_fit <- rpart(paste("cbind(epsq) ~ ",paste(variables,collapse="+"), collapse=""), data = DTdata, method = rsqulist, cp=0, 20=20, maxdepth=10)
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
      CV_cart[[paste0("loceff",depths)]] = CV_cart[[paste0("loceff",depths)]] + sum((preds[[depths]]-(ep_test)^2)^2)
      f_w_test[[paste0("loceff", depths)]][[k]] <- w_test[[depths]]*RD_test
    }

    # Unweighted estimator
    f_w_test[["unw"]][[k]] <- RD_test
  }

  min_theta <- function(estimor) {
    eval_theta_full <- function(theta) {
      mu <- sqrt_link_inverse( unlist(res_x_on_z_theta_k)*theta + unlist(fitted_Gy_on_z_theta_k) )
      eps <- sqrt_link_deriv(mu) * ( unlist(y_theta_k) - mu )
      score <- sum( unlist(f_w_test[[estimor]]) * eps )
      return(score)
    }
    theta_hat <- uniroot(eval_theta_full,interval=c(0,2))$root
    mu_theta <- sqrt_link_inverse( unlist(res_x_on_z_theta_k)*theta_hat + unlist(fitted_Gy_on_z_theta_k) )
    eps <- sqrt_link_deriv(mu_theta) * ( unlist(y_theta_k) - mu_theta )
    psi <- unlist(f_w_test[[estimor]]) * eps
    psi_sq <- psi^2
    d_theta_eps <- - unlist(res_x_on_z_theta_k)
    d_theta_psi <- unlist(f_w_test[[estimor]]) * d_theta_eps
    V_hat <- nrow(rdf) * sum(psi_sq) / (sum(d_theta_psi))^2
    allbetas_semieff = theta_hat
    vars_estimor = V_hat # sqrt(V_hat/nrow(rdf)) for standard deviation
    return(list(theta=theta_hat, var=vars_estimor))
  }
  outputs_list <- Map(min_theta, estimators)
  outputs_theta <- sapply(outputs_list, function(x) x$theta)
  betas_single <- as.data.frame(t(unlist(outputs_theta)))
  outputs_var <- sapply(outputs_list, function(x) x$var)
  vars_single <- as.data.frame(t(unlist(outputs_var)))
  betas = rbind(betas, betas_single)
  vars = rbind(vars, vars_single)
  cvs_cart = rbind(cvs_cart, as.data.frame(t(unlist(CV_cart))))

  write.csv(betas,paste0("~/Documents/Simulation 1b/betas.csv"))
  write.csv(vars,paste0("~/Documents/Simulation 1b/vars.csv"))

}



















for (i in seq_len(I.evalrand)) {
  n <- nis[i]
  index_tracker2 <- index_tracker1 + n - 1
  sub_bit = index_tracker1:index_tracker2
  nodes_i = nodesis[sub_bit]
  
  phi = theta/(1+theta*n)
  vec_sub = numeric(num_leaves)
  for (j in 1:n) {
    vec_sub[nodes_i[j]] = vec_sub[nodes_i[j]] + 1
  }
  unique_nodes_i = unique(nodes_i)
  for (k in unique_nodes_i) {
    for (kk in unique_nodes_i) {
      if (k==kk) {
        XWX[k,kk] = XWX[k,kk] + vec_sub[k]
      }
      XWX[k,kk] = XWX[k,kk] - phi*vec_sub[k]*vec_sub[kk]
    }
  }
  
  index_tracker1 = index_tracker2 + 1
}
