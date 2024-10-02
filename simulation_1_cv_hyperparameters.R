# Hyperparameter tuning for random forests in simulation 1
library(ranger)

# ------------------- Simulation 1a -------------------

# Data generating distribution (Simulation 1a)
covmat = toeplitz(ARMAacf(ar=c(0.9), lag.max=10-1))
expit = function(x)  exp(x)/(1+exp(x))
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

# Hyperparameter tuning grid
minimum.node.sizes <- c(10,20,50,100,200,500,1000,2000)
sample.fractions <- c(0.01,0.05,0.1,0.2,0.5,0.8,1.0)

# l(.) and m(.) cross-validation
CVl <- CVm <- matrix(0,length(minimum.node.sizes),length(sample.fractions))
colnames(CVl)=sample.fractions; rownames(CVl)=minimum.node.sizes
colnames(CVm)=sample.fractions; rownames(CVm)=minimum.node.sizes
for (runs in 1:100) {
  rdf <- generate_data(n=20000)
  K <- 2
  index <- rbinom(nrow(rdf), 1, 0.5)

  for (min.node.index in seq_len(length(minimum.node.sizes))) {
    for (sam.frac.index in seq_len(length(sample.fractions))) {
      min.node.here = minimum.node.sizes[min.node.index]
      sam.frac.here = sample.fractions[sam.frac.index]
      print(c(min.node.here,sam.frac.here))

      for (k in seq_len(K)) {
        rdf_train <- rdf[index==k-1,]
        rdf_test <- rdf[index!=k-1,]

        m_model <- ranger(D ~ X.1+X.2+X.3+X.4+X.5+X.6+X.7+X.8+X.9+X.10, data=rdf_train, min.node.size=min.node.here, sample.fraction=sam.frac.here, mtry=ceiling(sqrt(10)), num.trees = 500)
        l_model <- ranger(Y ~ X.1+X.2+X.3+X.4+X.5+X.6+X.7+X.8+X.9+X.10, data=rdf_train, min.node.size=min.node.here, sample.fraction=sam.frac.here, mtry=ceiling(sqrt(10)), num.trees = 500)

        CVm[min.node.index,sam.frac.index] <- CVm[min.node.index,sam.frac.index] + sum((rdf_test$D-predict(m_model,rdf_test)$predictions)^2)
        CVl[min.node.index,sam.frac.index] <- CVl[min.node.index,sam.frac.index] + sum((rdf_test$Y-predict(l_model,rdf_test)$predictions)^2)
      }
    }
  }
}
min.l = which(CVl == min(CVl), arr.ind = TRUE)
l.min.node.size = minimum.node.sizes[min.l[1]]
l.sample.fraction = sample.fractions[min.l[2]]
min.m = which(CVm == min(CVm), arr.ind = TRUE)
m.min.node.size = minimum.node.sizes[min.m[1]]
m.sample.fraction = sample.fractions[min.m[2]]

# sigma_sq(.) cross-validation
CVsigsq <- matrix(0,length(minimum.node.sizes),length(sample.fractions))
colnames(CVsigsq)=sample.fractions; rownames(CVsigsq)=minimum.node.sizes
for (runs in 1:100) {
  rdf <- generate_data(n=20000)
  K <- 2
  index <- rbinom(nrow(rdf), 1, 0.5)

  for (min.node.index in seq_len(length(minimum.node.sizes))) {
    for (sam.frac.index in seq_len(length(sample.fractions))) {
      min.node.here = minimum.node.sizes[min.node.index]
      sam.frac.here = sample.fractions[sam.frac.index]

      for (k in seq_len(K)) {
        rdf_train <- rdf[index==k-1,]
        rdf_test <- rdf[index!=k-1,]

        m_model <- ranger(D ~ X.1+X.2+X.3+X.4+X.5+X.6+X.7+X.8+X.9+X.10, data=rdf_train, min.node.size=m.min.node.size, sample.fraction=m.sample.fraction, mtry=ceiling(sqrt(10)), num.trees = 500)
        l_model <- ranger(Y ~ X.1+X.2+X.3+X.4+X.5+X.6+X.7+X.8+X.9+X.10, data=rdf_train, min.node.size=l.min.node.size, sample.fraction=l.sample.fraction, mtry=ceiling(sqrt(10)), num.trees = 500)

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
        epsq_test <- (RY_test - 1 * RD_test)^2

        w_fit <- ranger(epsq ~ D + X.1+X.2+X.3+X.4+X.5+X.6+X.7+X.8+X.9+X.10, data=cbind(epsq,rdf_train), min.node.size=min.node.here, sample.fraction=sam.frac.here, mtry=ceiling(sqrt(11)), num.trees = 500)
        sig_sq_test <- predict(w_fit, rdf_test)$predictions
        sig_sq_train <- predict(w_fit, rdf_train)$predictions
        CVsigsq[min.node.index,sam.frac.index] <- CVsigsq[min.node.index,sam.frac.index] + sum((epsq_test-sig_sq_test)^2)
      }
    }
  }
}
min.sig = which(CVsigsq == min(CVsigsq), arr.ind = TRUE)
sig.min.node.size = minimum.node.sizes[min.sig[1]]
sig.sample.fraction = sample.fractions[min.sig[2]]

# h(.) cross-validation
CVh <- matrix(0,length(minimum.node.sizes),length(sample.fractions))
colnames(CVh)=sample.fractions; rownames(CVh)=minimum.node.sizes
for (runs in 1:100) {
  print(opopop)
  set.seed(opopop)
  rdf <- generate_data(n=20000)
  K <- 2
  index <- rbinom(nrow(rdf), 1, 0.5)

  for (min.node.index in seq_len(length(minimum.node.sizes))) {
    for (sam.frac.index in seq_len(length(sample.fractions))) {
      min.node.here = minimum.node.sizes[min.node.index]
      sam.frac.here = sample.fractions[sam.frac.index]
      print(c(min.node.here,sam.frac.here))

      for (k in seq_len(K)) {
        rdf_train <- rdf[index==k-1,]
        rdf_test <- rdf[index!=k-1,]

        m_model <- ranger(D ~ X.1+X.2+X.3+X.4+X.5+X.6+X.7+X.8+X.9+X.10, data=rdf_train, min.node.size=m.min.node.size, sample.fraction=m.sample.fraction, mtry=ceiling(sqrt(10)), num.trees = 500)
        l_model <- ranger(Y ~ X.1+X.2+X.3+X.4+X.5+X.6+X.7+X.8+X.9+X.10, data=rdf_train, min.node.size=l.min.node.size, sample.fraction=l.sample.fraction, mtry=ceiling(sqrt(10)), num.trees = 500)

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
        epsq_test <- (RY_test - 1 * RD_test)^2

        w_fit <- ranger(epsq ~ D + X.1+X.2+X.3+X.4+X.5+X.6+X.7+X.8+X.9+X.10, data=cbind(epsq,rdf_train), min.node.size=sig.min.node.size, sample.fraction=sig.sample.fraction, mtry=ceiling(sqrt(11)), num.trees = 500)
        w_test <- 1/predict(w_fit, rdf_test)$predictions
        w_train <- 1/predict(w_fit, rdf_train)$predictions
        h <- predict(ranger(D ~ X.1+X.2+X.3+X.4+X.5+X.6+X.7+X.8+X.9+X.10, data=rdf_train, min.node.size=min.node.here, sample.fraction=sam.frac.here, mtry=ceiling(sqrt(10)), num.trees = 500, case.weights=w_train), rdf_test)$predictions
        CVh[min.node.index,sam.frac.index] <- CVh[min.node.index,sam.frac.index] + sum(w_test*(rdf_test$D-h)^2)
      }
    }
  }
}
min.h = which(CVh == min(CVh), arr.ind = TRUE)
h.min.node.size = minimum.node.sizes[min.h[1]]
h.sample.fraction = sample.fractions[min.h[2]]

# ------------------- Simulation 1b -------------------

# Data generating distribution (Simulation 1b)
covmat = toeplitz(ARMAacf(ar=c(0.9), lag.max=10-1))
expit = function(x)  exp(x)/(1+exp(x))
generate_data <- function(n) {
  g_0 = function(x) expit(x[,1]) + expit(x[,2]) + expit(x[,3]) + expit(x[,4]) + expit(x[,5])
  m_0 = function(x) expit(x[,1]) + expit(x[,2]) + expit(x[,3]) + expit(x[,4]) + expit(x[,5])
  beta = 1
  X = as.matrix(MASS::mvrnorm(n = n, rep(0,10), covmat))
  p <- (exp(aaa*X[,1])/(1+exp(aaa*X[,1])))
  p[p<0.01] <- 0.01
  B <- rbinom(n, 1, p)
  zeta <- B/p - 1
  {
    mu <- m_0(X)
    sigma_sq <- 0.1*mu^2 * B / p
    alpha <- mu^2/sigma_sq
    shape <- alpha
    beta <- mu/sigma_sq
    scale <- 1/beta
  }
  D <- rep(Inf,n)
  for (i in 1:n) {
    if (sigma_sq[i]==0) D[i] = mu[i]
    if (sigma_sq[i]>0) D[i] <- rgamin.node.indexa(1, shape = shape[i], scale = scale[i]) #statmod::rinvgauss(1,mean=mu[i], shape=lambda[i]) #
  }
  beta=1
  eta <- beta*D + 0 +1*g_0(X)
  hist(eta)
  mu <- (eta)^2
  sigma_sq <- 0.01 * mu^2 * B / sqrt(p)
  alpha <- mu^2/sigma_sq
  shape <- alpha
  beta <- mu/sigma_sq
  scale <- 1/beta
  Y <- rep(Inf,n)
  for (i in 1:n) {
    if (sigma_sq[i]==0) Y[i] = mu[i]
    if (sigma_sq[i]>0) Y[i] <- rgamin.node.indexa(1, shape = shape[i], scale = scale[i])
  }
  rdf <- data.frame(X=X, D=D, Y=Y)
  return(rdf)
}

# Hyperparameter tuning grid
minimum.node.sizes <- c(10,20,50,100,200,500,1000,2000)
sample.fractions <- c(0.01,0.05,0.1,0.2,0.5,0.8,1.0)

# E(X|Z) and E(Y|X,Z) cross-validation
CVy_on_xz <- CVm <- matrix(0,length(minimum.node.sizes),length(sample.fractions))
colnames(CVy_on_xz)=sample.fractions; rownames(CVy_on_xz)=minimum.node.sizes
colnames(CVm)=sample.fractions; rownames(CVm)=minimum.node.sizes
for (runs in 1:100) {
  rdf <- generate_data()
  K <- 2
  index <- rbinom(nrow(rdf), 1, 0.5)

  for (min.node.index in seq_len(length(minimum.node.sizes))) {
    for (sam.frac.index in seq_len(length(sample.fractions))) {
      min.node.here = minimum.node.sizes[min.node.index]
      sam.frac.here = sample.fractions[sam.frac.index]
      print(c(min.node.here,sam.frac.here))

      for (k in seq_len(K)) {
        rdf_train <- rdf[index==k-1,]
        rdf_test <- rdf[index!=k-1,]

        m_model <- ranger(D ~ X.1+X.2+X.3+X.4+X.5+X.6+X.7+X.8+X.9+X.10, data=rdf_train, min.node.size=min.node.here, sample.fraction=sam.frac.here, mtry=ceiling(sqrt(10)), num.trees = 500)
        Y_on_XZ_model <- ranger(Y ~ D+X.1+X.2+X.3+X.4+X.5+X.6+X.7+X.8+X.9+X.10, data=rdf_train, min.node.size=min.node.here, sample.fraction=sam.frac.here, mtry=ceiling(sqrt(11)), num.trees = 500)

        CVm[min.node.index,sam.frac.index] <- CVm[min.node.index,sam.frac.index] + sum((rdf_test$D-predict(m_model,rdf_test)$predictions)^2)
        CVy_on_xz[min.node.index,sam.frac.index] <- CVy_on_xz[min.node.index,sam.frac.index] + sum((rdf_test$Y-predict(Y_on_XZ_model,rdf_test)$predictions)^2)
      }
    }
  }
}
min.m = which(CVm == min(CVm), arr.ind = TRUE)
m.min.node.size = minimum.node.sizes[min.m[1]]
m.sample.fraction = sample.fractions[min.m[2]]
min.yonxz = which(CVy_on_xz == min(CVy_on_xz), arr.ind = TRUE)
yonxz.min.node.size = minimum.node.sizes[min.yonxz[1]]
yonxz.sample.fraction = sample.fractions[min.yonxz[2]]

sqrt_link_inverse = function(x) x^2
sqrt_link_deriv = function (x) 1/2 * x^(-1/2)
link_itself = function(x) sqrt(x)
# G(Z) = E[g(E[Y|X,Z])|Z] cross-validation
CVG <- matrix(0,length(minimum.node.sizes),length(sample.fractions))
colnames(CVG)=sample.fractions; rownames(CVG)=minimum.node.sizes
for (runs in 1:100) {
  rdf <- generate_data()
  K <- 2
  index <- rbinom(nrow(rdf), 1, 0.5)

  for (min.node.index in seq_len(length(minimum.node.sizes))) {
    for (sam.frac.index in seq_len(length(sample.fractions))) {
      min.node.here = minimum.node.sizes[min.node.index]
      sam.frac.here = sample.fractions[sam.frac.index]
      print(c(min.node.here,sam.frac.here))

      for (k in seq_len(K)) {
        rdf_train <- rdf[index==k-1,]
        rdf_test <- rdf[index!=k-1,]

        l_first_model <- ranger(Y ~ D+X.1+X.2+X.3+X.4+X.5+X.6+X.7+X.8+X.9+X.10, data=rdf_train, min.node.size=yonxz.min.node.size, sample.fraction=yonxz.sample.fraction, mtry=ceiling(sqrt(11)), num.trees = 500)
        G_dataset <- cbind(rdf_train, G_vals=link_itself(predict(l_first_model, rdf_train)$predictions))
        l_second_model <- ranger(G_vals ~ X.1+X.2+X.3+X.4+X.5+X.6+X.7+X.8+X.9+X.10, data=G_dataset, min.node.size=min.node.here, sample.fraction=sam.frac.here, mtry=ceiling(sqrt(10)), num.trees = 500)
        G_vals_test <- link_itself(predict(l_first_model, rdf_test)$predictions)

        fitted_G_on_z_train <- predict(l_second_model, rdf_train)$predictions
        fitted_G_on_z_test <- predict(l_second_model, rdf_test)$predictions

        CVG[min.node.index,sam.frac.index] <- CVG[min.node.index,sam.frac.index] + sum((G_vals_test-fitted_G_on_z_test)^2)
      }
    }
  }
  print(CVG-min(CVG))
}
min.G = which(CVG == min(CVG), arr.ind = TRUE)
G.min.node.size = minimum.node.sizes[min.G[1]]
G.sample.fraction = sample.fractions[min.G[2]]

# sigsq(.) cross-validation
CVsigsq <- matrix(0,length(minimum.node.sizes),length(sample.fractions))
colnames(CVsigsq)=sample.fractions; rownames(CVsigsq)=minimum.node.sizes
for (runs in 1:100) {
  rdf <- generate_data()
  K <- 2
  index <- rbinom(nrow(rdf), 1, 0.5)

  for (min.node.index in seq_len(length(minimum.node.sizes))) {
    for (sam.frac.index in seq_len(length(sample.fractions))) {
      min.node.here = minimum.node.sizes[min.node.index]
      sam.frac.here = sample.fractions[sam.frac.index]
      print(c(min.node.here,sam.frac.here))

      for (k in seq_len(K)) {
        rdf_train <- rdf[index==k-1,]
        rdf_test <- rdf[index!=k-1,]

        m_model <- ranger(D ~ X.1+X.2+X.3+X.4+X.5+X.6+X.7+X.8+X.9+X.10, data=rdf_train, min.node.size=m.min.node.size, sample.fraction=m.sample.fraction, mtry=ceiling(sqrt(10)), num.trees = 500)
        l_first_model <- ranger(Y ~ D+X.1+X.2+X.3+X.4+X.5+X.6+X.7+X.8+X.9+X.10, data=rdf_train, min.node.size=yonxz.min.node.size, sample.fraction=yonxz.sample.fraction, mtry=ceiling(sqrt(11)), num.trees = 500)
        G_dataset <- cbind(rdf_train, G_vals=link_itself(predict(l_first_model, rdf_train)$predictions))
        l_second_model <- ranger(G_vals ~ X.1+X.2+X.3+X.4+X.5+X.6+X.7+X.8+X.9+X.10, data=G_dataset, min.node.size=G.min.node.size, sample.fraction=G.sample.fraction, mtry=ceiling(sqrt(10)), num.trees = 500)
        fitted_G_on_z_train <- predict(l_second_model, rdf_train)$predictions
        fitted_G_on_z_test <- predict(l_second_model, rdf_test)$predictions

        RD_train <- rdf_train$D - predict(m_model, rdf_train)$predictions
        RD_test <- rdf_test$D - predict(m_model, rdf_test)$predictions
        eval_theta_k <- function(theta) {
          mu <- sqrt_link_inverse( RD_train*theta + fitted_G_on_z_train )
          eps <- sqrt_link_deriv(mu) * ( rdf_train$Y - mu )
          score <- sum( RD_train * eps )
          return(score)
        }
        beta_train <- uniroot(eval_theta_k,interval=c(-10,10))$root

        mu_nuisance <- sqrt_link_inverse( RD_train*beta_train + fitted_G_on_z_train )
        d_theta_eps_nuisance <- -RD_train
        eps_nuisance <- sqrt_link_deriv(mu_nuisance) * ( rdf_train$Y - mu_nuisance )
        eps_sq_nuisance <- eps_nuisance^2

        mu_test <- sqrt_link_inverse( RD_test*1 + fitted_G_on_z_test )
        d_theta_eps_nuisance <- -RD_train
        eps_test <- sqrt_link_deriv(mu_test) * ( rdf_test$Y - mu_test )
        eps_sq_test <- eps_test^2

        epsq <- eps_sq_nuisance
        xisq <- RD_train^2
        xiepsq <- RD_train^2 * eps_sq_nuisance

        w_fit <- ranger(epsq ~ D + X.1+X.2+X.3+X.4+X.5+X.6+X.7+X.8+X.9+X.10, data=cbind(epsq,rdf_train), min.node.size=min.node.here, sample.fraction=sam.frac.here, mtry=ceiling(sqrt(10)), num.trees = 500)
        w_test <- 1/predict(w_fit, rdf_test)$predictions
        w_train <- 1/predict(w_fit, rdf_train)$predictions

        CVsigsq[min.node.index,sam.frac.index] <- CVsigsq[min.node.index,sam.frac.index] + sum((eps_sq_test-predict(w_fit, rdf_test)$predictions)^2)
      }
    }
  }
}
min.sigsq = which(CVsigsq == min(CVsigsq), arr.ind = TRUE)
sig.min.node.size = minimum.node.sizes[min.sigsq[1]]
sig.sample.fraction = sample.fractions[min.sigsq[2]]

# h(.) cross-validation
CVh <- matrix(0,length(minimum.node.sizes),length(sample.fractions))
colnames(CVh)=sample.fractions; rownames(CVh)=minimum.node.sizes
for (runs in 21:100) {
  rdf <- generate_data()
  K <- 2
  index <- rbinom(nrow(rdf), 1, 0.5)

  for (min.node.index in seq_len(length(minimum.node.sizes))) {
    for (sam.frac.index in seq_len(length(sample.fractions))) {
      min.node.here = minimum.node.sizes[min.node.index]
      sam.frac.here = sample.fractions[sam.frac.index]
      print(c(min.node.here,sam.frac.here))

      for (k in seq_len(K)) {
        rdf_train <- rdf[index==k-1,]
        rdf_test <- rdf[index!=k-1,]

        m_model <- ranger(D ~ X.1+X.2+X.3+X.4+X.5+X.6+X.7+X.8+X.9+X.10, data=rdf_train, min.node.size=m.min.node.size, sample.fraction=m.sample.fraction, mtry=ceiling(sqrt(10)), num.trees = 500)
        l_first_model <- ranger(Y ~ D+X.1+X.2+X.3+X.4+X.5+X.6+X.7+X.8+X.9+X.10, data=rdf_train, min.node.size=yonxz.min.node.size, sample.fraction=yonxz.sample.fraction, mtry=ceiling(sqrt(11)), num.trees = 500)
        G_dataset <- cbind(rdf_train, G_vals=link_itself(predict(l_first_model, rdf_train)$predictions))
        l_second_model <- ranger(G_vals ~ X.1+X.2+X.3+X.4+X.5+X.6+X.7+X.8+X.9+X.10, data=G_dataset, min.node.size=G.min.node.size, sample.fraction=G.sample.fraction, mtry=ceiling(sqrt(10)), num.trees = 500)
        fitted_G_on_z_train <- predict(l_second_model, rdf_train)$predictions
        fitted_G_on_z_test <- predict(l_second_model, rdf_test)$predictions

        RD_train <- rdf_train$D - predict(m_model, rdf_train)$predictions
        RD_test <- rdf_test$D - predict(m_model, rdf_test)$predictions
        eval_theta_k <- function(theta) {
          mu <- sqrt_link_inverse( RD_train*theta + fitted_G_on_z_train )
          eps <- sqrt_link_deriv(mu) * ( rdf_train$Y - mu ) #G_link$linkderiv doesn't exist
          score <- sum( RD_train * eps )
          return(score)
        }
        beta_train <- uniroot(eval_theta_k,interval=c(-10,10))$root

        mu_nuisance <- sqrt_link_inverse( RD_train*beta_train + fitted_G_on_z_train )
        d_theta_eps_nuisance <- -RD_train
        eps_nuisance <- sqrt_link_deriv(mu_nuisance) * ( rdf_train$Y - mu_nuisance )
        eps_sq_nuisance <- eps_nuisance^2

        mu_test <- sqrt_link_inverse( RD_test*1 + fitted_G_on_z_test )
        d_theta_eps_nuisance <- -RD_train
        eps_test <- sqrt_link_deriv(mu_test) * ( rdf_test$Y - mu_test )
        eps_sq_test <- eps_test^2

        epsq <- eps_sq_nuisance
        xisq <- RD_train^2
        xiepsq <- RD_train^2 * eps_sq_nuisance

        w_fit <- ranger(epsq ~ D + X.1+X.2+X.3+X.4+X.5+X.6+X.7+X.8+X.9+X.10, data=cbind(epsq,rdf_train), min.node.size=sig.min.node.size, sample.fraction=sig.sample.fraction, mtry=ceiling(sqrt(11)), num.trees = 500)
        w_test <- 1/predict(w_fit, rdf_test)$predictions
        w_train <- 1/predict(w_fit, rdf_train)$predictions
        h <- predict(ranger(D ~ X.1+X.2+X.3+X.4+X.5+X.6+X.7+X.8+X.9+X.10, data=cbind(WD=w_train*rdf_train$D, rdf_train), min.node.size=min.node.here, sample.fraction=sam.frac.here, mtry=ceiling(sqrt(10)), num.trees = 500, case.weights=w_train), rdf_test)$predictions

        CVh[min.node.index,sam.frac.index] <- CVh[min.node.index,sam.frac.index] + sum(w_test*(rdf_test$D-h)^2)
      }
    }
  }
}
min.h = which(CVh == min(CVh), arr.ind = TRUE)
h.min.node.size = minimum.node.sizes[min.h[1]]
h.sample.fraction = sample.fractions[min.h[2]]
