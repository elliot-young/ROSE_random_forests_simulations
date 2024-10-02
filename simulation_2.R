# Reproducible code for Simulation 1a
# Note that one can fit a rose random forest estimator directly using [https://github.com/elliot-young/rose].

# Load relevant rose random forest functions
{
  {
    library(rpart)
    library(ranger)
    library(mgcv)
    library(glmnet)

    itemp <- function(y, offset, parms, wt) {
      if (ncol(y) != 2) {
        stop("Matrix of response must be a 2 column matrix")
      }
      if (!missing(parms) && length(parms) > 0){
        warning("parameter argument ignored")
      }
      if (length(offset)) y <- y - offset
      sfun <- function(weigh, avar, ylevel, digits) {
        paste(" weights=", format(signif(weigh, digits)), ", AsymVar=", format(signif(avar, digits)), sep='')
      }
      environment(sfun) <- .GlobalEnv
      list(y = y, parms = NULL, numresp = 1, numy = 2, summary=sfun)
    }
    etemp <- function(y, wt, parms) {
      sum_xisq = sum(y[,1])
      sum_xiepsq = sum(y[,2])
      weigh = sum_xisq/sum_xiepsq
      avar = length(y[,1]) - weigh*sum_xisq
      list(label = weigh, deviance = avar)
    }
    stemp <- function(y, wt, x, parms, continuous) {
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
    ulist <- list(eval = etemp, split = stemp, init = itemp)

    # For Rsq loss:
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
      rss <- sum(wt*(y-wmean)^2)
      list(label = wmean, deviance = rss)
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

    # RANDOM FOREST ALL
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
      avar = length(y[,1]) - sum_xisq*sum_xisq/sum_xiepsq
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
  {
    library(rpart)
    library(ranger)
    library(mgcv)
    library(glmnet)

    # RANDOM FOREST ALL
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
      avar = 1000*length(y[,1]) - sum_xisq*sum_xisq/sum_xiepsq
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

}

SEED.HERE=1
T=1000
NUMTREES=100;
VARS_ALL_eachn <- BETAS_ALL_eachn <- CVS_ALL_eachn <- list()
maxdepths = 1:15

n = 12800
{
  ours_vs_unw = Rsq_vs_unw = ours_vs_Rsq = numeric(0)
  BETAS_ALL = VARS_ALL = CVS_ALL = list()
  A = 0
  for (maxdepth in maxdepths) {
    A=A+1
    MAXDEPTH=maxdepth
    cat("MAXDEPTH = ",maxdepth)
    BETAS = VARS = CVS = data.frame()
    b_Rsq = numeric(0)
    b_ours = numeric(0)
    {
      MAXDEPTH=MAXDEPTH; SUBSIZEPROP = 1; MINBUCKET=10
      allbetas_wrong <- allbetas_unw <- allbetas_oursRsq <- allbetas_ours <- numeric(0)
      allbetas_ours1=allbetas_ours2=allbetas_ours3=allbetas_ours4 <- numeric(0)
      allbetas_oursRsq1=allbetas_oursRsq2=allbetas_oursRsq3=allbetas_oursRsq4 <- numeric(0)
      vars_wrong <- vars_unw <- vars_oursRsq <- vars_ours <- numeric(0)
      vars_ours1=vars_ours2=vars_ours3=vars_ours4 <- numeric(0)
      vars_oursRsq1=vars_oursRsq2=vars_oursRsq3=vars_oursRsq4 <- numeric(0)
    }

    set.seed(SEED.HERE)
    while (dim(BETAS)[1]<T) {
      print_progress(dim(BETAS)[1], T)
      CV1=CV2=CV3=CV4=0
      CV0 = 0

      allbetas_ours=allbetas_oursRsq=vars_ours=vars_oursRsq <- numeric(0)

      beta=1

      f_errY = function(D,X1,X2,X3) 1+2*expit(1*(X1+X2)) + 2*expit(X2+X3)

      X1 = rnorm(n,0,4)
      X2 = rnorm(n,0,4)
      X3 = rnorm(n,0,4)
      D = rnorm(n,0,4)
      err_Y = f_errY(D,X1,X2,X3)*rnorm(n)
      Y = beta*D + 1*err_Y

      rdf <- data.frame(X1=X1, X2=X2, X3=X3, D=D, Y=Y)

      K <- 1
      index <- rbinom(nrow(rdf), 1, 0.5)

      # One that'll go wrong
      beta_hat_final_num_wrong <- beta_hat_final_den_wrong <- numeric(K)
      beta_hat_final_num_unw <- beta_hat_final_den_unw <- numeric(K)
      beta_hat_final_num_ours <- beta_hat_final_den_ours <- numeric(K)
      beta_hat_final_num_oursRsq <- beta_hat_final_den_oursRsq <- numeric(K)
      beta_hat_final_num_oracle <- beta_hat_final_den_oracle <- numeric(K)

      beta_hat_final_num_ours1 <- beta_hat_final_den_ours1 <- numeric(K)
      beta_hat_final_num_oursRsq1 <- beta_hat_final_den_oursRsq1 <- numeric(K)
      beta_hat_final_num_ours2 <- beta_hat_final_den_ours2 <- numeric(K)
      beta_hat_final_num_oursRsq2 <- beta_hat_final_den_oursRsq2 <- numeric(K)
      beta_hat_final_num_ours3 <- beta_hat_final_den_ours3 <- numeric(K)
      beta_hat_final_num_oursRsq3 <- beta_hat_final_den_oursRsq3 <- numeric(K)
      beta_hat_final_num_ours4 <- beta_hat_final_den_ours4 <- numeric(K)
      beta_hat_final_num_oursRsq4 <- beta_hat_final_den_oursRsq4 <- numeric(K)


      V_num_wrong <- V_num_unw <- V_num_ours <- V_num_oursRsq <- numeric(K)
      V_num_ours1 <- V_num_oursRsq1 <- numeric(K)
      V_num_ours2 <- V_num_oursRsq2 <- numeric(K)
      V_num_ours3 <- V_num_oursRsq3 <- numeric(K)
      V_num_ours4 <- V_num_oursRsq4 <- numeric(K)

      for (k in seq_len(K)) {
        rdf_train <- rdf[index==k-1,]
        rdf_test <- rdf[index!=k-1,]


        RD_train <- rdf_train$D
        RD_test <- rdf_test$D
        RY_train <- rdf_train$Y
        RY_test <- rdf_test$Y

        beta_train <- sum(RY_train*RD_train)/sum(RD_train*RD_train)

        epsq <- (RY_train-beta_train*RD_train)^2#; colnames(epsq) <- "epsq"
        xisq <- RD_train*RD_train#; colnames(xisq) <- "xisq"
        xiepsq <- RD_train^2*(RY_train-beta_train*RD_train)^2#; colnames(xiepsq) <- "xiepsq"




        MAXDEPTH=MAXDEPTH
        w_test_randomforest_NUM <- w_test_randomforest_DEN <- matrix(0,dim(rdf_test)[1],NUMTREES)
        for (fsts in seq_len(NUMTREES)) {
          DTdata <- cbind(xisq=xisq,xiepsq=xiepsq,rdf_train)
          DTdata <- DTdata[sample(seq_len(dim(rdf_train)[1]), floor(SUBSIZEPROP*dim(rdf_train)[1]), replace=T),]#bootstrap
          w_fit <- rpart(cbind(xisq, xiepsq) ~ X1+X2+X3, data = DTdata, method = ulist_rf, cp=0, minbucket=MINBUCKET, maxdepth=MAXDEPTH)
          w_test_randomforest_NUM[,fsts] <- predict(w_fit, rdf_test, type="matrix")[,1]
          w_test_randomforest_DEN[,fsts] <- predict(w_fit, rdf_test, type="matrix")[,2]
        }
        w_test <- rowMeans(w_test_randomforest_NUM)/rowMeans(w_test_randomforest_DEN)
        beta_hat_final_num_ours3[k] <- sum(w_test*RY_test*RD_test)
        beta_hat_final_den_ours3[k] <- sum(w_test*RD_test*RD_test)
        V_num_ours3[k] <- sum(w_test^2*(RY_test-beta_train*RD_test)^2*RD_test^2)
        #4



        # OTHER Rsq - naive forest
        #1
        MAXDEPTH=MAXDEPTH
        w_test_randomforest <- matrix(0,dim(rdf_test)[1],NUMTREES)
        for (fsts in seq_len(NUMTREES)) {
          DTdata <- cbind(xisq=xisq,xiepsq=xiepsq,rdf_train)
          DTdata <- DTdata[sample(seq_len(dim(rdf_train)[1]), floor(SUBSIZEPROP*dim(rdf_train)[1]), replace=T),]#bootstrap
          rsq_fit_test <- rpart(epsq ~ X1+X2+X3, data = cbind(epsq=epsq,rdf_train), method = rsqulist, cp=0, minbucket=MINBUCKET, maxdepth=MAXDEPTH) #minsplit = sqrt(n)
          w_test_randomforest[,fsts] <- predict(rsq_fit_test, rdf_test)
        }
        w_test <- 1/rowMeans(w_test_randomforest)
        preds <- rowMeans(w_test_randomforest)
        CV1=CV1+sum((preds-(RY_test-beta_train*RD_test)^2)^2)
        beta_hat_final_num_oursRsq1[k] <- sum(w_test*RY_test*RD_test)
        beta_hat_final_den_oursRsq1[k] <- sum(w_test*RD_test*RD_test)
        V_num_oursRsq1[k] <- sum(w_test^2*(RY_test-beta_train*RD_test)^2*RD_test^2)

        preds0 = mean(epsq)
        CV0 = CV0 + sum((preds0-(RY_test-beta_train*RD_test)^2)^2)

        # unweighted
        beta_hat_final_num_unw[k] <- sum(RY_test*RD_test)
        beta_hat_final_den_unw[k] <- sum(RD_test*RD_test)
        V_num_unw[k] <- sum((RY_test-beta_train*RD_test)^2*RD_test^2)

        # oracle
        beta_hat_final_num_wrong[k] <- sum(1/(f_errY(rdf_test$D,rdf_test$X1,rdf_test$X2,rdf_test$X3)^2)*RY_test*RD_test)
        beta_hat_final_den_wrong[k] <- sum(1/(f_errY(rdf_test$D,rdf_test$X1,rdf_test$X2,rdf_test$X3)^2)*RD_test*RD_test)
        V_num_wrong[k] <- sum(1/(f_errY(rdf_test$D,rdf_test$X1,rdf_test$X2,rdf_test$X3)^4)*(RY_test-beta_train*RD_test)^2*RD_test^2)
      }

      #oracle
      beta_final <- sum(beta_hat_final_num_wrong)/sum(beta_hat_final_den_wrong)
      allbetas_wrong <- rep(beta_final,4)
      vars_wrong <- rep( dim(rdf)[1]*sum(V_num_wrong)/(sum(beta_hat_final_den_wrong))^2,  4 )
      #
      beta_final <- sum(beta_hat_final_num_unw)/sum(beta_hat_final_den_unw)
      allbetas_unw <- rep(beta_final,4)
      vars_unw <- rep( dim(rdf)[1]*sum(V_num_unw)/(sum(beta_hat_final_den_unw))^2,  4 )
      #
      beta_final <- sum(beta_hat_final_num_oursRsq1)/sum(beta_hat_final_den_oursRsq1)
      allbetas_oursRsq <- append(allbetas_oursRsq, beta_final)
      vars_oursRsq <- append(vars_oursRsq, dim(rdf)[1]*sum(V_num_oursRsq1)/(sum(beta_hat_final_den_oursRsq1))^2)
      #
      beta_final <- sum(beta_hat_final_num_ours3)/sum(beta_hat_final_den_ours3)
      allbetas_ours <- append(allbetas_ours, beta_final)
      vars_ours <- append(vars_ours, dim(rdf)[1]*sum(V_num_ours3)/(sum(beta_hat_final_den_ours3))^2)

      BETAS = rbind(BETAS, data.frame(UNW=allbetas_unw[1], ours3=allbetas_ours[1], ourRsq1=allbetas_oursRsq[1], ORACLE=allbetas_wrong[1] ))
      VARS = rbind(VARS, data.frame(UNW=vars_unw[1], ours3=vars_ours[1], ourRsq1=vars_oursRsq[1], ORACLE=vars_wrong[1] ))
      CVS = rbind(CVS, data.frame(ourRsq0=CV0, ourRsq1=CV1))


      aaa = 1:(dim(VARS)[1])
      b_Rsq = append(b_Rsq, 100*(1-mean(VARS$ourRsq1)/mean(VARS$UNW)))
      b_ours = append(b_ours, 100*(1-mean(VARS$ours3)/mean(VARS$UNW)))
    }
    c(100-100*mean(VARS$ours3)/mean(VARS$UNW),     100-100*mean(VARS$ourRsq1)/mean(VARS$UNW),    100-100*mean(VARS$ours3)/mean(VARS$ourRsq1))
    ours_vs_unw=append(ours_vs_unw,   100-100*mean(VARS$ours3)/mean(VARS$UNW))
    Rsq_vs_unw=append(Rsq_vs_unw,    100-100*mean(VARS$ourRsq1)/mean(VARS$UNW))
    ours_vs_Rsq=append(ours_vs_Rsq,    100-100*mean(VARS$ours3)/mean(VARS$ourRsq1))
    VARS_ALL[[A]] = VARS
    BETAS_ALL[[A]] = BETAS
    CVS_ALL[[A]] = CVS
    print(paste0("COMPLETED n = ",n,", MAXDEPTH=", maxdepth))
  }
  {
    ourimprovements = theirimprovements = oraimprovements = numeric(length(maxdepths)); A=0
    for (mmm in maxdepths) {
      A=A+1
      ourimprovements[A] = mean(VARS_ALL[[A]]$ours3)/mean(VARS_ALL[[1]]$UNW)
      theirimprovements[A] = mean(VARS_ALL[[A]]$ourRsq1)/mean(VARS_ALL[[1]]$UNW)
      oraimprovements[A] = mean(VARS_ALL[[A]]$ORACLE)/mean(VARS_ALL[[1]]$UNW)
    }
    plot(maxdepths, ourimprovements, col="dark orange", ylim=c(min(ourimprovements,theirimprovements),max(ourimprovements,theirimprovements)))
    points(maxdepths, theirimprovements, col="dark green")
    points(maxdepths, oraimprovements, col="magenta")
    points(maxdepths,rep(1,length(maxdepths)))
  }
}
VARS_ALL_eachn$n_12800 = VARS_ALL
BETAS_ALL_eachn$n_12800 = BETAS_ALL
CVS_ALL_eachn$n_12800 = CVS_ALL
n = 25600
{
  ours_vs_unw = Rsq_vs_unw = ours_vs_Rsq = numeric(0)
  BETAS_ALL = VARS_ALL = CVS_ALL = list()
  A = 0
  for (maxdepth in maxdepths) {
    A=A+1
    MAXDEPTH=maxdepth
    cat("MAXDEPTH = ",maxdepth)
    BETAS = VARS = CVS = data.frame()
    b_Rsq = numeric(0)
    b_ours = numeric(0)
    {
      MAXDEPTH=MAXDEPTH; SUBSIZEPROP = 1; MINBUCKET=10
      allbetas_wrong <- allbetas_unw <- allbetas_oursRsq <- allbetas_ours <- numeric(0)
      allbetas_ours1=allbetas_ours2=allbetas_ours3=allbetas_ours4 <- numeric(0)
      allbetas_oursRsq1=allbetas_oursRsq2=allbetas_oursRsq3=allbetas_oursRsq4 <- numeric(0)
      vars_wrong <- vars_unw <- vars_oursRsq <- vars_ours <- numeric(0)
      vars_ours1=vars_ours2=vars_ours3=vars_ours4 <- numeric(0)
      vars_oursRsq1=vars_oursRsq2=vars_oursRsq3=vars_oursRsq4 <- numeric(0)
    }

    set.seed(SEED.HERE)
    while (dim(BETAS)[1]<T) {
      print_progress(dim(BETAS)[1], T)
      CV1=CV2=CV3=CV4=0
      CV0 = 0

      allbetas_ours=allbetas_oursRsq=vars_ours=vars_oursRsq <- numeric(0)

      beta=1

      f_errY = function(D,X1,X2,X3) 1+2*expit(1*(X1+X2)) + 2*expit(X2+X3)

      X1 = rnorm(n,0,4)
      X2 = rnorm(n,0,4)
      X3 = rnorm(n,0,4)
      D = rnorm(n,0,4)
      err_Y = f_errY(D,X1,X2,X3)*rnorm(n)
      Y = beta*D + 1*err_Y

      rdf <- data.frame(X1=X1, X2=X2, X3=X3, D=D, Y=Y)

      K <- 1
      index <- rbinom(nrow(rdf), 1, 0.5)

      # One that'll go wrong
      beta_hat_final_num_wrong <- beta_hat_final_den_wrong <- numeric(K)
      beta_hat_final_num_unw <- beta_hat_final_den_unw <- numeric(K)
      beta_hat_final_num_ours <- beta_hat_final_den_ours <- numeric(K)
      beta_hat_final_num_oursRsq <- beta_hat_final_den_oursRsq <- numeric(K)
      beta_hat_final_num_oracle <- beta_hat_final_den_oracle <- numeric(K)

      beta_hat_final_num_ours1 <- beta_hat_final_den_ours1 <- numeric(K)
      beta_hat_final_num_oursRsq1 <- beta_hat_final_den_oursRsq1 <- numeric(K)
      beta_hat_final_num_ours2 <- beta_hat_final_den_ours2 <- numeric(K)
      beta_hat_final_num_oursRsq2 <- beta_hat_final_den_oursRsq2 <- numeric(K)
      beta_hat_final_num_ours3 <- beta_hat_final_den_ours3 <- numeric(K)
      beta_hat_final_num_oursRsq3 <- beta_hat_final_den_oursRsq3 <- numeric(K)
      beta_hat_final_num_ours4 <- beta_hat_final_den_ours4 <- numeric(K)
      beta_hat_final_num_oursRsq4 <- beta_hat_final_den_oursRsq4 <- numeric(K)


      V_num_wrong <- V_num_unw <- V_num_ours <- V_num_oursRsq <- numeric(K)
      V_num_ours1 <- V_num_oursRsq1 <- numeric(K)
      V_num_ours2 <- V_num_oursRsq2 <- numeric(K)
      V_num_ours3 <- V_num_oursRsq3 <- numeric(K)
      V_num_ours4 <- V_num_oursRsq4 <- numeric(K)

      for (k in seq_len(K)) {
        rdf_train <- rdf[index==k-1,]
        rdf_test <- rdf[index!=k-1,]


        RD_train <- rdf_train$D
        RD_test <- rdf_test$D
        RY_train <- rdf_train$Y
        RY_test <- rdf_test$Y

        beta_train <- sum(RY_train*RD_train)/sum(RD_train*RD_train)

        epsq <- (RY_train-beta_train*RD_train)^2#; colnames(epsq) <- "epsq"
        xisq <- RD_train*RD_train#; colnames(xisq) <- "xisq"
        xiepsq <- RD_train^2*(RY_train-beta_train*RD_train)^2#; colnames(xiepsq) <- "xiepsq"




        MAXDEPTH=MAXDEPTH
        w_test_randomforest_NUM <- w_test_randomforest_DEN <- matrix(0,dim(rdf_test)[1],NUMTREES)
        for (fsts in seq_len(NUMTREES)) {
          DTdata <- cbind(xisq=xisq,xiepsq=xiepsq,rdf_train)
          DTdata <- DTdata[sample(seq_len(dim(rdf_train)[1]), floor(SUBSIZEPROP*dim(rdf_train)[1]), replace=T),]#bootstrap
          w_fit <- rpart(cbind(xisq, xiepsq) ~ X1+X2+X3, data = DTdata, method = ulist_rf, cp=0, minbucket=MINBUCKET, maxdepth=MAXDEPTH) #maxdepth=MAXDEPTH
          w_test_randomforest_NUM[,fsts] <- predict(w_fit, rdf_test, type="matrix")[,1]
          w_test_randomforest_DEN[,fsts] <- predict(w_fit, rdf_test, type="matrix")[,2]
        }
        w_test <- rowMeans(w_test_randomforest_NUM)/rowMeans(w_test_randomforest_DEN)
        beta_hat_final_num_ours3[k] <- sum(w_test*RY_test*RD_test)
        beta_hat_final_den_ours3[k] <- sum(w_test*RD_test*RD_test)
        V_num_ours3[k] <- sum(w_test^2*(RY_test-beta_train*RD_test)^2*RD_test^2)
        #4



        # OTHER Rsq - naive forest
        #1
        MAXDEPTH=MAXDEPTH
        w_test_randomforest <- matrix(0,dim(rdf_test)[1],NUMTREES)
        for (fsts in seq_len(NUMTREES)) {
          DTdata <- cbind(xisq=xisq,xiepsq=xiepsq,rdf_train)
          DTdata <- DTdata[sample(seq_len(dim(rdf_train)[1]), floor(SUBSIZEPROP*dim(rdf_train)[1]), replace=T),]#bootstrap
          rsq_fit_test <- rpart(epsq ~ X1+X2+X3, data = cbind(epsq=epsq,rdf_train), method = rsqulist, cp=0, minbucket=MINBUCKET, maxdepth=MAXDEPTH) #minsplit = sqrt(n)
          w_test_randomforest[,fsts] <- predict(rsq_fit_test, rdf_test)
        }
        w_test <- 1/rowMeans(w_test_randomforest)
        preds <- rowMeans(w_test_randomforest)
        CV1=CV1+sum((preds-(RY_test-beta_train*RD_test)^2)^2)
        beta_hat_final_num_oursRsq1[k] <- sum(w_test*RY_test*RD_test)
        beta_hat_final_den_oursRsq1[k] <- sum(w_test*RD_test*RD_test)
        V_num_oursRsq1[k] <- sum(w_test^2*(RY_test-beta_train*RD_test)^2*RD_test^2)

        preds0 = mean(epsq)
        CV0 = CV0 + sum((preds0-(RY_test-beta_train*RD_test)^2)^2)

        # unweighted
        beta_hat_final_num_unw[k] <- sum(RY_test*RD_test)
        beta_hat_final_den_unw[k] <- sum(RD_test*RD_test)
        V_num_unw[k] <- sum((RY_test-beta_train*RD_test)^2*RD_test^2)

        # oracle
        beta_hat_final_num_wrong[k] <- sum(1/(f_errY(rdf_test$D,rdf_test$X1,rdf_test$X2,rdf_test$X3)^2)*RY_test*RD_test)
        beta_hat_final_den_wrong[k] <- sum(1/(f_errY(rdf_test$D,rdf_test$X1,rdf_test$X2,rdf_test$X3)^2)*RD_test*RD_test)
        V_num_wrong[k] <- sum(1/(f_errY(rdf_test$D,rdf_test$X1,rdf_test$X2,rdf_test$X3)^4)*(RY_test-beta_train*RD_test)^2*RD_test^2)
      }

      #oracle
      beta_final <- sum(beta_hat_final_num_wrong)/sum(beta_hat_final_den_wrong)
      allbetas_wrong <- rep(beta_final,4)
      vars_wrong <- rep( dim(rdf)[1]*sum(V_num_wrong)/(sum(beta_hat_final_den_wrong))^2,  4 )
      #
      beta_final <- sum(beta_hat_final_num_unw)/sum(beta_hat_final_den_unw)
      allbetas_unw <- rep(beta_final,4)
      vars_unw <- rep( dim(rdf)[1]*sum(V_num_unw)/(sum(beta_hat_final_den_unw))^2,  4 )
      #
      beta_final <- sum(beta_hat_final_num_oursRsq1)/sum(beta_hat_final_den_oursRsq1)
      allbetas_oursRsq <- append(allbetas_oursRsq, beta_final)
      vars_oursRsq <- append(vars_oursRsq, dim(rdf)[1]*sum(V_num_oursRsq1)/(sum(beta_hat_final_den_oursRsq1))^2)
      #
      beta_final <- sum(beta_hat_final_num_ours3)/sum(beta_hat_final_den_ours3)
      allbetas_ours <- append(allbetas_ours, beta_final)
      vars_ours <- append(vars_ours, dim(rdf)[1]*sum(V_num_ours3)/(sum(beta_hat_final_den_ours3))^2)

      BETAS = rbind(BETAS, data.frame(UNW=allbetas_unw[1], ours3=allbetas_ours[1], ourRsq1=allbetas_oursRsq[1], ORACLE=allbetas_wrong[1] ))
      VARS = rbind(VARS, data.frame(UNW=vars_unw[1], ours3=vars_ours[1], ourRsq1=vars_oursRsq[1], ORACLE=vars_wrong[1] ))
      CVS = rbind(CVS, data.frame(ourRsq0=CV0, ourRsq1=CV1))

      aaa = 1:(dim(VARS)[1])
      b_Rsq = append(b_Rsq, 100*(1-mean(VARS$ourRsq1)/mean(VARS$UNW)))
      b_ours = append(b_ours, 100*(1-mean(VARS$ours3)/mean(VARS$UNW)))
    }
    c(100-100*mean(VARS$ours3)/mean(VARS$UNW),     100-100*mean(VARS$ourRsq1)/mean(VARS$UNW),    100-100*mean(VARS$ours3)/mean(VARS$ourRsq1))
    ours_vs_unw=append(ours_vs_unw,   100-100*mean(VARS$ours3)/mean(VARS$UNW))
    Rsq_vs_unw=append(Rsq_vs_unw,    100-100*mean(VARS$ourRsq1)/mean(VARS$UNW))
    ours_vs_Rsq=append(ours_vs_Rsq,    100-100*mean(VARS$ours3)/mean(VARS$ourRsq1))
    VARS_ALL[[A]] = VARS
    BETAS_ALL[[A]] = BETAS
    CVS_ALL[[A]] = CVS
    print(paste0("COMPLETED n = ",n,", MAXDEPTH=", maxdepth))
  }
  {
    ourimprovements = theirimprovements = oraimprovements = numeric(length(maxdepths)); A=0
    for (mmm in maxdepths) {
      A=A+1
      ourimprovements[A] = mean(VARS_ALL[[A]]$ours3)/mean(VARS_ALL[[1]]$UNW)
      theirimprovements[A] = mean(VARS_ALL[[A]]$ourRsq1)/mean(VARS_ALL[[1]]$UNW)
      oraimprovements[A] = mean(VARS_ALL[[A]]$ORACLE)/mean(VARS_ALL[[1]]$UNW)
    }
    plot(maxdepths, ourimprovements, col="dark orange", ylim=c(min(ourimprovements,theirimprovements),max(ourimprovements,theirimprovements)))
    points(maxdepths, theirimprovements, col="dark green")
    points(maxdepths, oraimprovements, col="magenta")
    points(maxdepths,rep(1,length(maxdepths)))
  }
}
VARS_ALL_eachn$n_25600 = VARS_ALL
BETAS_ALL_eachn$n_25600 = BETAS_ALL
CVS_ALL_eachn$n_25600 = CVS_ALL
