# Cross validation for Siumulation 3 (see simulation3.R)

NUMBEROFTREES=500
MNSs = c(10, 20, 30, 40, 50, 100, 200, 500, 1000, 2000, 5000)
SPs = c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 1.0)
n=10000
NUMFORCV = 10


CV = matrix(0,length(MNSs),length(SPs));  colnames(CV)=SPs;  rownames(CV)=MNSs
for (seq_repeats in 1:NUMFORCV) {
  set.seed(seq_repeats)
  rdf = generate_data()
  rdf = cbind(rdf, Xzero = rdf$X==0)
  rdf_cv = rdf[1:(n/2),]
  rdf_cvtest = rdf[(n/2+1):n,]
  a=0; b=0
  for (MNS in MNSs) {
    a = a+1
    b=0
    for (SP in SPs) {
      b = b+1
      set.seed(seq_repeats)
      X_model <- ranger(X ~ Z.1+Z.2+Z.3, data=rdf_cv, min.node.size=MNS, sample.fraction = SP, num.trees = 500, mtry=3)
      CV[a,b] <- CV[a,b]  +  sum(  ( rdf_cvtest$X - predict(X_model, rdf_cvtest)$predictions )^2  )
    }}
  print(seq_repeats)
}
print(CV-min(CV))
X_model_MNS = MNSs[which(CV==min(CV), arr.ind=TRUE)[1]]
X_model_SP = SPs[which(CV==min(CV), arr.ind=TRUE)[2]]


# CV for E(Y|Z)
MNSs = c(10, 20, 30, 40, 50, 100, 200, 500, 1000, 2000); SPs = c(0.01, 0.05, 0.10, 0.20, 0.50, 0.80, 1.00)
CV = matrix(0,length(MNSs),length(SPs));  colnames(CV)=SPs;  rownames(CV)=MNSs
for (seq_repeats in 1:NUMFORCV) {
  print(seq_repeats)
  set.seed(seq_repeats)
  rdf = generate_data()
  rdf = cbind(rdf, Xzero = rdf$X==0)
  rdf_cv = rdf[1:(n/2),]
  rdf_cvtest = rdf[(n/2+1):n,]
  a=0; b=0
  for (MNS in MNSs) {
    a = a+1
    b=0
    for (SP in SPs) {
      b = b+1
      set.seed(seq_repeats)
      Y_model <- ranger(Y ~ Z.1+Z.2+Z.3, data=rdf_cv, min.node.size=MNS, sample.fraction = SP, num.trees = 500, mtry=3)
      CV[a,b] <- CV[a,b]  +  sum(  ( rdf_cvtest$Y - predict(Y_model, rdf_cvtest)$predictions )^2  )
    }}
}
print(CV-min(CV))
Y_model_MNS = MNSs[which(CV==min(CV), arr.ind=TRUE)[1]]
Y_model_SP = SPs[which(CV==min(CV), arr.ind=TRUE)[2]]
# CV for sigmasq(X,Z)
CV = matrix(0,length(MNSs),length(SPs));  colnames(CV)=SPs;  rownames(CV)=MNSs
for (seq_repeats in 1:NUMFORCV) {
  a=0; b=0
  print(seq_repeats)
  set.seed(seq_repeats)
  rdf = generate_data()
  rdf = cbind(rdf, Xzero = rdf$X==0)
  rdf_cv = rdf[1:(n/2),]
  rdf_cvtest = rdf[(n/2+1):n,]
  set.seed(seq_repeats)
  X_model <- ranger(X ~ Z.1+Z.2+Z.3, data=rdf_cv, min.node.size=X_model_MNS, sample.fraction = X_model_SP, num.trees = 500, mtry=3)
  set.seed(seq_repeats)
  Y_model <- ranger(Y ~ Z.1+Z.2+Z.3, data=rdf_cv, min.node.size=Y_model_MNS, sample.fraction = Y_model_SP, num.trees = 500, mtry=3)
  RX_cv = rdf_cv$X - predict(X_model, rdf_cv)$predictions
  RX_cvtest = rdf_cvtest$X - predict(X_model, rdf_cvtest)$predictions
  RY_cv = rdf_cv$Y - predict(Y_model, rdf_cv)$predictions
  RY_cvtest = rdf_cvtest$Y - predict(Y_model, rdf_cvtest)$predictions
  epssq_cv = (RY_cv - RX_cv)^2
  epssq_cvtest = (RY_cvtest - RX_cvtest)^2
  for (MNS in MNSs) {
    a = a+1
    b = 0
    for (SP in SPs) {
      b = b+1
      set.seed(seq_repeats)
      epssq_model <- ranger(epssq ~ X+Z.1+Z.2+Z.3, data=cbind(rdf_cv,epssq=epssq_cv), min.node.size=MNS, sample.fraction = SP, num.trees = 500, mtry=3)
      CV[a,b] <- CV[a,b]  +  sum(  ( epssq_cvtest - predict(epssq_model, rdf_cvtest)$predictions )^2  )
    }}
}
print(CV-min(CV))
sigsq_model_MNS = MNSs[which(CV==min(CV), arr.ind=TRUE)[1]]
sigsq_model_SP = SPs[which(CV==min(CV), arr.ind=TRUE)[2]]


# CV for h_0(Z)
CV = matrix(0,length(MNSs),length(SPs));  colnames(CV)=SPs;  rownames(CV)=MNSs
for (seq_repeats in 1:NUMFORCV) {
  a=0; b=0
  print(seq_repeats)
  set.seed(seq_repeats)
  rdf = generate_data()
  rdf = cbind(rdf, Xzero = rdf$X==0)
  rdf_cv = rdf[1:(n/2),]
  rdf_cvtest = rdf[(n/2+1):n,]
  set.seed(seq_repeats)
  X_model <- ranger(X ~ Z.1+Z.2+Z.3, data=rdf_cv, min.node.size=X_model_MNS, sample.fraction = X_model_SP, num.trees = 500, mtry=3)
  set.seed(seq_repeats)
  Y_model <- ranger(Y ~ Z.1+Z.2+Z.3, data=rdf_cv, min.node.size=Y_model_MNS, sample.fraction = Y_model_SP, num.trees = 500, mtry=3)
  RX_cv = rdf_cv$X - predict(X_model, rdf_cv)$predictions
  RX_cvtest = rdf_cvtest$X - predict(X_model, rdf_cvtest)$predictions
  RY_cv = rdf_cv$Y - predict(Y_model, rdf_cv)$predictions
  RY_cvtest = rdf_cvtest$Y - predict(Y_model, rdf_cvtest)$predictions
  epssq_cv = (RY_cv - RX_cv)^2
  epssq_cvtest = (RY_cvtest - RX_cvtest)^2
  set.seed(seq_repeats)
  epssq_model <- ranger(epssq ~ X+Z.1+Z.2+Z.3, data=cbind(rdf_cv,epssq=epssq_cv), min.node.size=sigsq_model_MNS, sample.fraction = sigsq_model_SP, num.trees = 500, mtry=3)
  
  for (MNS in MNSs) {
    a = a+1
    b = 0
    for (SP in SPs) {
      b=b+1
      set.seed(seq_repeats)
      h_model <- ranger(X ~ Z.1+Z.2+Z.3, data=rdf_cv, min.node.size=MNS, sample.fraction = SP, num.trees = 500, case.weights = 1/predict(epssq_model, rdf_cv)$predictions, mtry=3)
      
      CV[a,b] <- CV[a,b]  +  sum(  1/predict(epssq_model, rdf_cvtest)$predictions * ( rdf_cvtest$X - predict(h_model, rdf_cvtest)$predictions )^2  )
    }}
  print(CV-min(CV))
}
print(CV-min(CV))
h_model_MNS = MNSs[which(CV==min(CV), arr.ind=TRUE)[1]]
h_model_SP = SPs[which(CV==min(CV), arr.ind=TRUE)[2]]


# CV for P(X=0|Z)
MNSs = c(10, 20, 30, 40, 50, 100, 200, 500, 1000, 2000); SPs = c(0.01, 0.05, 0.10, 0.20, 0.50)
CV = matrix(0,length(MNSs),length(SPs));  colnames(CV)=SPs;  rownames(CV)=MNSs
for (seq_repeats in 1:NUMFORCV) {
  set.seed(seq_repeats)
  rdf = generate_data()
  rdf = cbind(rdf, Xzero = as.numeric(rdf$X==0))
  rdf_cv = rdf[1:(n/2),]
  rdf_cvtest = rdf[(n/2+1):n,]
  a=0; b=0
  for (MNS in MNSs) {
    a = a+1
    b=0
    for (SP in SPs) {
      b = b+1
      set.seed(seq_repeats)
      Xzero_model <- grf::probability_forest(X=rdf_cv[,1:3], Y=as.factor(rdf_cv[,"Xzero"]), num.trees=500, sample.fraction=SP, mtry=3, honesty=FALSE, min.node.size=MNS)
      CV[a,b] <- CV[a,b]  +  sum(  ( rdf_cvtest$Xzero - predict(Xzero_model, rdf_cvtest[,1:3])$predictions[,2] )^2  )
    }}
  print(seq_repeats)
}
print(CV-min(CV))
Xzero_model_MNS = MNSs[which(CV==min(CV), arr.ind=TRUE)[1]]
Xzero_model_SP = SPs[which(CV==min(CV), arr.ind=TRUE)[2]]
