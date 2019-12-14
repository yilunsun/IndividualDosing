sz=1000 # sample size
job.id=1 # seed number when simulating data

set.seed(job.id)
x1<-runif(sz)
x2<-runif(sz)
x3<-runif(sz)
Ind_x<-NULL
Ind_x<-ifelse(x1>0.5,0.8,0.2) #0.19
Ind_x<-ifelse(x1<=0.5&x2>0.3,0.5,Ind_x)
table(Ind_x)
alpha<-1/(exp(0.3*x1+0.2*x2+x3*.1)+1)
beta<-1-alpha
d1<-rbeta(sz,beta,alpha)
reward<-(Ind_x==0.8)*dbeta(d1,24,6.75)*1.83+(Ind_x==0.5)*dbeta(d1,24,24)*1.82+(Ind_x==0.2)*dbeta(d1,6.75,24)*1.83+0.3*x1+0.4*x2+0.5*x3+rnorm(sz,sd=0.5)

quality_assessment = function(dose) {
  predictoins = sapply(dose, function(x){
    Ind_x = x
    (Ind_x==0.8)*dbeta(d1,24,6.75)*1.83+(Ind_x==0.5)*dbeta(d1,24,24)*1.82+(Ind_x==0.2)*dbeta(d1,6.75,24)*1.83+0.3*x1+0.4*x2+0.5*x3+rnorm(sz,sd=0.5)
  })
  
  return(apply(predictoins, 1, function(x){max(x) - min(x)}))
}

hist(quality_assessment(c(0.2,0.5,0.8)))

y=reward
a=d1
x=cbind(x1,x2,x3,x4=runif(sz),x5=runif(sz),x6=runif(sz))
a.out=seq(0.1,0.9,length.out=9)
#dd=sample(rep(c(0:8),each=112)[1:1000])
dd=match(Ind_x,a.out)-1

library(BART)
X1.pred = NULL
for (i in 1:9){
  X1.pred = rbind(X1.pred, cbind(i/10,x))
}
bartfit1 = wbart(cbind(a,x), y, X1.pred, sparse=T, k=qnorm(1.95/2), nskip=500, ndpost=200)$yhat.test

post_mean1 = apply(bartfit1,2,mean)
V_train_1 = matrix(post_mean1, nrow = nrow(x))
y_initial = apply(V_train_1, 1, which.max)-1
table(y_initial)
library(splines)
library(gbm)
library(BayesTreeMCMC)
library(xgboost)
library(KernSmooth)
predict.BayesClassTree = function (tree, test_data) {
  
  name = NULL;
  
  for (i in 1:ncol(test_data)){
    name = c(name, paste0("X",i))
  }
  
  #set name of test_data, suppose it is standardized
  colnames(test_data) = name;
  
  #extract MAP tree
  MAPtree = tree[[1]][["OptTree"]];
  
  #currently use majority as prediction in each leaf???
  vote = NULL;
  for (k in 1:length(MAPtree[["table"]])){
    vote = cbind(vote, which.max(MAPtree[["table"]][[k]])-1);
  }
  
  #start mimic run-down process
  #all in root
  node = rep(MAPtree[["node_name"]][[1]],nrow(test_data));
  
  #running down tree
  for (j in 1:length(MAPtree[[1]])){
    node_current = MAPtree[["node_name"]][[j]];
    if (MAPtree[["split_var"]][[j]]=="leaf") next;
    index = which(node==node_current);
    splitvar = test_data[index,MAPtree[["split_var"]][[j]]];
    splivalue = MAPtree[["split"]][[j]];
    leftson = MAPtree[["leftson"]][[j]];
    rightson = MAPtree[["rightson"]][[j]];
    if (j==1) {
      node[index] = ifelse( splitvar <= splivalue, leftson, rightson);
    } else {
      node[index] = ifelse(node[index] %in% MAPtree[["node_name"]][which(MAPtree[["split_var"]]=="leaf")], node[index],
                           ifelse( splitvar <= splivalue, leftson, rightson));
    }
  }
  
  result = NULL;
  
  for (l in 1:length(MAPtree[[1]])){
    result[which(node==MAPtree[["node_name"]][l])] = vote[l];
  }
  
  result;
}
## save a copy of original data!

#BayesianCART(NumericMatrix x, // baseline covariate
#             IntegerVector y, // categorized initial dose level: 0, 1, 2,...,8 -> 0.1, ..., 0.9; can be random or previous estimation of optimal dose level
#             NumericVector V, // V is observed reward
#             NumericVector a, // a is observed dose level: 0.1 - 0.9
#             NumericVector candidate_dose, // the dose levels to choose from
#             int cat_num, bool standardization, int burnin, int Length, int every, int nChain,
#             double size, double shape, 
#             String prior_leaf = "uniform", int MinimumLeafSize=1, unsigned int seed=123, int MinLeaf=30)
#y_initial<-sample(y_initial,1000)
y_initial2 = sample(0:8, 1000, replace = T)
tree_1 = BayesianCART(x=x, y=y_initial2, V=y,a=a, candidate_dose = a.out, 0, F, 0,10,3,1,8,0.5,"uniform",1,1123,50) 
table(predict.BayesClassTree(tree_1,test_data=x),Ind_x)
# standardization within [0,1]
# every 3 
#n chain number of mc 
# size of poisson value 
# shape  probability to grow a tree to  
#10 
#50 for minsplit 
#acceptance rate 

#line 147 
