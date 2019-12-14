#function that makes prediction
predict.BayesClassTree = function (tree, test_data) {

  name = NULL;
  
  for (i in 1:ncol(test_data)){
    name = c(name, paste0("X",i))
  }
  
  #set name of test_data, suppose it is standardized
  colnames(test_data) = name;
  
  #extract MAP tree
  MAPtree = tree[[1]][["MAPtree"]];
  
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