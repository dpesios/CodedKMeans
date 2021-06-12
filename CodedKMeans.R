
library(igraph)
library(Matrix)
library(lattice)
library(vcd)
library(brotli)

# -- AUXILIARY FUNCTIONS -- 
findClusterIdentifierInGi <- function(u, TT) {
  # stopifnot(class(TT)=="matrix")
  
  identifier_cl <- NA;
  
  stop_flag <- FALSE;
  for(i in 1:as.integer(dim(TT)[[1]])) {
    for(j in 1:as.integer(dim(TT)[[2]])) {
      if(u %in% TT[[i, j]]) {
        identifier_cl <- i;
        stop_flag <- TRUE;
        break
      }
    }
    
    if (stop_flag) { break }
  }
  
  identifier_cl
}

findConnectedCompIdentifierInGi <- function(u, TT) {
  # stopifnot(class(TT)=="matrix")
  
  identifier_cc <- NA;
  
  stop_flag <- FALSE;
  for(i in 1:as.integer(dim(TT)[[1]])) {
    for(j in 1:as.integer(dim(TT)[[2]])) {
      if(u %in% TT[[i, j]]) {
        identifier_cc <- j;
        stop_flag <- TRUE;
        break
      }
    }
    
    if (stop_flag) { break }
  }
  
  identifier_cc
}

getPossibleWorldOfUGraph <- function(uGraph, randomNum) {
  possibleWorldGraph <- delete.edges(
    uGraph, 
    E(uGraph)[E(uGraph)$Probabilities <= randomNum]
  );
  
  possibleWorldGraph
}

# Calculate probability of a random possible world (subgraph)
calculateProbabilityOfPW <- function(edgeListInitial, edgeListPW) {
  term <- 1
  # for(i in 1:length(edgeListPW)) {
  #   term <- ifelse(edgeListPW[i] %in% edgeListInitial, 
  #                  edgeListInitial[i]$Probabilities, 
  #                  1 - edgeListInitial[i]$Probabilities
  #   )
  #   term <- term*term
  # }
  for(i in 1:length(edgeListInitial)) {
    term <- ifelse(edgeListInitial[i] %in% edgeListInitial, 
                   edgeListInitial[i]$Probabilities, 
                   1 - edgeListInitial[i]$Probabilities
    )
    
    term <- as.numeric(term)
    term <- term*term
  }
  term
}

# Calculates the current iteration value of the objective function F_s
calculateObjFunctionFs <- function(T_i, K, initialUGraph, possibleWorld) {
  pr_i <- calculateProbabilityOfPW(E(initialUGraph), E(possibleWorld));
  
  productSum <- 0;
  
  for(k in 1:K) {
    # print(norm_C_k(k, T_i));
    # print(H_ij_k(k, T_i, gorder(initialUGraph)));
    
    productSum <- productSum + norm_C_k(k, T_i)*H_ij_k(k, T_i, gorder(initialUGraph));
  }
  
  pr_i*productSum
}

# Counts the number of cluster members by means of auxiliary table T_i
norm_C_k <- function(row_identifier, T_i) {
  number_of_members <- 0;
  
  for(c in 1:dim(T_i)[[2]]) {
    number_of_members <- number_of_members + length(T_i[[row_identifier, c]]);
  }
  
  number_of_members
}

# Calculates the entropy of a fragment (CC) of a possible world
H_ij_k <- function(row_identifier, T_i, max_elements) {
  sum_quantity <- 0.0;
  
  for(r in 1:dim(T_i)[[1]]) {
    numerators_list <- vector(mode = "list", length = max_elements);
    denominator <- 0;
    
    for(c in 1:dim(T_i)[[2]]) {
      numerators_list[[c]] <- append(numerators_list[[c]], c(length(T_i[[r,c]])))
      denominator <- denominator + length(T_i[[r,c]])
    }
    
    for(i in 1:length(numerators_list)) {
      if(!is.null(numerators_list[[i]])) {
        sum_quantity <- sum_quantity - (ifelse(numerators_list[[i]] == 0, 0, (numerators_list[[i]]/denominator)*log10(numerators_list[[i]]/denominator)));
      } else {
        break;
      }
    }
    
    numerators_list <- NULL;
  }
  
  sum_quantity
}
# END OF AUXILIARY FUNCTIONS

# probGraphDataDF <- read.csv('toy_example.txt', sep = "\t", header=FALSE);
probGraphDataDF <- read.csv('nature04670-s7.csv', sep = ",", header=FALSE);
colnames(probGraphDataDF)[3] <- c("Probabilities");
(probGraphDataDF);

# Create graph from data frame
pbGraph <- graph.data.frame(probGraphDataDF, directed = FALSE);

# plot(pbGraph, layout = layout_nicely(pbGraph),
#      vertex.color = "grey",
#      vertex.label = V(pbGraph)$names,
#      edge.label = E(pbGraph)$Probabilities,
#      edge.label.font = 1,
#      vertex.label.font = 33,
#      vertex.size = 15,
#      edge.arrow.size = 0.1,
#      edge.lty = 'solid',
#      edge.width = 1,
#      edge.curved = FALSE);

# Declare the number of clusters to return
K <- 6;
# Declare the number of samples and the max number of iterations
N <- 150; max_iterations <- 100;

# Create a possible world instantiation of the uncertain graph
randomNum <- runif(1);
possibleWorldGraph <- getPossibleWorldOfUGraph(pbGraph, randomNum);
  
pwGraphFragments <- components(possibleWorldGraph);
num_of_fragments <- pwGraphFragments$no;
  
# Create an auxialiary table for current instantiation
T_i <- matrix(rep(list(), K*num_of_fragments), nrow = K, ncol = num_of_fragments)

for(u in 1:length(pwGraphFragments$membership)) {
  column_identifier <- pwGraphFragments$membership[[u]] # get connected component identifier
  row_identifier <- sample.int(K, 1, replace = FALSE); # randomly choose the cluster to assign to
    
  T_i[[row_identifier, column_identifier]] <- append(T_i[[row_identifier, column_identifier]], c(names(pwGraphFragments$membership)[u]))
}

F_value_next <- 0;
F_value_prev <- calculateObjFunctionFs(T_i, K, pbGraph, possibleWorldGraph);

tolerance <- 0.1;
iterator <- 0;

finalClusterConfiguration <- vector(mode = "list", length = K);

while( ((F_value_prev - F_value_next) > tolerance) && (iterator < max_iterations)) {
  F_value_next <- F_value_prev;
  
  # STAGE 1
  keepTiMatrices <- list(); # list of T_i matrices
  keepTiMatricesCoded <- list(); # list of T_i matrices coded
  possbileWorldList <- list();
  
  for(i in 1:N) {
    possibleWorldGraph <- getPossibleWorldOfUGraph(pbGraph, runif(1));
    possbileWorldList[[i]] <- possibleWorldGraph;
    
    pwGraphFragments <- components(possibleWorldGraph);
    num_of_fragments <- pwGraphFragments$no;
    
    # Create an auxiliary table for current instantiation
    T_i <- matrix(rep(list(), K*num_of_fragments), nrow = K, ncol = num_of_fragments) # matrix of lists
    TiMatrixEncoded <- matrix(rep(list(), K*num_of_fragments), nrow = K, ncol = num_of_fragments) # matrix of encodings for G_i
    
    for(u in 1:length(pwGraphFragments$membership)) {
      column_identifier <- pwGraphFragments$membership[[u]] # get connected component identifier
      row_identifier <- sample.int(K, 1, replace = FALSE); # randomly choose the cluster to assign to
      
      T_i[[row_identifier, column_identifier]] <- append(T_i[[row_identifier, column_identifier]], c(names(pwGraphFragments$membership)[u]));
    }
    
    keepTiMatrices[[i]] <- T_i; # possible world G_i lies in i position

    for(k in 1:dim(T_i)[[1]]) {
      for(u in 1:dim(T_i)[[2]]) {
        TiMatrixEncoded[[k, u]] <- brotli_compress(charToRaw(toString(T_i[[k, u]]))); # instead of pure huffman encoding
      }
    }
    
    keepTiMatricesCoded[[i]] <- TiMatrixEncoded; 
  }
  
  # STAGE 2
  F_value_prev <- 0.0;
  
  for(x in 1:length(V(pbGraph))) {
    u <- names(V(pbGraph))[x];
    
    codelengths_list <- vector(mode = "numeric", length = K);
    sum_codeLength_of_u <- 0;
    
    for(i in 1:N) {
      T_i_mat <- keepTiMatrices[[i]];
      T_i_coded_mat <- keepTiMatricesCoded[[i]]; # same index as above

      cc_identifier_of_u <- findConnectedCompIdentifierInGi(u, T_i_mat); # dim(T_i_mat)[[2]]
      cluster_identifier_of_u <- findClusterIdentifierInGi(u, T_i_mat);
      
      for(k in 1:K) {
        sum_codeLength_of_u <- sum_codeLength_of_u + length(T_i_coded_mat[[k, cc_identifier_of_u]]);
      }
      
      codelengths_list[[cluster_identifier_of_u]] <- codelengths_list[[cluster_identifier_of_u]] + sum_codeLength_of_u;
      
      F_value_prev <- F_value_prev + calculateObjFunctionFs(T_i_mat, K, pbGraph, possbileWorldList[[i]]);
    }
    
    # Assign u to cluster where code length sum is minimized
    index <- which.min(codelengths_list);
    
    finalClusterConfiguration[[index]] <- append(finalClusterConfiguration[[index]], u);
  }
  print(F_value_prev)
  
  iterator <- iterator + 1;
}

print("Final Cluster Configuration:");
print(finalClusterConfiguration)

capture.output(
  cat(format(finalClusterConfiguration), sep = "\n"), 
  file = paste("CConfiguration", toString(K), format(Sys.time(), "%d-%b-%Y %H.%M"), ".txt", sep = "_"),
  append = FALSE
)
