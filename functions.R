cond_order_m <- function(X, Y, Z, sample_size, orientation = "in", nredo = 500, kernel_type = "epanechnikov",
                         nmulti=5, ftol=1.19e-07){
  # Conditional Order-m Efficiency Estimation
  #
  # This function computes conditional order-m efficiency scores (Cazals et al., 2002;Daraio and Simar, 2005)
  #
  #   PARAMETERS
  #   X: (n, p) matrix of inputs with n units and p input variables
  #   Y: (n, q) matrix of outputs with n units and q output variables
  #   Z: (n, r) matrix of contextual with n units and r contextual variables
  #   sample_size: Number of observations to benchmark in each Monte Carlo iteration 
  #   orientation:  orientation: "in" for input orientation, "out" for output orientation ("in" by default)
  #   nredo: Number of Monte Carlo replications (500 by default)
  #   kernel_type: Type of kernel (only continuous variables) to use for density estimation ("epanechnikov" by default)
  #   nmulti: Number of multi-starts for bandwidth selection (5 by default)
  #   ftol: Tolerance level for convergence in bandwidth selection (1.19e-07 by default)
  #
  #   RETURNS
  #   A data frame containing:
  #   eff: Efficiency scores for each observation
  #   error: Standard errors of the efficiency estimates
  #   band_contin: Bandwidths used for each observation and environmental variables

  
  kernel_value <- function(x, kernel = "epanechnikov"){
    # return the value of different kernel functions
    
    if (!kernel %in% c("gaussian", "epanechnikov", "biweight", "triangular"))
      stop("Incorrect name of kernel function")
    
    if (kernel == "gaussian")
      return(dnorm(x, mean = 0, sd = 1))
    if (kernel == "epanechnikov")
      return(ifelse(abs(x) <= 1, 3/4 * (1 - x^2), 0))
    if (kernel == "biweight")
      return(ifelse(abs(x) <= 1, 15/16 * (1 - x^2)^2, 0))
    if (kernel == "triangular")
      return(ifelse(abs(x) <= 1, 1 - abs(x), 0))
  }
  
  filter_of_sample <- function(matrix, number, orient = ">="){
    
    t <- c()
    if (orient == ">="){
      
      for (i in 1:nrow(matrix)){
        if (FALSE %in% (matrix[i, ] >= matrix[number, ]))
          t <- c(t, i)
      }
    } else if (orient == "<="){
      
      for (i in 1:nrow(matrix)){
        if (FALSE %in% (matrix[i, ] <= matrix[number, ]))
          t <- c(t, i)
      }
    } 
    else stop("Error! Please, try again.")
    return (t) # return the numbers of observations that then will be excluded
  }
  
  
  # Input validation checks
  if(is.matrix(X) == FALSE | is.matrix(Y) == FALSE | is.matrix(Z) == FALSE )
    stop("Error! Incorrect type of variables")
  if(nrow(X) != nrow(Y) | nrow(X) != nrow(Z)) 
    stop("All matrices must have the same number of rows")
  if(!orientation %in% c("in", "out"))
    stop("orientation must be either 'in' or 'out'")
  
  # Initialize result vectors
  eff_score <- c()       # To store efficiency scores
  error_score <- c()     # To store standard errors  
  band_contin <- c()     # To store bandwidths used
  
  # Input orientation case
  if (orientation == "in" & nrow(X) == nrow(Y) ){ 
    for (i in 1:nrow(X)){
      
      # Get reference set where Y_i >= y (dominating observations)
      t <- filter_of_sample(Y, i) 
      
      # Handle different cases for reference set size
      if (is.null(t)){
        X_sub <- X
        Z_sub <- Z
      } 
      else if (length(t) == nrow(X) - 1){
        X_sub <- t(X[-t, ])
        Z_sub <- t(Z[-t, ])
      }
      else{
        X_sub <- as.matrix(X[-t, ])
        Z_sub <- as.matrix(Z[-t, ])
      }
      
      # Only proceed if we have enough observations
      if (nrow(Z_sub) > 1){
        # Estimate bandwidths for kernel density
        kerz <- np::npcdensbw(xdat=Z_sub, 
                              ydat=X_sub, 
                              bwmethod="cv.ls", 
                              ckertype = "epanechnikov", 
                              ukertype = "liracine",
                              nmulti=nmulti, 
                              ftol=ftol)
        
        # Calculate adaptive bandwidth
        h <- kerz$xbw*nrow(X)**(-ncol(Y)/((4+ncol(Y)+ncol(Z))*(4+ncol(Z) )) )
        
        eff_score_j <- c() # Temporary storage for Monte Carlo results
        prob_kern <- c()   # Kernel probabilities
        
        # Calculate kernel weights for reference set
        for (k in 1:nrow(X_sub)){
          z_sub_k <- 1
          for (l in 1:ncol(Z_sub))
            z_sub_k <- z_sub_k * (1/h[l])**(ncol(Z)) * prod(kernel_value((Z[i, l] - Z_sub[k, l]) / h[l], kernel = kernel_type))
          prob_kern <- c(prob_kern, z_sub_k)
        }
        
        # Calculate kernel weights for full sample (for normalization)
        kernel <- c()
        for (k in 1:nrow(Y)){
          z_k <- 1
          for (l in 1:ncol(Z))
            z_k <- z_k * (1/h[l])**(ncol(Z)) * prod(kernel_value((Z[i, l] - Z[k, l]) / h[l], kernel = kernel_type))
          kernel <- c(kernel, z_k)
            }
         }
      else{ # Case when not enough observations
        eff_score <- c(eff_score, 1) # Default efficiency score: no units to compare
        error_score <- c(error_score, NA)
        band_contin <- rbind(band_contin, rep(NA, ncol(Z)) )
        next
      } 
      
      # Normalize probabilities
      prob_kern <- prob_kern / sum(kernel)
      # Monte Carlo simulation
      for (j in 1:nredo){
        # Resample with probabilities based on kernel weights
        X_rand <- X_sub[sample(nrow(X_sub), sample_size, replace = TRUE, prob = prob_kern), ]
        # Calculate efficiency score for this replication
        eff_score_j <- c(eff_score_j, min( matrixStats::rowMaxs(t(t(X_rand) / X[i, ])) ))
      }
      
      # Store results for this observation
      eff_score <- c(eff_score, mean(eff_score_j)) 
      error_score <- c(error_score, sd(eff_score_j) / sqrt(nredo))
      band_contin <- rbind(band_contin, h)
      cat(i ," ", eff_score[i]," ", h, "\n") # print to follow the procedure
    }
  } 
  # Output orientation case
  else if (orientation == "out" & nrow(X) == nrow(Y) ){ # start of output orientation
    for (i in 1:nrow(Y)){
      
      # Get reference set where X_i <= x (dominated observations)
      t <- filter_of_sample(X, i, orient = "<=")
      
      # Handle different cases for reference set size
      if (is.null(t)){
        Y_sub <- Y
        Z_sub <- Z
      } 
      else if (length(t) == nrow(Y) - 1){
        Y_sub <- t(Y[-t, ])
        Z_sub <- t(Z[-t, ])
      }
      else{
        Y_sub <- as.matrix(Y[-t, ])
        Z_sub <- as.matrix(Z[-t, ])
      }
      
      # Only proceed if we have enough observations
      if (nrow(Z_sub) > 1){
        
        # Estimate bandwidths for kernel density
        kerz <- np::npcdensbw(xdat=Z_sub, 
                              ydat=Y_sub, 
                              bwmethod="cv.ls", 
                              ckertype = "epanechnikov", 
                              ukertype = "liracine",
                              nmulti=nmulti,
                              ftol=ftol)               
        
        # Calculate adaptive bandwidth
        h <- kerz$xbw*nrow(X)**(-ncol(Y)/((4+ncol(Y)+ncol(Z))*(4+ncol(Z) )) )
        
        eff_score_j <- c() # Temporary storage for Monte Carlo results
        prob_kern <- c()   # Kernel probabilities
        
        # Calculate kernel weights for reference set
        for (k in 1:nrow(Y_sub)){
          z_sub_k <- 1
          for (l in 1:ncol(Z_sub))
            z_sub_k <- z_sub_k * (1/h[l])**(ncol(Z)) * prod(kernel_value((Z[i, l] - Z_sub[k, l]) / h[l], kernel = kernel_type))
          prob_kern <- c(prob_kern, z_sub_k)
        }
        
        # Calculate kernel weights for full sample (for normalization)
        kernel <- c()
        for (k in 1:nrow(Y)){
          z_k <- 1
          for (l in 1:ncol(Z))
            z_k <- z_k * (1/h[l])**(ncol(Z)) * prod(kernel_value((Z[i, l] - Z[k, l]) / h[l], kernel = kernel_type))
          kernel <- c(kernel, z_k)
        }
      }
      else{ # Case when not enough observations
        eff_score <- c(eff_score, 1) # Default efficiency score: no units to compare
        error_score <- c(error_score, NA)
        band_contin <- rbind(band_contin, rep(NA, ncol(Z)) )
        next
      }
      
      # Normalize probabilities
      prob_kern <- prob_kern / sum(kernel)
      
      # Monte Carlo simulation
      for (j in 1:nredo){
        # Resample with probabilities based on kernel weights
        Y_rand <- Y_sub[sample(nrow(Y_sub), sample_size, replace = TRUE, prob = prob_kern), ]
        # Calculate efficiency score for this replication
        eff_score_j <- c(eff_score_j, max( matrixStats::rowMins(t(t(Y_rand) / Y[i, ])) ))
      }
      # Store results for this observation
      eff_score <- c(eff_score, mean(eff_score_j)) 
      error_score <- c(error_score, sd(eff_score_j) / sqrt(nredo))
      band_contin <- rbind(band_contin, h)
      cat(i ," ", eff_score[i]," ", h, "\n")
    }
  } 
  
  # Return results as a data frame
  return(data.frame(eff = eff_score, error = error_score, band_contin=band_contin))
}


order_m <- function(X, Y, sample_size, orientation = "in", nredo = 2000){
  # Order-m Efficiency Estimation
  #
  # This function computes order-m efficiency scores (Cazals et al., 2002)
  #
  #   PARAMETERS
  #   X: (n, p) matrix of inputs with n units and p input variables
  #   Y: (n, q) matrix of outputs with n units and q output variables
  #   sample_size: Number of observations to benchmark in each Monte Carlo iteration 
  #   orientation:  orientation: "in" for input orientation, "out" for output orientation ("in" by default)
  #   nredo: Number of Monte Carlo replications (500 by default)
  #
  #   RETURNS
  #   A data frame containing:
  #   eff: Efficiency scores for each observation
  #   error: Standard errors of the efficiency estimates

  
  eff_score <- c()
  error_score <- c()

  filter_of_sample <- function(matrix, number, orient = ">="){
    
    t <- c()
    if (orient == ">="){
      
      for (i in 1:nrow(matrix)){
        if (FALSE %in% (matrix[i, ] >= matrix[number, ]))
          t <- c(t, i)
      }
    } else if (orient == "<="){
      
      for (i in 1:nrow(matrix)){
        if (FALSE %in% (matrix[i, ] <= matrix[number, ]))
          t <- c(t, i)
      }
    } 
    else stop("Error! Please, try again.")
    return (t) # return the numbers of observations that then will be excluded
  }
  
  if (orientation == "in" & nrow(X) == nrow(Y)){
    for (i in 1:nrow(X)){
      
      t <- filter_of_sample(Y, i) # inlude input observations such that Y_i >= y
      if (is.null(t)) 
        X_sub <- X
      else if (length(t) == nrow(X) - 1) 
        X_sub <- t(X[-t, ])
      else 
        X_sub <- as.matrix(X[-t, ])
      
      eff_score_j <- c()
      for (j in 1:nredo){
        
        X_rand <- X_sub[sample(nrow(X_sub), sample_size, replace = TRUE), ]
        eff_score_j <- c(eff_score_j, min( matrixStats::rowMaxs(t(t(X_rand) / X[i, ])) ))
      }
      eff_score <- c(eff_score, mean(eff_score_j))
      error_score <- c(error_score, sd(eff_score_j) / sqrt(nredo))
    }
  }
  else if (orientation == "out" & nrow(X) == nrow(Y)){
    for (i in 1:nrow(Y)){
      
      t <- filter_of_sample(X, i, orient = "<=") # inlude input observations such that X_i <= x
      if (is.null(t)) 
        Y_sub <- Y
      else if (length(t) == nrow(Y) - 1) 
        Y_sub <- t(Y[-t, ])
      else 
        Y_sub <- as.matrix(Y[-t, ])
      
      eff_score_j <- c()
      for (j in 1:nredo){
        
        Y_rand <- Y_sub[sample(nrow(Y_sub), sample_size, replace = TRUE), ]
        eff_score_j <- c(eff_score_j, max( matrixStats::rowMins(t(t(Y_rand) / Y[i, ])) ))
      }
      eff_score <- c(eff_score, mean(eff_score_j))
      error_score <- c(error_score, sd(eff_score_j) / sqrt(nredo))
    }
  }
  else stop("Error! Please, try again.")
  return(data.frame(eff = eff_score, error = error_score))
}




ROC_curve <- function(X, Y, grid, Z = NULL, h = NULL, orientation = "in", nredo = 1000){
  
  # Calculate ROC curve showing the percentage of points outside order-m frontier (both conditional and unconditional)
  # as a function of m (the number of objects against which unit is benchmarked)
  # The approach is described in (Daraio and Simar, 2005)
  #
  # PARAMETERS
  # X: (n, p) matrix with n units and p input variables
  # Y: (n, q) matrix with n units and q output variables 
  # grid: vector of m values to be estimated 
  # Z: (n, r) matrix with n units and r exogenous variables for conditional order-m procedure. Defaults to NULL
  # h: r-dimensional vector of bandwidthes of exogenous factors for conditional order-m procedure. Defaults to NULL
  # orientation: model orientation. Can be "in" (input orientation, default) or "out (output orientation)"
  #
  # RETURNS
  # Grid plot where m values are on x-axis and the percentages of units outside order-m frontier (depend on m) are on y-axis
  
  require(matrixStats)
  
  #if (!is.matrix(X) | !is.matrix(Y) | !is.vector(grid))
  #  stop("X and Y must be matrixes and grid must be vector")
  #if (!is.null(Z) & !is.matrix(Z))
  #  stop("Z must be matrix")
  #if (!orientation %in% c("in", "out"))
  #  stop("Orientation must be 'in'(input) or 'out' (output)")
  
  percent <- c()
  if (is.null(Z)){ # unconditional order-m 
    if (orientation == "in"){
      for (i in 1:length(grid))
        percent <- c(percent, length(which(order_m(X,Y,sample_size=grid[i],orientation=orientation,nredo=nredo)$eff > 1)) / nrow(X))
    }
    if (orientation == "out"){
      for (i in 1:length(grid))
        percent <- c(percent, length(which(order_m(X,Y,sample_size=grid[i],orientation=orientation,nredo=nredo)$eff < 1)) / nrow(X))
    }
  } else{ # conditional order-m
    if (orientation == "in"){
      for (i in 1:length(grid))
        percent <- c(percent, length(which(cond_order_m(X,Y,Z,sample_size=grid[i],orientation=orientation,nredo=nredo)$eff > 1)) / nrow(X))
    }
    if (orientation == "out"){
      for (i in 1:length(grid))
        percent <- c(percent, length(which(cond_order_m(X,Y,Z,sample_size=grid[i],orientation=orientation,nredo=nredo)$eff < 1)) / nrow(X))
    }
  }
  
  return(plot(grid, percent, xlab = "m", ylab = "Percentage of units outside frontier"))  
}
