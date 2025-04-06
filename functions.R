


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

filter_of_sample_fdh <- function(matrix, number, Z, h, orient = ">="){
  
  t <- c()
  if (orient == ">="){
    
    for (i in 1:nrow(matrix)){
      if ((FALSE %in% (matrix[i, ] >= matrix[number, ])) & 
          (FALSE %in% (abs(Z[i, ] - Z[number, ]) <= h)) )
        t <- c(t, i)
    }
  } else if (orient == "<="){
    
    for (i in 1:nrow(matrix)){
      if ((FALSE %in% (matrix[i, ] <= matrix[number, ])) &
          (FALSE %in% (abs(Z[i, ] - Z[number, ]) <= h)) )
        t <- c(t, i)
    }
  }
  else stop("Error! Please, try again.")
  return (t) # return the numbers of observations that then will be excluded
}


fdh <- function(X, Y, orientation = "in"){
  
  if(is.matrix(X) == FALSE | is.matrix(Y) == FALSE)
    stop("Error! Incorrect type of variables")
  
  eff_score <- c()
  if (orientation == "in" & nrow(X) == nrow(Y)){
    for (i in 1:nrow(X)){
      
      t <- filter_of_sample(Y, i) # inlude input observations such that Y_i >= y
      if (is.null(t)) 
        X_sub <- X
      else if (length(t) == nrow(X) - 1) 
        X_sub <- t(X[-t, ])
      else 
        X_sub <- as.matrix(X[-t, ])
      
      eff_score <- c(eff_score, min( rowMaxs(t(t(X_sub) / X[i, ])) ))
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
      
      eff_score <- c(eff_score, max( rowMins(t(t(Y_sub) / Y[i, ])) ))
    }
  }
  else stop("Error! Please, try again.")
  return(eff_score)
}


order_m <- function(X, Y, sample_size, orientation = "in", nredo = 1000){
  
 # if(is.matrix(X) == FALSE | is.matrix(Y) == FALSE)
  #  stop("Error! Incorrect type of variables")
  
  require(matrixStats)
  
  eff_score <- c()
  error_score <- c()
  
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
        eff_score_j <- c(eff_score_j, min( rowMaxs(t(t(X_rand) / X[i, ])) ))
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
        eff_score_j <- c(eff_score_j, max( rowMins(t(t(Y_rand) / Y[i, ])) ))
      }
      eff_score <- c(eff_score, mean(eff_score_j))
      error_score <- c(error_score, sd(eff_score_j) / sqrt(nredo))
    }
  }
  else stop("Error! Please, try again.")
  return(data.frame(eff = eff_score, error = error_score))
}


cond_fdh <- function(X, Y, Z, h, orientation = "in"){
  
  eff_score <- c()
  if (orientation == "in" & nrow(X) == nrow(Y)){
    for (i in 1:nrow(X)){
      
      t <- filter_of_sample_fdh(Y, i, Z, h)
      if (is.null(t))
        X_sub <- X
      else if (length(t) == nrow(X) - 1)
        X_sub <- t(X[-t, ])
      else
        X_sub <- as.matrix(X[-t, ])
      
      eff_score <- c(eff_score, min( rowMaxs(t(t(X_sub) / X[i, ])) ))
    }
  }
  else if (orientation == "out" & nrow(X) == nrow(Y)){
    for (i in 1:nrow(Y)){
      
      t <- filter_of_sample_fdh(X, i, Z, h, orient = "<=")
      if (is.null(t))
        Y_sub <- Y
      else if (length(t) == nrow(X) - 1)
        Y_sub <- t(Y[-t, ])
      else
        Y_sub <- as.matrix(Y[-t, ])
      
      eff_score <- c(eff_score, max( rowMins(t(t(Y_sub) / Y[i, ])) ))
    }
  }
  else stop("Error! Please, try again.")
  return(eff_score)
}


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

ordkernel <- function(h, z_k, z_i){
  # ordered Li and Racine (2007) discrete kernel 
  return(h**(abs(z_k - z_i)))
}


cond_order_m <- function(X, Y, Z, h, sample_size, orientation = "in", nredo = 500, kernel_type = "epanechnikov"){
  
  if(is.matrix(X) == FALSE | is.matrix(Y) == FALSE | is.matrix(Z) == FALSE | is.vector(h) == FALSE)
    stop("Error! Incorrect type of variables")
  
  require(matrixStats)
  
  eff_score <- c()
  error_score <- c()
  if (orientation == "in" & nrow(X) == nrow(Y) & ncol(Z) == length(h)){
    for (i in 1:nrow(X)){
      
      t <- filter_of_sample(Y, i) # include input observations such that Y_i >= y
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
      
      eff_score_j <- c()
      prob_kern <- c()
      for (k in 1:nrow(X_sub)){
        z_sub_k <- 1
        for (l in 1:ncol(Z_sub))
          z_sub_k <- z_sub_k * kernel_value((Z[i, l] - Z_sub[k, l]) / h[l], kernel = kernel_type)
        prob_kern <- c(prob_kern, z_sub_k)
      }
      
      kernel <- c()
      for (k in 1:nrow(X)){
        z_k <- 1
        for (l in 1:ncol(Z))
          z_k <- z_k * kernel_value((Z[i, l] - Z[k, l]) / h[l], kernel = kernel_type)
        kernel <- c(kernel, z_k)
      }
      
      prob_kern <- prob_kern / sum(kernel)
      for (j in 1:nredo){
        
        X_rand <- X_sub[sample(nrow(X_sub), sample_size, replace = TRUE, prob = prob_kern), ]
        eff_score_j <- c(eff_score_j, min( rowMaxs(t(t(X_rand) / X[i, ])) ))
      }
      eff_score <- c(eff_score, mean(eff_score_j)) 
      error_score <- c(error_score, sd(eff_score_j) / sqrt(nredo))
    }
  }
  else if (orientation == "out" & nrow(X) == nrow(Y) & ncol(Z) == length(h)){
    for (i in 1:nrow(Y)){
      
      t <- filter_of_sample(X, i, orient = "<=") # include input observations such that X_i <= x
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
      
      eff_score_j <- c()
      prob_kern <- c()
      for (k in 1:nrow(Y_sub)){
        z_sub_k <- 1
        for (l in 1:ncol(Z_sub))
          z_sub_k <- z_sub_k * kernel_value((Z[i, l] - Z_sub[k, l]) / h[l], kernel = kernel_type)
        prob_kern <- c(prob_kern, z_sub_k)
      }
      #cat(str(Y_sub), ' ', i, '\n')
      kernel <- c()
      for (k in 1:nrow(Y)){
        z_k <- 1
        for (l in 1:ncol(Z))
          z_k <- z_k * kernel_value((Z[i, l] - Z[k, l]) / h[l], kernel = kernel_type)
        kernel <- c(kernel, z_k)
      }
      
      prob_kern <- prob_kern / sum(kernel)
      for (j in 1:nredo){
        
        Y_rand <- Y_sub[sample(nrow(Y_sub), sample_size, replace = TRUE, prob = prob_kern), ]
        eff_score_j <- c(eff_score_j, max( rowMins(t(t(Y_rand) / Y[i, ])) ))
      }
      eff_score <- c(eff_score, mean(eff_score_j)) 
      error_score <- c(error_score, sd(eff_score_j) / sqrt(nredo))
    }
  }
  else stop("Error! The number of rows of matrixes must be equal")
  
  return(data.frame(eff = eff_score, error = error_score))
}


cond_order_m2 <- function(X, Y, Z, sample_size, orientation = "in", nredo = 500, kernel_type = "epanechnikov"){
  
  #if(is.matrix(X) == FALSE | is.matrix(Y) == FALSE | is.matrix(Z) == FALSE | is.vector(h) == FALSE)
  #  stop("Error! Incorrect type of variables")
  
  require(matrixStats)
  
  eff_score <- c()
  error_score <- c()
  band_contin <- c()
  if (orientation == "in" & nrow(X) == nrow(Y) ){
    for (i in 1:nrow(X)){
      
      t <- filter_of_sample(Y, i) # include input observations such that Y_i >= y
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
      
      eff_score_j <- c()
      prob_kern <- c()
      for (k in 1:nrow(X_sub)){
        z_sub_k <- 1
        for (l in 1:ncol(Z_sub))
          z_sub_k <- z_sub_k * kernel_value((Z[i, l] - Z_sub[k, l]) / h[l], kernel = kernel_type)
        prob_kern <- c(prob_kern, z_sub_k)
      }
      
      kernel <- c()
      for (k in 1:nrow(X)){
        z_k <- 1
        for (l in 1:ncol(Z))
          z_k <- z_k * kernel_value((Z[i, l] - Z[k, l]) / h[l], kernel = kernel_type)
        kernel <- c(kernel, z_k)
      }
      
      prob_kern <- prob_kern / sum(kernel)
      for (j in 1:nredo){
        
        X_rand <- X_sub[sample(nrow(X_sub), sample_size, replace = TRUE, prob = prob_kern), ]
        eff_score_j <- c(eff_score_j, min( rowMaxs(t(t(X_rand) / X[i, ])) ))
      }
      eff_score <- c(eff_score, mean(eff_score_j)) 
      error_score <- c(error_score, sd(eff_score_j) / sqrt(nredo))
    }
  }
  else if (orientation == "out" & nrow(X) == nrow(Y) ){
    for (i in 1:nrow(Y)){
      
      t <- filter_of_sample(X, i, orient = "<=") # include input observations such that X_i <= x
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
      
      if (nrow(Z_sub) > 1){
        #kerz <- np::npcdensbw(xdat=Z_sub, ydat=Y_sub, bwmethod="cv.ls", ckertype = "epanechnikov", ukertype = "liracine",
        #                             remin=FALSE, nmulti=0, ftol=0.1, tol=0.1)
        kerz <- np::npcdensbw(xdat=Z_sub, ydat=Y_sub, bwmethod="cv.ls", ckertype = "epanechnikov", ukertype = "liracine",
                              nmulti=3)                                     
        
        h <- kerz$xbw
        h <- h*nrow(X)**(-ncol(Y)/(4+ncol(Y)+ncol(Z)*(4+ncol(Z))))
        
        eff_score_j <- c()
        prob_kern <- c()
        
        for (k in 1:nrow(Y_sub)){
          z_sub_k <- 1
          z_sub_k <- z_sub_k * (1/h)**(ncol(Z)) * prod(kernel_value((Z[i, ] - Z_sub[k, ]) / h, kernel = kernel_type))
          prob_kern <- c(prob_kern, z_sub_k)
        }
        
        kernel <- c()
        for (k in 1:nrow(Y)){
          z_k <- 1
          z_k <- z_k * (1/h)**(ncol(Z)) * prod(kernel_value((Z[i, ] - Z[k, ]) / h, kernel = kernel_type))
          kernel <- c(kernel, z_k)
        }
      }
      else{
        eff_score <- c(eff_score, 1)
        error_score <- c(error_score, NA)
        band_contin <- c(band_contin, NA)
        next
      }
      
      prob_kern <- prob_kern / sum(kernel)
      for (j in 1:nredo){
        
        Y_rand <- Y_sub[sample(nrow(Y_sub), sample_size, replace = TRUE, prob = prob_kern), ]
        eff_score_j <- c(eff_score_j, max( rowMins(t(t(Y_rand) / Y[i, ])) ))
      }
      eff_score <- c(eff_score, mean(eff_score_j)) 
      error_score <- c(error_score, sd(eff_score_j) / sqrt(nredo))
      band_contin <- c(band_contin, h)
      cat(i ," ", eff_score[i]," ", h, "\n")
    }
  }
  else stop("Error! The number of rows of matrixes must be equal")
  
  return(data.frame(eff = eff_score, error = error_score, bandwidth = band_contin))
}


ROC_curve <- function(X, Y, grid, Z = NULL, h = NULL, orientation = "in"){
  
  # Calculate ROC curve showing the percentage of points outside order-m frontier (both conditional and unconditional)
  # as a function of m (the number of objects against which unit is benchmarked)
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
        percent <- c(percent, length(which(order_m(X,Y,sample_size=grid[i],orientation=orientation)$eff > 1)) / nrow(X))
    }
    if (orientation == "out"){
      for (i in 1:length(grid))
        percent <- c(percent, length(which(order_m(X,Y,sample_size=grid[i],orientation=orientation)$eff < 1)) / nrow(X))
    }
  } else{ # conditional order-m
    if (orientation == "in"){
      for (i in 1:length(grid))
        percent <- c(percent, length(which(cond_order_m2(X,Y,Z,sample_size=grid[i],orientation=orientation)$eff > 1)) / nrow(X))
    }
    if (orientation == "out"){
      for (i in 1:length(grid))
        percent <- c(percent, length(which(cond_order_m2(X,Y,Z,sample_size=grid[i],orientation=orientation)$eff < 1)) / nrow(X))
    }
  }
  
  #return(plot(grid, percent, xlab = "m", ylab = "Percentage of units outside frontier"))  
  return(percent)
}

cond_order_m_v2 <- function(X, Y, Z, sample_size, orientation = "in", nredo = 200, kernel_type = "epanechnikov"){
  
  #if(is.matrix(X) == FALSE | is.matrix(Y) == FALSE | is.matrix(Z) == FALSE | is.matrix(h) == FALSE)
  #  stop("Error! Incorrect type of variables")
  
  require(matrixStats)
  
  n_contin_z <- sum(sapply(Z, is.numeric))
  n_ordered_z <- sum(sapply(Z, is.ordered))
  n_unordered_z <- ncol(Z) - n_contin_z - n_ordered_z
  cat(n_contin_z, " ", n_ordered_z, "\n")
  
  Z <- as.matrix(Z)
  Z <- apply(Z, 2, as.numeric)
  
  eff_score <- c()
  error_score <- c()
  band_contin <- c()
  band_ordered <- c()
  if (orientation == "in" & nrow(X) == nrow(Y)){
    for (i in 1:nrow(X)){
      
      t <- filter_of_sample(Y, i) # include input observations such that Y_i >= y
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
      
      kerz <- npudens(bws=npudensbw(dat=Z_sub, bwmethod="cv.ls"),
                      cykertype="epanechnikov",cxkertype="epanechnikov",
                      tdat=Z_sub, edat=Z_sub)
      h <- kerz$bw
      
      eff_score_j <- c()
      prob_kern <- c()
      for (k in 1:nrow(X_sub)){
        z_sub_k <- 1
        for (l in 1:ncol(Z_sub))
          z_sub_k <- z_sub_k * kernel_value((Z[i, l] - Z_sub[k, l]) / h[l], kernel = kernel_type)
        prob_kern <- c(prob_kern, z_sub_k)
      }
      
      kernel <- c()
      for (k in 1:nrow(X)){
        z_k <- 1
        for (l in 1:ncol(Z))
          z_k <- z_k * kernel_value((Z[i, l] - Z[k, l]) / h[l], kernel = kernel_type)
        kernel <- c(kernel, z_k)
      }
      
      prob_kern <- prob_kern / sum(kernel)
      for (j in 1:nredo){
        
        X_rand <- X_sub[sample(nrow(X_sub), sample_size, replace = TRUE, prob = prob_kern), ]
        eff_score_j <- c(eff_score_j, min( rowMaxs(t(t(X_rand) / X[i, ])) ))
      }
      eff_score <- c(eff_score, mean(eff_score_j)) 
      error_score <- c(error_score, sd(eff_score_j) / sqrt(nredo))
      band <- c(band, h)
    }
  }
  else if (orientation == "out" & nrow(X) == nrow(Y)){
    for (i in 1:nrow(Y)){
      
      t <- filter_of_sample(X, i, orient = "<=") # include input observations such that X_i <= x
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
      
      if (nrow(Z_sub) > 1){
        #kerz <- np::npcdensbw(xdat=Z_sub, ydat=Y_sub, bwmethod="cv.ls", ckertype = "epanechnikov", ukertype = "liracine",
        #                             remin=FALSE, nmulti=0, ftol=0.1, tol=0.1)
        kerz <- np::npcdensbw(xdat=Z_sub, ydat=Y_sub, bwmethod="cv.ls", ckertype = "epanechnikov", ukertype = "liracine",
                              nmulti=2)                                     
                          
        h <- kerz$xbw
        h[1:n_contin_z] <- h[1:n_contin_z]*nrow(X)**(-ncol(Y)/(4+ncol(Y)+ncol(Z)*(4+ncol(Z))))
      
        eff_score_j <- c()
        prob_kern <- c()
        
        for (k in 1:nrow(Y_sub)){
          z_sub_k <- 1
          z_sub_k <- z_sub_k * (1/h[1:n_contin_z])**(n_contin_z) * prod(kernel_value((Z[i, 1:n_contin_z] - Z_sub[k, 1:n_contin_z]) / h[1:n_contin_z], kernel = kernel_type))
          z_sub_k <- z_sub_k * prod(ordkernel(h[(n_contin_z+1):(n_contin_z+n_ordered_z)], Z[i, (n_contin_z+1):(n_contin_z+n_ordered_z)], 
                                                                      Z_sub[k, (n_contin_z+1):(n_contin_z+n_ordered_z)]) )
          prob_kern <- c(prob_kern, z_sub_k)
        }
        
        kernel <- c()
        for (k in 1:nrow(Y)){
          z_k <- 1
          z_k <- z_k * (1/h[1:n_contin_z])**(n_contin_z) * prod(kernel_value((Z[i, 1:n_contin_z] - Z[k, 1:n_contin_z]) / h[1:n_contin_z], kernel = kernel_type))
          z_k <- z_k * prod(ordkernel(h[(n_contin_z+1):(n_contin_z+n_ordered_z)], Z[i, (n_contin_z+1):(n_contin_z+n_ordered_z)], 
                                                              Z[k, (n_contin_z+1):(n_contin_z+n_ordered_z)]) )
  
          kernel <- c(kernel, z_k)
        }
      }
      else{
        eff_score <- c(eff_score, 1)
        error_score <- c(error_score, NA)
        band_contin <- c(band_contin, NA)
        band_ordered <- c(band_ordered, NA)
        next
      }
      
      prob_kern <- prob_kern / sum(kernel)
      cat(length(prob_kern)," ", nrow(Z_sub), "\n")
      for (j in 1:nredo){
        
        Y_rand <- Y_sub[sample(nrow(Y_sub), sample_size, replace = TRUE, prob = prob_kern), ]
        eff_score_j <- c(eff_score_j, max( rowMins(t(t(Y_rand) / Y[i, ])) ))
      }
      
      
      eff_score <- c(eff_score, mean(eff_score_j)) 
      error_score <- c(error_score, sd(eff_score_j) / sqrt(nredo))
      band_contin <- c(band_contin, h[1:n_contin_z])
      band_ordered <- c(band_ordered, h[(n_contin_z+1):(n_contin_z+n_ordered_z)])
      cat(i ," ", eff_score[i]," ", h[1:n_contin_z], " ", h[(n_contin_z+1):(n_contin_z+n_ordered_z)] ,"\n")
      
    }
  }
  else stop("Error! The number of rows of matrixes must be equal")
  
  return(data.frame(eff = eff_score, error = error_score, bandwidth_contin = band_contin, 
                    bandwidth_ordered = band_ordered))
}


bandwidth_CI_v2 <- function(x, indic_col, ngood, nbad, Q=NULL, Q_ord=NULL)
{
  x_num = x[, indic_col]
  n_indic <- dim(as.matrix(x_num))[2]
  n_unit <- dim(as.matrix(x_num))[1]
  for (i in seq(1, n_indic)) {
    if (!is.numeric(x_num[, i])) {
      stop(paste("Data set not numeric at column:", i))
    }
  }
  for (i in seq(1, n_unit)) {
    for (j in seq(1, n_indic)) {
      if (is.na(x_num[i, j])) {
        message(paste("Pay attention: NA values at column:", 
                      i, ", row", j, ". Composite indicator has been computed, but results may be misleading, Please refer to OECD handbook, pg. 26."))
      }
    }
  }
  
  
  x_num = as.matrix(x_num)
  
  g <- ngood
  b <- nbad
  
  if(!is.null(Q) & is.null(Q_ord)){ 
    Q = as.matrix(Q)
    number_exogenous <- dim(Q)[2]
    nq_cont <- dim(Q)[2]
    nq_ord  <- 0
  }
  
  if(is.null(Q) & !is.null(Q_ord)){ 
    Q_ord = as.matrix(Q_ord)
    number_exogenous <- dim(Q_ord)[2]
    nq_cont <- 0
    nq_ord  <- dim(Q_ord)[2]
  }
  
  if(!is.null(Q) & !is.null(Q_ord)){ 
    Q = as.matrix(Q)
    Q_ord = as.matrix(Q_ord)
    number_exogenous <- dim(Q)[2] + dim(Q_ord)[2] 
    nq_cont <- dim(Q)[2]
    nq_ord  <- dim(Q_ord)[2]
  }
  
  # Create (n x (b+g)) matrix with only FALSE-values
  flag_ii<-matrix(FALSE,nrow=length(x_num[,1]),ncol=(b+g))
  # Flag all data values which satisfy the conditions below with TRUE-values
  
  flag_gi <- apply(matrix(x_num[,1:g],ncol=g), 2, function(x) x<=median(x))
  flag_bi <- apply(matrix(x_num[,(g+1):(g+b)],ncol=b), 2, function(x) x>=median(x))
  flag_ii <- cbind(flag_gi,flag_bi)
  
  # Flag the NUTS II regions which satisfy the (b+g) conditions
  flagg <- matrix(nrow=length(x_num[,1]),ncol=1)
  for (j in (1:length(x_num[,1])))
  {
    flagg[j] <- all(flag_ii[j,])
  }
  # Flag for the observations who satisfy all criteria the data and exogenous variables with TRUE-values 
  if(!is.null(Q) & is.null(Q_ord)){ 
    flag_QQ <- subset(Q,subset=flagg,drop=TRUE)   # the continuous exogenous variables
    data_fr1 <- data.frame(flag_QQ)
  }
  
  if(is.null(Q) & !is.null(Q_ord)){ 
    flag_ordd <- subset(Q_ord,subset=flagg,drop=TRUE)  # the discrete (ordered) exogenous variables
    data_fr1 <- data.frame(flag_ordd)
  }
  
  if(!is.null(Q) & !is.null(Q_ord)){ 
    flag_QQ <- subset(Q,subset=flagg,drop=TRUE)   # the continuous exogenous variables
    flag_ordd <- subset(Q_ord,subset=flagg,drop=TRUE)  # the discrete (ordered) exogenous variables
    # Collect the concerned data in a dataframe
    data_fr1 <- data.frame(flag_QQ,flag_ordd)
  }
  
  flag_yy <- subset(x_num,subset=flagg,drop=TRUE)   # the performance variables
  
  # Compute the bandwiths using the performance data and exogenous data
  bw1 <-  npcdensbw(ydat=flag_yy,xdat=data_fr1,cykertype="epanechnikov",cxkertype="epanechnikov",bwmethod="cv.ls",oxkertype="liracine",itmax = 1000)   ##set itmax to 1000 nmulti to 4 and tol high
  # Take the bandwiths for the exogenous variables
  bwfull <- bw1$xbw
  
  # Create an empty (nx1)-matrix
  aha<- matrix(nrow=length(x_num[,1]),ncol=1)
  # Create an empty matrix to store the endogenous bandwiths for the exogenous variables
  bw_cx <- matrix(nrow=length(x_num[,1]),ncol=number_exogenous)
  # Set i equal to the start value 1 (start loop with DMU 1)
  i <- 1
  # Define the loop to indicate which NUTS II regions dominate region i
  for (i in (1:length(x_num[,1])))
  {
    #Define a (n x (b+g)) matrix with all FALSE values
    flag_i<-matrix(FALSE,nrow=length(x_num[,1]),ncol=(b+g))
    # Flag the data values that dominate the values of DMU i
    flag_g <- apply(matrix(x_num[,1:g],ncol=g), 2, function(x) x<=x[i])
    flag_b <- apply(matrix(x_num[,(g+1):(g+b)],ncol=b), 2, function(x) x>=x[i])
    flag_i <- cbind(flag_g,flag_b)
    #Define a (n x 1) matrix 
    flag <- matrix(nrow=length(x_num[,1]),ncol=1)
    # Define a loop to indicate which NUTS II-regions dominate region i
    for (j in (1:length(x_num[,1])))
    {
      flag[j] <- all(flag_i[j,])
    }
    # print(i)
    # Take the subset of data values that correspond to the flagged NUTS II-regions 
    
    # Flag for the observations who satisfy all criteria the data and exogenous variables with TRUE-values 
    if(!is.null(Q) & is.null(Q_ord)){ 
      flag_Q <- matrix(subset(Q,subset=flag,drop=TRUE), ncol=nq_cont)   # the continuous exogenous variables
      # Collect the concerned data (continuous and discrete exogenous variables and performance data of dominating regions) in a dataframe 
      data_fr <- data.frame(flag_Q)
    }
    
    if(is.null(Q) & !is.null(Q_ord)){ 
      flag_ord <- matrix(subset(Q_ord,subset=flag,drop=TRUE), ncol=nq_ord)  # the discrete (ordered) exogenous variables
      # Collect the concerned data (continuous and discrete exogenous variables and performance data of dominating regions) in a dataframe 
      data_fr <- data.frame(flag_ord)
    }
    
    if(!is.null(Q) & !is.null(Q_ord)){ 
      flag_Q <- matrix(subset(Q,subset=flag,drop=TRUE), ncol=nq_cont)   # the continuous exogenous variables
      flag_ord <-  matrix(subset(Q_ord,subset=flag,drop=TRUE), ncol=nq_ord)  # the discrete (ordered) exogenous variables
      # Collect the concerned data (continuous and discrete exogenous variables and performance data of dominating regions) in a dataframe 
      data_fr <- data.frame(flag_Q, flag_ord)
    }
    
    flag_y <- subset(x_num,subset=flag,drop=TRUE)      # the performance variables
    
    #print(sum(flag))     # display the number of NUTS-regions dominating region i
    aha[i]=sum(flag)
    # Define a loop to determine the proper bandwith based on the number of dominating regions
    if((sum(flag))<=40)               # if number of dominating regions lower than or equal to threshold
    {
      bw_cx[i,] <-bwfull              # bandwith equal to what was derived endogenously above
    } else
    {                               # if number of dominating regions higher than threshold
      
      # Compute the bandwiths using the performance data and exogenous data
      bw <-  npcdensbw(ydat=flag_y,xdat=data_fr,cykertype="epanechnikov",cxkertype="epanechnikov",bwmethod="cv.ls",oxkertype="liracine",itmax = 10000,nmulti=4)   ##set itmax to 1000 nmulti to 4 and tol high
      # Take the bandwiths for the exogenous variables
      bw_cx[i,] <- bw$xbw
      #      print(bw_cx[i,])
    }
  }
  
  # save the output of the first step
  #write.matrix(bw_cx, file = "welbandwidth.txt", sep ="\t")
  
  # before opening the .txt file, first reset the options for system separators
  # Open the file with the endogenously estimated bandwiths (for the exogenous variables)
  #bw_cx<-read.table("welbandwidth.txt", sep = "\t")
  bw_cx2=bw_cx
  # Adjust the estimated bandwiths for the discrete (ordered) exogenous variables
  if(!is.null(Q) & is.null(Q_ord)){ 
    bw_cx2 <- bw_cx
  }
  
  if(is.null(Q) & !is.null(Q_ord)){ 
    bw_cx2[apply(bw_cx, 2, function(x) x>0.9999999)]<-0.9999
  }
  
  if(!is.null(Q) & !is.null(Q_ord)){ 
    bw_cx2[,nq_cont+1:(nq_cont+nq_ord-2)][apply(bw_cx[,nq_cont+1:(nq_cont+nq_ord-2)], 2, function(x) x>0.9999999)]<-0.9999
  }
  
  # Specify the final version of the matrix with endogenous bandwiths
  bw_cx<-bw_cx2
  
  #write.matrix(bw_cx, file = "welbandwidth_adj.txt", sep ="\t")
  
  r <- list(bandwidth = bw_cx,ci_method = "bandwidth_CI")
  r$call <- match.call()
  class(r) <- "CI"
  r
  
}


ci_rbod_constr_bad_Q_v2 <- function (x, indic_col, ngood=1, nbad=1, low_w=0, pref=NULL, M, B, Q=NULL, Q_ord=NULL, bandwidth) 
{
  x_num = x[, indic_col]
  n_indic <- dim(as.matrix(x_num))[2]
  n_unit <- dim(as.matrix(x_num))[1]
  for (i in seq(1, n_indic)) {
    if (!is.numeric(x_num[, i])) {
      stop(paste("Data set not numeric at column:", i))
    }
  }
  for (i in seq(1, n_unit)) {
    for (j in seq(1, n_indic)) {
      if (is.na(x_num[i, j])) {
        message(paste("Pay attention: NA values at column:", 
                      i, ", row", j, ". Composite indicator has been computed, but results may be misleading, Please refer to OECD handbook, pg. 26."))
      }
    }
  }
  if(!is.null(pref)){
    if (!is.character(pref)) {
      stop(paste("The preference vector (pref) is not of character type"))
    }
  }
  if(n_indic!=(ngood+nbad)){
    stop(paste("The number of simple indicators is not equal to the sum of the
               number of good and bad outputs, please insert the right number
               of simple indicators in x or choose the right number of columns
               in indic_col or indicate the right number of good and bad outputs
               in ngood and nbad."))
  }
  if(ngood==0){
    stop(paste("The number of good outputs (ngood) has to be greater than 0."))
  }
  if(nbad==0){
    stop(paste("The number of bad outputs (nbad) has to be greater than 0."))
  }
  
  if(!is.null(Q) & is.null(Q_ord)){ 
    Q = as.matrix(Q)
    #   number_exogenous <- dim(Q)[2]
    #   nq_cont <- dim(Q)[2]
    #   nq_ord  <- 0
  }
  
  if(is.null(Q) & !is.null(Q_ord)){ 
    Q_ord = as.matrix(Q_ord)
    #     number_exogenous <- dim(Q_ord)[2]
    #     nq_cont <- 0
    #     nq_ord  <- dim(Q_ord)[2]
  }
  
  if(!is.null(Q) & !is.null(Q_ord)){ 
    Q = as.matrix(Q)
    Q_ord = as.matrix(Q_ord)
    #     number_exogenous <- dim(Q)[2] + dim(Q_ord)[2] 
    #     nq_cont <- dim(Q)[2]
    #     nq_ord  <- dim(Q_ord)[2]
  }
  
  
  x_num = as.matrix(x_num)
  
  g <- ngood
  b <- nbad
  p <- length(pref)
  
  bw_cx <- bandwidth
  
  Avg <- apply(x_num, 2, mean)
  
  
  score <- rep(0, times = n_unit)
  eff <- rep(0, times = n_unit)
  w <- cbind(matrix(0, nrow = n_unit, ncol = (n_indic)))
  final <- cbind(matrix(0, nrow = n_unit, ncol = (1+n_indic)))
  kerzi <- matrix(nrow = n_unit, ncol=1)
  
  f.sign <- c(rep(-1, g), rep(1, b))
  
  for (i in 1:n_unit) {
    
    if(!is.null(pref))
    {
      # display the number of the NUTS II regions under assessment
      #print(i)
      # Collect the data in a dataframe
      if(!is.null(Q) & is.null(Q_ord)){ 
        dat <- data.frame(Q)
      }
      
      if(is.null(Q) & !is.null(Q_ord)){ 
        dat <- data.frame(Q_ord)
      }
      
      if(!is.null(Q) & !is.null(Q_ord)){ 
        dat <- data.frame(Q,Q_ord)
      }
      
      tdata <- dat[i,]
      # Estimate the densities
      kerz <- npksum(bws=t(bw_cx[i,]),txdat=dat, exdat=tdata, return.kernel.weights=TRUE,cykertype="epanechnikov",cxkertype="epanechnikov",oxkertype="liracine") 
      
      #kerz <- npudens(bws=t(bw_cx[i,]),cykertype="epanechnikov",cxkertype="epanechnikov",oxkertype="liracine",tdat=tdata,edat=dat)
      kerz<-kerz$kw
      kerzi[i,]<-sum(kerz)
      # collect the performance data
      y_aggr <- data.frame(matrix(x_num[,(g+1):(g+b)],ncol=b),matrix(x_num[,1:g],ncol=g),kerz=kerz)
      yk <- y_aggr[i,c(1:(b+g))]
      q_y <- (b+g)
      subsety<-y_aggr
      score_B <- matrix(nrow=B,ncol=1)
      w_B <- cbind(matrix(0, nrow = B, ncol = (n_indic)))
      final_B <- cbind(matrix(0, nrow = B, ncol = (1+n_indic)))
      fl <- seq(1:n_unit)
      k <- 1
      
      for(k in 1:B)
      {
        draw <- sample(fl, size=M, replace=TRUE, prob=subsety$kerz)
        x_num_sample <- x_num[draw,]
        
        # Different cases of ordVWR
        # Case 1: the preference vector has only 2 elements => only ordVWR_type1  
        if(length(pref)==2)
        {
          f.obj <- c(f.sign * x_num[i, ])
          
          f.dir <- c("==",rep(">=", times = M), rep(">=", times = g+b), 
                     rep(">=",times=(n_indic-1)),
                     rep(">=", times = g+b))
          
          f.rhs <- c(1, rep(0, times = M), rep(0, times = g+b), 
                     rep(0,times=(n_indic-1)),
                     rep(0, times = g+b))  
          
          Constr1 <- c(x_num[i,])
          Constr2 <- cbind(t(f.sign * t(x_num_sample))) 
          
          I = diag(x = 1, n_indic)
          VWRg = t((I - low_w)*as.vector(Avg)) 
          
          ### OrdVWR_type1
          I = diag(x = 1, n_indic)
          colnames(I) <- colnames(x_num)
          I[,pref[1]] <- 0
          
          M_imp =  matrix(0,nrow=(n_indic), ncol=n_indic)
          colnames(M_imp) <- colnames(x_num)
          M_imp[,pref[1]] <- 1
          
          ordVWR_type1 <- (-I*Avg) + (M_imp*Avg[pref[1]])
          ordVWR_type1 <- ordVWR_type1[-which(colnames(M_imp)==pref[1]),]
          
          ### OrdVWR_type2
          I2 = diag(x = 1, n_indic)
          colnames(I2) <- colnames(x_num)
          I2[,pref[length(pref)]] <- 0
          
          M_imp2 =  matrix(0,nrow=(n_indic), ncol=n_indic)
          colnames(M_imp2) <- colnames(x_num)
          M_imp2[,which(colnames(M_imp)==pref[length(pref)])] <- -Avg[pref[length(pref)]]
          ordVWR_type2 <- (I2*Avg)+ M_imp2
          
          '%!in%' <- function(x,y)!('%in%'(x,y))
          ordVWR_type2 <- ordVWR_type2[-which(colnames(M_imp) %!in% pref),]
          ordVWR_type2 <- ordVWR_type2[-which(colnames(M_imp)==pref[1]),]
          ordVWR_type2 <- ordVWR_type2[-which(colnames(M_imp)==pref[length(pref)]),]
          
          
          Pweight_g <- cbind(diag(1, nrow = g, ncol = g), matrix(0,nrow=g,ncol=b))
          Pweight_b <- cbind(matrix(0,nrow=b,ncol=g),diag(1, nrow = b, ncol = b))
          
          f.con <- rbind(Constr1, Constr2, Pweight_g, Pweight_b, VWRg, ordVWR_type1)
          
          jj <- lp("min", f.obj, f.con, f.dir, f.rhs)
          score_B[k] <- jj$objval
          w_B[k, ] <- rbind(jj$solution)
          final_B[k, ] <- c(score_B[k], w_B[k, ])
        }
        
        # Case 2: the preference vector has all indicators         
        
        if(length(pref)==n_indic){    
          
          f.obj <- c(f.sign * x_num[i, ])
          
          f.dir <- c("==",rep(">=", times = M), rep(">=", times = g+b), 
                     rep(">=",times=(n_indic-1)), rep(">=",times=(p-2)),
                     rep(">=", times = g+b))
          
          f.rhs <- c(1, rep(0, times = M), rep(0, times = g+b), 
                     rep(0,times=(n_indic-1)), rep(0,times=(p-2)),
                     rep(0, times = g+b))  
          
          Constr1 <- c(x_num[i,])
          Constr2 <- cbind(t(f.sign * t(x_num_sample))) 
          
          I = diag(x = 1, n_indic)
          VWRg = t((I - low_w)*as.vector(Avg)) 
          
          ### OrdVWR_type1
          I = diag(x = 1, n_indic)
          colnames(I) <- colnames(x_num)
          I[,pref[1]] <- 0
          
          M_imp =  matrix(0,nrow=(n_indic), ncol=n_indic)
          colnames(M_imp) <- colnames(x_num)
          M_imp[,pref[1]] <- 1
          
          ordVWR_type1 <- (-I*Avg) + (M_imp*Avg[pref[1]])
          ordVWR_type1 <- ordVWR_type1[-which(colnames(M_imp)==pref[1]),]
          
          ### OrdVWR_type2
          I2 = diag(x = 1, n_indic)
          colnames(I2) <- colnames(x_num)
          I2[,pref[length(pref)]] <- 0
          
          M_imp2 =  matrix(0,nrow=(n_indic), ncol=n_indic)
          colnames(M_imp2) <- colnames(x_num)
          rownames(M_imp2) <- colnames(x_num)
          M_imp2[,which(colnames(M_imp2)==pref[length(pref)])] <- -Avg[pref[length(pref)]]
          ordVWR_type2 <- (I2*Avg)+ M_imp2
          colnames(ordVWR_type2) <- colnames(x_num)
          rownames(ordVWR_type2) <- colnames(x_num)
          
          #'%!in%' <- function(x,y)!('%in%'(x,y))
          #ordVWR_type2 <- ordVWR_type2[-which(colnames(ordVWR_type2) %!in% pref),]
          ordVWR_type2 <- ordVWR_type2[-which(rownames(ordVWR_type2)==pref[length(pref)]),]
          ordVWR_type2 <- ordVWR_type2[-which(rownames(ordVWR_type2)==pref[1]),]
          
          Pweight_g <- cbind(diag(1, nrow = g, ncol = g), matrix(0,nrow=g,ncol=b))
          Pweight_b <- cbind(matrix(0,nrow=b,ncol=g),diag(1, nrow = b, ncol = b))
          
          f.con <- rbind(Constr1, Constr2, Pweight_g, Pweight_b, VWRg, ordVWR_type1, ordVWR_type2)
          
          jj <- lp("min", f.obj, f.con, f.dir, f.rhs)
          score_B[k] <- jj$objval
          w_B[k, ] <- rbind(jj$solution)
          final_B[k, ] <- c(score_B[k], w_B[k, ])
        }
        
        # Case 3: the preference vector has a number of elements >2 and < of the number of indicators 
        
        if(length(pref)>2 & length(pref)<n_indic){
          
          f.obj <- c(f.sign * x_num[i, ])
          
          f.dir <- c("==",rep(">=", times = M), rep(">=", times = g+b), 
                     rep(">=",times=(n_indic-1)), rep(">=",times=(p-2)),
                     rep(">=", times = g+b))
          
          f.rhs <- c(1, rep(0, times = M), rep(0, times = g+b), 
                     rep(0,times=(n_indic-1)), rep(0,times=(p-2)),
                     rep(0, times = g+b))  
          
          Constr1 <- c(x_num[i,])
          Constr2 <- cbind(t(f.sign * t(x_num_sample))) 
          
          I = diag(x = 1, n_indic)
          VWRg = t((I - low_w)*as.vector(Avg)) 
          
          ### OrdVWR_type1
          I = diag(x = 1, n_indic)
          colnames(I) <- colnames(x_num)
          I[,pref[1]] <- 0
          
          M_imp =  matrix(0,nrow=(n_indic), ncol=n_indic)
          colnames(M_imp) <- colnames(x_num)
          M_imp[,pref[1]] <- 1
          
          ordVWR_type1 <- (-I*Avg) + (M_imp*Avg[pref[1]])
          ordVWR_type1 <- ordVWR_type1[-which(colnames(M_imp)==pref[1]),]
          
          ### OrdVWR_type2
          I2 = diag(x = 1, n_indic)
          colnames(I2) <- colnames(x_num)
          I2[,pref[length(pref)]] <- 0
          
          M_imp2 =  matrix(0,nrow=(n_indic), ncol=n_indic)
          colnames(M_imp2) <- colnames(x_num)
          rownames(M_imp2) <- colnames(x_num)
          M_imp2[,which(colnames(M_imp2)==pref[length(pref)])] <- -Avg[pref[length(pref)]]
          ordVWR_type2 <- (I2*Avg)+ M_imp2
          colnames(ordVWR_type2) <- colnames(x_num)
          rownames(ordVWR_type2) <- colnames(x_num)
          
          '%!in%' <- function(x,y)!('%in%'(x,y))
          ordVWR_type2 <- ordVWR_type2[-which(rownames(ordVWR_type2) %!in% pref),]
          ordVWR_type2 <- ordVWR_type2[-which(rownames(ordVWR_type2)==pref[length(pref)]),]
          ordVWR_type2 <- ordVWR_type2[-which(rownames(ordVWR_type2)==pref[1]),]
          
          Pweight_g <- cbind(diag(1, nrow = g, ncol = g), matrix(0,nrow=g,ncol=b))
          Pweight_b <- cbind(matrix(0,nrow=b,ncol=g),diag(1, nrow = b, ncol = b))
          
          f.con <- rbind(Constr1, Constr2, Pweight_g, Pweight_b, VWRg, ordVWR_type1, ordVWR_type2)
          
          jj <- lp("min", f.obj, f.con, f.dir, f.rhs)
          score_B[k] <- jj$objval
          w_B[k, ] <- rbind(jj$solution)
          final_B[k, ] <- c(score_B[k], w_B[k, ])
        }
      }
    }else{
      # display the number of the NUTS II regions under assessment
      #print(i)
      # Collect the data in a dataframe
      if(!is.null(Q) & is.null(Q_ord)){ 
        dat <- data.frame(Q)
      }
      
      if(is.null(Q) & !is.null(Q_ord)){ 
        dat <- data.frame(Q_ord)
      }
      
      if(!is.null(Q) & !is.null(Q_ord)){ 
        dat <- data.frame(Q, Q_ord)
      }
      
      tdata <- dat[i,]
      # Estimate the densities
      kerz <- npksum(bws=t(bw_cx[i,]),txdat=dat, exdat=tdata, return.kernel.weights=TRUE,cykertype="epanechnikov",cxkertype="epanechnikov",oxkertype="liracine") 
      
      #kerz <- npudens(bws=t(bw_cx[i,]),cykertype="epanechnikov",cxkertype="epanechnikov",oxkertype="liracine",tdat=tdata,edat=dat)
      kerz<-kerz$kw
      kerzi[i,]<-sum(kerz)
      # collect the performance data
      y_aggr <- data.frame(matrix(x_num[,(g+1):(g+b)],ncol=b),matrix(x_num[,1:g],ncol=g),kerz=kerz)
      yk <- y_aggr[i,c(1:(b+g))]
      q_y <- (b+g)
      subsety<-y_aggr
      score_B <- matrix(nrow=B,ncol=1)
      w_B <- cbind(matrix(0, nrow = B, ncol = (n_indic)))
      final_B <- cbind(matrix(0, nrow = B, ncol = (1+n_indic)))
      fl <- seq(1:n_unit)
      k <- 1
      
      for(k in 1:B)
      {
        draw <- sample(fl, size=M, replace=TRUE, prob=subsety$kerz)
        x_num_sample <- x_num[draw,]
        
        f.obj <- c(f.sign * x_num[i, ])
        
        f.dir <- c("==",rep(">=", times = M), rep(">=", times = g+b), 
                   rep(">=", times = g+b))
        f.rhs <- c(1, rep(0, times = M), rep(0, times = g+b), 
                   rep(0, times = g+b))  
        
        Constr1 <- c(x_num[i,])
        Constr2 <- cbind(t(f.sign * t(x_num_sample))) 
        
        I = diag(x = 1, n_indic)
        VWRg = t((I - low_w)*as.vector(Avg)) 
        Pweight_g <- cbind(diag(1, nrow = g, ncol = g), matrix(0,nrow=g,ncol=b))
        Pweight_b <- cbind(matrix(0,nrow=b,ncol=g),diag(1, nrow = b, ncol = b))
        
        f.con <- rbind(Constr1, Constr2, VWRg, Pweight_g, Pweight_b)
        jj <- lp("min", f.obj, f.con, f.dir, f.rhs)
        score_B[k] <- jj$objval
        w_B[k, ] <- rbind(jj$solution)
        final_B[k, ] <- c(score_B[k], w_B[k, ])
      }
    }
    score[i]<- apply(score_B,2,mean) 
    w[i,]<- apply(w_B,2,mean)  
    final[i,]<- apply(final_B,2,mean)   
  }
  eff <- (1/(1+score))
  
  colnames(w) <- colnames(x_num)
  
  # Write the composite scores and the indicator target values in one result matrix
  
  target <- matrix(nrow=n_unit,ncol=(b+g)) 
  for(j in 1:g)
  {
    target[,j] <- x_num[,j] + (Avg[j]*score)
  }
  for(j in (g+1):(g+b))
  {
    target[,j] <- x_num[,j] - (Avg[j]*score)
  }
  colnames(target) <- colnames(x_num)
  
  r <- list(ci_rbod_constr_bad_Q_est = eff, ci_rbod_constr_bad_Q_weights = w, 
            ci_rbod_constr_bad_Q_target = target, ci_method = "rbod_constr_bad_Q")
  r$call <- match.call()
  class(r) <- "CI"
  r
}

