

library(MASS)

# The EM function 


vcem = function(X, Y, Z, beta, u, tol=1e-5, maxiter=10000){
  
  # Check the dimensions 
  n = dim(X)[1]; 
  p = dim(X)[2]; 
  q = dim(Z)[2];

  if (length(Y) != n){
    print('The rows of X does not equal to length of Y .. '); 
    return(0); 
  }
  
  if (dim(Z)[1] != n){
    print('The rows of Z does not equal to length of Y .. '); 
    return(0); 
  }
  
  if (length(beta) != p){
    print('The length of beta does not equal to columns of X .. '); 
    return(0); 
  }
  
  if (length(u) != q+1){
    print('The length of u does not equal to columns of Z .. '); 
    return(0); 
  }
  
  
  # Some variables
  d = 1.0; # difference in log-likelihood
  beta0 = beta; # The beta container
  u0 = u; # The u container
  P = X %*% ginv(t(X) %*% X) %*% t(X); # Projection matrix
  z = cbind(rep(1.0, n), Z);
  counter = 0;
  
  
  while((d > tol) & (counter < maxiter)){
    
    counter = counter + 1;
    sigma = diag(u0) %*% Z %*% t(Z); # Calculate sigma
    invsigma = ginv(sigma); # Inverse sigma

        
    # Calculate random effects update
    for(i in 1:length(u0)){
  
      b = invsigma %*% (Y - X %*% beta0); 
      v = t(Z[, i]) %*% b;
      u0[i] = u0[i] - u0[i]^2 * (t(Z[, i]) %*% invsigma %*% Z[, i]) + u0[i]^2 * t(v) %*% v; 
    }
    
    
    # Calculate fixed effect update
    s = X %*% beta0 + u0[1] * b;
    beta0 = solve(X, P %*% s);
    
    
    # Calculate the log-likelihood
    
    
    
    # Print the fixed/random effect iteration and the statistics of this iteration
    
    
  }
  
  
  # return the random effect u0 and fixed effect beta0
  
  
}