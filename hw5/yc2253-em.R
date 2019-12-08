

library(MASS)

# The EM function for variance components: Y = X beta + Zu + E
# X: regression matrix
# Y: response
# Zlist: variance components, i.e. z vectors
# beta: Initial guess for the fixed effects
# u: initial guess for the random effects
# tol: convergence criterion for log-likelihood improvement
# maxiter: max number of iteration 

vcem = function(X, Y, zlist, beta, u, tol=1e-5, maxiter=10000){
  
  n = length(Y);                        # total dimension
  q = length(zlist);                    # number of variance components
  r = rep(n, q+1);                      # rank for matrix Z's
  Zlist = list();
  for(i in 1:q){
    r[i+1] = dim(zlist[i])[2];
    Zlist[i] = zlist[i] %*% t(zlist[i]);
  }
  
  
  
  
  
  
  
  
  # Initialize some variables
  d = 1.0;                              # difference in log-likelihood
  beta0 = beta;                         # The beta container
  u0 = u;                               # The u container
  P = X %*% ginv(t(X) %*% X) %*% t(X);  # Projection matrix
  z = cbind(rep(1.0, n), Z);            # The Z container with random error, E
  counter = 0;                          # Counter of the iteration
  ll0 = -Inf                            # Calculate the loglikelihood
  
  while((d > tol) & (counter < maxiter)){
    
    counter = counter + 1;
    sigma = diag(rep(u0[1], n));        # Calculate sigma
    for(i in 1:q){
      sigma = sigma + u0[i+1] * Zlist[i]; 
    }
    invsigma = ginv(sigma);             # Inverse sigma

        
    # Random effects update
    for(i in 1:length(u0)){
  
      b = invsigma %*% (Y - X %*% beta0);  # This term seems to show up many times, avoid repetitive computation 
      v = t(Z[, i]) %*% b;
      t = u0[i] - u0[i]^2 * (t(Z[, i]) %*% invsigma %*% Z[, i]) + u0[i]^2 * t(v) %*% v; 
      u0[i] = t / r[i];
      
    }
    
    
    # Fixed effect update
    s = X %*% beta0 + u0[1] * b;
    beta0 = solve(X, P %*% s);
    
    
    # Log-likelihood
    ll1 = -0.5 * sum(log(u0) - u0);
    d = ll1 - ll0;
    ll0 = ll1; 
    
    # Print the fixed/random effect iteration and the statistics of this iteration
    sprintf(fmt = "-- Iteration # %5d ...", counter); 
    sprintf(fmt = "   Random effect variable u =  ..", ); 
    sprintf(fmt = "   Fixed effect variable beta =  ..", );
    sprintf(fmt = "   Improvement on LL d = %0.4f ..", d);
    
    
  }
  
  
  # return the random effect u0 and fixed effect beta0
  # Return the stuff we got from the em algorithm
  resultlist = list("random" = u0, "fixed" = beta0, "ll" = finallikelihood); 
  return(resultlist); 
  
}

