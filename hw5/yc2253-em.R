

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
  r = rep(0, q);                        # rank for matrix Z's 
  Zlist = list();                       # Constructing the covariance matrices Z = zz' here
  for(i in 1:q){
    z = matrix(unlist(zlist[i]), nrow=n);
    r[i] = dim(z)[2];
    Zlist[i] = list(z %*% t(z));
  }
  
  
  # Initialize some variables
  d = 1.0;                              # difference in log-likelihood
  beta0 = beta;                         # The beta container
  u0 = u;                               # The u container
  P = X %*% ginv(t(X) %*% X) %*% t(X);  # Projection matrix
  counter = 0;                          # Counter of the iteration
  ll0 = -10                            # Calculate the loglikelihood
  
  
  
  # The EM Algorithm
  while((d > tol) & (counter < maxiter)){
    
    counter = counter + 1;
    sigma = matrix(0.0, n, n);          # Calculate sigma
    for(i in 1:q) sigma = sigma + u0[i] * matrix(unlist(Zlist[i]), ncol=n); 
    invsigma = solve(sigma);             # Inverse sigma

        
    # Random effects update
    b = invsigma %*% (Y - X %*% beta0); # This term seems to show up many times, avoid repetitive computation 
    for(i in 1:q){
      z = matrix(unlist(zlist[i]), ncol=n);
      Z = matrix(unlist(Zlist[i]), ncol=n);
      u0[i] = ( u0[i] * r[i] - sum(diag(u0[i]^2 * (z %*% invsigma %*% t(z)))) 
                        + u0[i]^2 * t(b) %*% (Z %*% b) )/r[i]; 
    }
    
    
    # Fixed effect update
    s = X %*% beta0 + u0[q] * b;
    beta0 = solve(t(X) %*% X, t(X) %*% P %*% s);
    
    
    # Log-likelihood
    ll1 = -0.5 * sum(r * log(abs(u0)) - u0);
    d = ll1 - ll0;
    ll0 = ll1; 
    
    # Print the fixed/random effect iteration and the statistics of this iteration
    print(sprintf(fmt = "-- Iteration # %5d ...", counter)); 
    print(u0);
    print(beta0);
    #sprintf(fmt = "   Random effect variable u =  ..", ); 
    #sprintf(fmt = "   Fixed effect variable beta =  ..", );
    print(sprintf(fmt = "   Improvement on LL d = %0.4f ..", d));
    
    
  }
  
  
  print("--------------- EM Done ---------------");
  print(sprintf(fmt="- Number of iterations = %5d.", counter));
  print(sprintf(fmt="- Final Log-Likelihood = %5f.", ll1)); 
  
  # return the random effect u0 and fixed effect beta0
  # Return the stuff we got from the em algorithm
  resultlist = list("random" = u0, "fixed" = beta0, "ll" = ll1); 
  return(resultlist); 
  
}

