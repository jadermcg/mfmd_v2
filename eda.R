eda = function(Y, theta, n_pop, stop_loop) {
  
  N = n1 = n2 = len(Y)
  
  
  for(i in 1:stop_loop) {
    
    #before = information_content(pssm_matrixC(theta$theta_1, theta$theta_2))
    
    # Inicial population
    P_theta_1 = foreach(i = 1:n_pop) %do% rdirichlet(1, theta$theta_1*n1)
    P_theta_2 = foreach(i = 1:n_pop) %do% rdirichlet(w, theta$theta_2*n2)
    P_priori = foreach(i = 1:n_pop) %do% rdirichlet(1, theta$priori*N)
    
    
    # Expectation maximization
    params = em(Y, P_theta_1, P_theta_2, P_priori, stop_loop)
    P_theta_1 = params$theta_1
    P_theta_2 = params$theta_2
    P_priori = params$priori
    
    idx = which.max(mapply(function(a,b) information_content(pssm_matrixC(a, b)), a = P_theta_1, b = P_theta_2, SIMPLIFY = F))
    theta$theta_1 = P_theta_1[[idx]]
    theta$theta_2 = P_theta_2[[idx]]
    theta$priori = P_priori[[idx]]
    n1 = params$n1[[idx]]
    n2 = params$n2[[idx]]
    
    #after = information_content(pssm_matrixC(theta$theta_1, theta$theta_2))
    
    #fprintf('Before %f\nAfter %f\n\n', before, after)
    
  }
  
  return (theta)
  
}
