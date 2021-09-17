em = function(Y, theta, stop_loop) {
  
  N = len(Y)
  L = len(Y[[1]])
  w = len(Y[[1]][[1]])
  n = N*L
  log_like = c()
  inf_cont = c()
  
  pb <- progress_bar$new(total = stop_loop)
  for (i in 1:stop_loop) {
    pb$tick()
    
    # Expectation
    E = expectation(Y = Y, theta = theta, w = w)
    
    # Maximizarion priori
    n1 = sum(E$rn1)
    n2 = sum(E$rn2)
    priori = maximization_priori(n = n, n1 = n1, n2 = n2)
    
    # Maximization theta_1 e theta_2
    theta_1 = maximization_theta_1(Y = Y, w = w, rn1 = E$rn1, n1 = n1)
    theta_2 = maximization_theta_2(Y = Y, w = w, rn2 = E$rn2, n2 = n2)
    theta = list(theta_1 = theta_1, theta_2 = theta_2, priori = priori[[1]])
    
    
    # Log_like and information content
    log_like[i] = sum(log(E$lhood))
    inf_cont[i] = information_content(pssm_matrixC(theta_1, theta_2))
    
  }
  
  return (list(theta = theta, log_like = log_like, inf_cont = inf_cont))
  
}

maximization_theta_1 = function(Y, w, rn1, n1) {
  M = maximization_v2(Y = Y, rn = rn1, w = w)
  return ( colMeans(M/n1) )
}

maximization_theta_2 = function(Y, w, rn2, n2) {
  M = maximization_v2(Y = Y, rn = rn2, w = w)
  return ( M/n2  )
}

maximization_priori = function(n, n1, n2) {
  
  p1 = n1/n
  p2 = n2/n
  return ( mapply(function(a,b) cbind(a,b), a = p1, b = p2, SIMPLIFY = F)  )
}

expectation = function(Y, theta, w) {
  
  priori = theta$priori
  theta_1 = theta$theta_1
  theta_2 = theta$theta_2
  p1 = priori[1]
  p2 = priori[2]
  
  lhood_1 = likelihood_t1C_v3(Y = Y, w = w, theta = theta_1, priori = p1)
  lhood_2 = likelihood_t2C_v3(Y = Y, w = w, theta = theta_2, priori = p2)
  
  lhood = lhood_1 + lhood_2
  rn1 =   lhood_1 / lhood
  rn2 =   lhood_2 / lhood
  
  return(list(rn1 = rn1, rn2 = rn2, lhood_1 = lhood_1, lhood_2 = lhood_2, lhood = lhood))
}