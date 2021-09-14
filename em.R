em = function(Y, theta_1, theta_2, priori, stop_loop) {
  log_likes = list()
  inf_cont = list()
  N = len(Y)
  L = len(Y[[1]])
  n = N*L
  # Main loop
 
  pb <- progress_bar$new(total = stop_loop)
  for(i in 1:stop_loop) {
    pb$tick()
    
    # Expection
    E = expectation(Y, theta_1, theta_2, priori)
    rn1 = E$rn1
    rn2 = E$rn2
    lhood = E$lhood
    n1 = N1(rn1)
    n2 = N2(rn2)
    
    # Priori maximization
    priori = maximization_priori(n, n1, n2)
    
    # Theta_2 maximization
    theta_2 = maximization_theta_2(Y, rn2, n2)
    
    # Theta_1 maximization
    theta_1 = maximization_theta_1(Y, rn1, n1)
    
    # Log-likelihood evolution
    #log_likes[[i]] = lapply(lhood, function(x) -sum(log(x)))
    
    # Pssm scores
    #inf_cont[[i]] = mapply(function(t1, t2) information_content(pssm_matrixC(t1, t2)), t1=theta_1, t2=theta_2, SIMPLIFY = F)
    
  }
  
  return ( list(theta_1 = theta_1, theta_2 = theta_2, priori = priori, 
                rn1 = rn1, rn2 = rn2, log_likes = log_likes, inf_cont = inf_cont) )
  
}


expectation = function(Y, theta_1, theta_2, priori) {
  lhood_1 = mapply(function(a, b) likelihood_t1C_v3(Y, w, a, b[1]), a = theta_1, b = priori, SIMPLIFY = F)
  lhood_2 = mapply(function(a, b) likelihood_t2C_v3(Y, w, a, b[2]), a = theta_2, b = priori, SIMPLIFY = F)
  
  lhood = mapply(function(l1, l2) l1 + l2, l1 = lhood_1, l2 = lhood_2, SIMPLIFY = F)
  rn1 = mapply(function(l1, l) l1/l, l1 = lhood_1, l = lhood, SIMPLIFY = F)
  rn2 = mapply(function(l2, l) l2/l, l2 = lhood_2, l = lhood, SIMPLIFY = F)
  
  return(list(rn1 = rn1, rn2 = rn2, lhood = lhood))
}

maximization_priori = function(n, n1, n2) {
  
  P1 = lapply(n1, function(x) x/n )
  P2 = lapply(n2, function(x) x/n )
  return ( mapply(function(a,b) cbind(a,b), a = P1, b = P2, SIMPLIFY = F)  )
}

maximization_theta_1 = function(Y, rn1, n1) {
  M = mapply(function(x) maximization_v2(Y, x, w), x = rn1, SIMPLIFY = F)
  return ( mapply(function(x,y) colMeans(x/y) ,x = M, y = n1, SIMPLIFY = F) )
}

maximization_theta_2 = function(Y, rn2, n2) {
  M = mapply(function(x) maximization_v2(Y, x, w), x = rn2, SIMPLIFY = F)
  return ( mapply(function(x,y) x/y ,x = M, y = n2, SIMPLIFY = F) )
}

N1 = function(rn1) {
  return ( lapply(rn1, sum)  )
}

N2 = function(rn2) {
  return ( lapply(rn2, sum) )
}
