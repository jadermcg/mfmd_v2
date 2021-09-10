em = function(Y, theta_1, theta_2, priori) {
  log_likes = list()
  inf_cont = list()
  N = len(Y)
  L = len(Y$seq1)
  
  # Main loop
  # stop_loop = round( log(N) * log(L) )
  stop_loop = 50
  pb <- progress_bar$new(total = stop_loop)
  for(i in 1:stop_loop) {
    pb$tick()
    
    # Expection
    Z = expectation(Y, theta_1, theta_2, priori)
    rn1 = Z$rn1
    rn2 = Z$rn2
    lhood_1 = Z$l1
    lhood_2 = Z$l2
    n1 = N1(rn1)
    n2 = N2(rn2)
    
    # Priori maximization
    priori = maximization_priori(Z, n, n1, n2)
    
    # Theta_2 maximization
    theta_2 = maximization_theta_2(Y, rn2, n2)
    
    # Theta_1 maximization
    theta_1 = maximization_theta_1(theta_2)
    
    # Log-likelihood evolution
    log_likes[[i]] = mapply(function(l1, l2) -sum(log(l1+l2)), l1=lhood_1, l2=lhood_2, SIMPLIFY = F)
    
    # Pssm scores
    inf_cont[[i]] = mapply(function(t1, t2) information_content(pssm_matrixC(t1, t2)), t1=theta_1, t2=theta_2, SIMPLIFY = F)
    
  }
  
  return ( list(theta_1 = theta_1, theta_2 = theta_2, priori = priori, 
                rn1 = rn1, rn2 = rn2, log_likes = log_likes, inf_cont = inf_cont) )
  
}


expectation = function(Y, theta_1, theta_2, priori) {
  lhood_1 = mapply(function(a,b) likelihood_t1C_v2(Y, w, a, b), a = theta_1, b = priori, SIMPLIFY = F)
  lhood_2 = mapply(function(a,b) likelihood_t2C_v2(Y, w, a, b), a = theta_2, b = priori, SIMPLIFY = F)
  
  lhood = mapply(function(l1, l2) l1 + l2, l1 = lhood_1, l2 = lhood_2, SIMPLIFY = F)
  rn1 = mapply(function(l1, l) l1/l, l1 = lhood_1, l = lhood, SIMPLIFY = F)
  rn2 = mapply(function(l2, l) l2/l, l2 = lhood_2, l = lhood, SIMPLIFY = F)
  
  return(list(rn1 = rn1, rn2 = rn2, l1 = lhood_1, l2 = lhood_2))
}

maximization_priori = function(Z, n, n1, n2) {
  rn1 = Z[[1]]
  rn2 = Z[[2]]
  
  P1 = lapply(n1, function(x) x/n )
  P2 = lapply(n2, function(x) x/n )
  return ( mapply(function(a,b) cbind(a,b), a = P1, b = P2, SIMPLIFY = F)  )
}

maximization_theta_1 = function(theta_2) {
  return ( lapply(theta_2, colMeans)  )
}

maximization_theta_2 = function(Y, rn2, n2) {
  M = mapply(function(x) maximization_t2C_v2(Y, x, w), x = rn2, SIMPLIFY = F)
  return ( mapply(function(x,y) x/y ,x = M, y = n2, SIMPLIFY = F) )
}

N1 = function(rn1) {
  return ( lapply(rn1, sum)  )
}

N2 = function(rn2) {
  return ( lapply(rn2, sum) )
}

"
# Results
cbind(which( rn_2  > 0.7 ) %/% L + 1, which( rn_2 > 0.7 ) %% L)

# Plots
par(mfrow = c(6,3), mar = c(1,1,1,1) )
f = foreach (i = 1:N) %do% plot( rn_2[ (((i-1)*L)+1) : (i*L)  ] , type='l')
pssm = pssm_matrix(theta_1, theta_2)
scores = foreach(i = 1:n, .combine = 'rbind') %dopar% { foreach(j = 1:w, .combine = 'c') %dopar%{ pssm[j, Y[i,j]] } }
dimnames(scores) = NULL
scores = rowSums(scores)
scores = Reshape(scores, n)
q = quantile(scores, .975)
cbind(which( scores  > q ) %/% L + 1, which( scores > q ) %% L)

par(mfrow = c(1,1), mar = c(1,1,1,1) )
plot(sort(scores), type = 'l')
abline(v = n - sum(scores > q), col='red')
"
