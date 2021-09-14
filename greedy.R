greedy_inicialization = function(Y, N, w, b) {
  # Y: sample data
  # w: width of motif
  # N: number of sequences
  
  # Inicial theta_1 values
  theta_1 = relative_nucleotide_freq(Y)
  # theta_1 = rep(.25, 4)
  
  # Do cartesian product into first and second sequences
  P = Y[[1]]
  Q = Y[[2]]
  TEMP = cartesian_prod(P, Q)
  
  # Get b first best motifs
  pfm = lapply(TEMP, function(temp) { lapply(temp, pfm_matrixC) })
  ppm = lapply(pfm, function(t) { lapply(t, ppm_matrixC) })
  pssm = lapply(ppm, function(x) { lapply(x,function(y){pssm_matrixC(theta_1, y)})})
  IC = lapply(pssm, function(x) lapply(x, information_content))
  
  POS = lapply(IC, which.max)
  R = foreach(i = 1:len(P)) %do% pfm[[i]][[POS[[i]]]]
  IC = foreach(i = 1:len(IC), .combine = 'c') %do% IC[[i]][[POS[[i]]]]
  P = R[order(IC, decreasing = T)[1:b]]
  
  # Main loop
  pb <- progress_bar$new(total = len(Y) - 2)
  for (k in 3:N )  {
    pb$tick()
    
    Q = Y[[k]]
    
    pfm = lapply(P, function(p) lapply(Q, function(q) update_pfmC(p,q) ) )
    ppm = lapply(pfm, function(t) { lapply(t, ppm_matrixC) })
    pssm = lapply(ppm, function(x) { lapply(x,function(y){pssm_matrixC(theta_1, y)})})
    IC = lapply(pssm, function(x) lapply(x, information_content))
    
    POS = lapply(IC, which.max)
    R = foreach(i = 1:len(pfm)) %do% pfm[[i]][[POS[[i]]]]
    IC = foreach(i = 1:len(IC), .combine = 'c') %do% IC[[i]][[POS[[i]]]]
    P = R[order(IC, decreasing = T)[1:b]]
  }
  
  # Return b best ppm's matrices found
  return ( lapply(P, ppm_matrixC) )
  
}


