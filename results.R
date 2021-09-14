get_pos = function(rn1, rn2, L) {
  line = mapply(function(a, b) (which(a > b) %/% L) + 1 , a = rn2, b = rn1, SIMPLIFY = F)
  column = mapply(function(a, b) (which(a > b) %% L) , a = rn2, b = rn1, SIMPLIFY = F)
  return ( mapply(function(seq,b) cbind(seq,b), seq = line, b = column, SIMPLIFY = F) )
}

get_pssms = function(theta_1, theta_2) {
  return ( mapply(function(t1, t2) pssm_matrixC(t1, t2), 
                  t1 = theta_1, t2 = theta_2, SIMPLIFY = F))
}


get_pssm_scores = function(pssms) {
  return ( do.call(c, lapply(pssms, information_content) ))
}
