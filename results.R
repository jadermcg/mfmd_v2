get_pos = function(rn1, rn2, N, L) {
  
  g = as.factor(foreach(i = 1:N, .combine = 'c') %do% rep(i, L))
  rn1 = split(rn1, g)
  rn2 = split(rn2, g)
  
  return ( mapply(function(a,b) which(a >= b), a = rn2, b = rn1, SIMPLIFY = F) )
}

get_pssms = function(theta_1, theta_2) {
  return ( mapply(function(t1, t2) pssm_matrixC(t1, t2), 
                  t1 = theta_1, t2 = theta_2, SIMPLIFY = F))
}


get_pssm_scores = function(pssms) {
  return ( do.call(c, lapply(pssms, information_content) ))
}

get_metrics = function(y_true, y_pred) {
  # y_true: vector 0's and 1's
  # y_pred: vector contain positions prediction
  
  yt_pos = which(y_true == 1)
  npos = len(yt_pos)
  nneg = len(y_true) - npos
  
  
  tp = fp = 0
  for (i in seq_along(y_pred)) {
    a = y_pred[i]
    found = F
    for (j in seq_along(yt_pos)) {
      b = yt_pos[j]
      if (a %near% b) { found = T; break}
    }
    
    if (found) tp = tp + 1 else fp = fp + 1
  }
  
  fn = npos - tp
  tn = nneg - fp
  
  precison = tp / (tp + fp)
  recall = tp / (tp + fn)
  fscore = 2 * ( (precison * recall) / (precison+recall))
  
  return (list(tp = tp, fp = fp, tn = tn, fn = fn, 
               precison = precison, recall = recall, fscore = fscore))
}
