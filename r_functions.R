"%==%" = function(x,y) { if(all.equal(target = x, current = y, tolerance=1e-5) == T) { T } else { F } }

load_dataset = function(file_name) {
  X = read.fasta(file = file_name, seqtype = 'DNA', as.string = F)
  return ( lapply(X, swap_letter_to_number) )
}

get_samples = function(X, w) {
  Y = w_mers(X, w)
  Y = lapply(Y, function(x) as.list(data.frame(t(x))))
  return (Y)
}

swap_letter_to_number = function (x) {
  x[which(x == 'a')] = 1
  x[which(x == 'c')] = 2
  x[which(x == 'g')] = 3
  x[which(x == 't')] = 4
  return (as.integer(x))
}


w_mers = function(X, w) {
  return ( lapply(X, function(x) w_mersC(x,w)) )
}

complement = function(x) {
  if (x == 1) return (4)
  if (x == 4) return (1)
  if (x == 2) return (3)
  if (x == 3) return (2)
}

information_content = function(pssm) {
  return (sum(pssm))
}


count_nucleotides = function(Y) {
  w = len(Y[[1]][[1]])
  
  return ( colSums(do.call(rbind, lapply(Y, function(y) count_nucleotidesC(y,w)))) + 1 )
  
}

relative_nucleotide_freq = function(Y) {
  n = len(Y$seq1) * len(Y)
  w = len(Y[[1]][[1]])
  
  return( count_nucleotides(Y) / (n*w+4) )
}

cartesian_prod = function(P, Q) {
  return (lapply(P, function(x) lapply(Q, function(y) rbind(x,y))))
}

likelihood_t1 = function(Y, w, theta, priori) {
  return ( lapply(Y, function(y) likelihood_t1C(y, w, theta, priori)) )
}

likelihood_t2 = function(Y, w, theta, priori) {
  return ( lapply(Y, function(y) likelihood_t2C(y, w, theta, priori)) )
}

