# Load libraries and user functions
source('libs.R')
source('greedy.R')
source('em.R')
source('r_functions.R')

# load raw dataset
file = 'datasets/GATA2.fa'
X = load_dataset(file)

# Parameters
# w: width of motif
# N: number of lines
# M: number of columns
# L: number of samples per line
# n: number of total samples
w = 14
N = len(X)
M = len(X$seq1)
L = M - w + 1
n = N * L
b = 10

# Sample dataset
Y = get_samples(X, w)

# Theta_1
theta_1 = foreach(i = 1:b) %do% relative_nucleotide_freq(Y)

# theta_2: Greedy inicialization
theta_2 = greedy_inicialization(Y, N, w, b)

# Set priori probabilities
priori = foreach(i = 1:b) %do% c(.9, .1)

tic()
# Main loop
pb <- progress_bar$new(total = 10)
for(i in 1:10) {
  pb$tick()
  
  # Expection
  Z = expectation(Y, theta_1, theta_2, priori)
  rn1 = Z[[1]]
  rn2 = Z[[2]]
  n1 = N1(rn1)
  n2 = N2(rn2)
  
  # Priori maximization
  priori = maximization_priori(Z, n, n1, n2)
  
  # Theta_2 maximization
  theta_2 = maximization_theta_2(Y, rn2, n2)
  
  # Theta_1 maximization
  theta_1 = maximization_theta_1(theta_2)
  
}
toc()

# stopCluster(cl)
