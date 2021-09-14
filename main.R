# Load libraries and user functions
source('libs.R')
source('greedy.R')
source('em.R')
source('r_functions.R')
source('results.R')

# load raw dataset
file_name = 'datasets/MA1235.1.fa'
X = load_dataset(file_name)
X = X[sample(1:len(X), len(X))]

# Parameters
# w: width of motif
# N: number of lines
# M: number of columns
# L: number of samples per line
# n: number of total samples
tic()
w = 11
N = len(X)
M = len(X$`1`)
L = M - w + 1
#n = N * L
beta  = L
# stop_loop = round( log(N) * log(L) )
stop_em_loop = 50
# Sample dataset
Y = get_samples(X, w)

# Theta_1
theta_1 = (foreach(i = 1:beta) %do% relative_nucleotide_freq(Y))
#theta_1 = (foreach(i = 1:beta) %do% rep(.25, 4))

# theta_2: Greedy inicialization
theta_2 = greedy_inicialization(Y[1:200], 200, w, beta)

# Set priori probabilities
priori = (foreach(i = 1:beta) %do% c(.5, .5))


# Expectation maximization
params = em(Y[1:200], theta_1, theta_2, priori, stop_em_loop)

theta_1 = params$theta_1
theta_2 = params$theta_2
priori = params$priori
rn1 = params$rn1
rn2 = params$rn2
inf_cont = params$inf_cont

pssms = get_pssms(theta_1, theta_2)
scores = get_pssm_scores(pssms)
best = which.max(scores)
pos = get_pos(rn1, rn2, L)
print(pos[[best]])
fprintf('%s -> %f', "O maior score é:", max(scores))
toc()


# stopCluster(cl)
