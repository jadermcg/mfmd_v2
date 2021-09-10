# Load libraries and user functions
source('libs.R')
source('greedy.R')
source('em.R')
source('r_functions.R')
source('results.R')

# load raw dataset
file = 'datasets/CRP.fa'
X = load_dataset(file)

# Parameters
# w: width of motif
# N: number of lines
# M: number of columns
# L: number of samples per line
# n: number of total samples
tic()
w = 22
N = len(X)
M = len(X$seq1)
L = M - w + 1
n = N * L
beta = round (log(L) * sqrt(L))
beta = L

# Sample dataset
Y = get_samples(X, w)

# Theta_1
theta_1 = foreach(i = 1:beta) %do% relative_nucleotide_freq(Y)

# theta_2: Greedy inicialization
theta_2 = greedy_inicialization(Y, N, w, beta)

# Set priori probabilities
priori = foreach(i = 1:beta) %do% c(.9, .1)


# Expectation maximization
params = em(Y, theta_1, theta_2, priori)

theta_1 = params$theta_1
theta_2 = params$theta_2
priori = params$priori
rn1 = params$rn1
rn2 = params$rn2
log_likes = params$log_likes
inf_cont = params$inf_cont

toc()

# Result analysis
pos = get_pos(rn1, rn2, L)
pssms = get_pssms(theta_1, theta_2)
pssm_scores = get_pssm_scores(pssms)
best = which.max(pssm_scores)
pssm = pssms[[best]]
theta_1 = theta_1[[best]]
theta_2 = theta_2[[best]]
priori = priori[[best]]
sample_scores = score_samples(Y, pssm, w)
g = factor(do.call(c,lapply(1:N, function(x) rep(x, L))))
sample_scores = split(sample_scores, g)
pos = pos[[best]]
log_likes = do.call(c, lapply(log_likes , function(x) x[[best]]))
inf_cont = do.call(c, lapply(inf_cont , function(x) x[[best]]))


# Plots
plot(log_likes, type = 'l', xlab = 'EM Iteration', ylab = 'Log likelihood', 
     main = 'Log-likelihood evolution')

plot(inf_cont, type = 'l', xlab = 'EM Iteration', ylab = 'Information content', 
     main = 'Information content evolution')

# stopCluster(cl)
