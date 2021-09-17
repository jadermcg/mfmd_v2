source('libs.R')
source('greedy.R')
source('eda.R')
source('em.R')
source('results.R')

# load raw dataset
file_name = 'datasets/MA1235.1.fa'
X = load_dataset(file_name)

# Parameters
# w: width of motif
# N: number of lines
# M: number of columns
# L: number of samples per line
# n: number of total samples
w = 11
N = len(X)
M = len(X$`1`)
L = M - w + 1
n = N * L
beta  = sqrt(L)
#beta = L
# Sample dataset
Y = get_samples(X, w)

# Theta_1
theta_1 = (foreach(i = 1:beta) %do% relative_nucleotide_freq(Y))

# theta_2: Greedy inicialization
theta_2 = greedy_inicialization(Y, N, w, beta)

# Set priori probabilities
priori = priori = (foreach(i = 1:beta) %do% c(.95, .05))

# Join theta parameters
theta = mapply(function(a,b,c) list(theta_1 = a, theta_2 = b, priori = c), a = theta_1, b = theta_2, c = priori, SIMPLIFY = F)[1:3]

# Run EDA
new_theta =  lapply(theta, function(x) eda(Y, x, 20, 3))

# Analysis
rn = lapply(new_theta, function(t) responsibilities(Y, w, t))
POS = lapply(rn, function(r) get_pos(r$rn1, r$rn2, N, L))
seq = rep(0, L)
seq[51] = 1
y_true = rep(seq, N)
y_pred = which( rn[[1]]$rn2 > rn[[1]]$rn1 )

metrics = get_metrics(y_true, y_pred)

