# Load libraries and user functions
source('libs.R')
source('greedy.R')
source('em.R')
source('results.R')

# load raw dataset
file_name = 'datasets/CRP.fa'
X = load_dataset(file_name)

# Parameters
# w: width of motif
# N: number of lines
# M: number of columns
# L: number of samples per line
# n: number of total samples
tic()
w = 22
N = len(X)
M = len(X$`1`)
L = M - w + 1
n = N * L
beta  = sqrt(L)
stop_loop = 15
# Sample dataset
Y = get_samples(X, w)

# Theta_1
theta_1 = (foreach(i = 1:beta) %do% relative_nucleotide_freq(Y))

# theta_2: Greedy inicialization
theta_2 = greedy_inicialization(Y, N, w, beta)

# Set priori probabilities
priori = (foreach(i = 1:beta) %do% c(1 - (1/L), 1/L))

# Thetas
thetas = mapply(function(t1, t2, p) list(theta_1 = t1, theta_2 = t2, priori = p), 
                t1 = theta_1, t2 = theta_2, p = priori, SIMPLIFY = F)

theta = thetas[[1]]
fprintf('IC before = %.4f\n', information_content(pssm_matrixC(theta$theta_1, theta$theta_2)))
fprintf('LL before = %.4f\n\n', log_likelihood(Y, theta, w)$log_like)

# Expectation maximization
params = mapply(function(t) em(Y = Y, theta = t, stop_loop = stop_loop), 
                t = thetas, SIMPLIFY = F)


# Analysis
best = order(sapply(params, function(x) information_content(
                            pssm_matrixC(x$theta$theta_1, x$theta$theta_2))), decreasing = T)[1]

theta = params[[best]]$theta
E = log_likelihood(Y, theta, w)
pos = cbind(which(E$rn2 > E$rn1) %/% L + 1, which(E$rn2 > E$rn1) %% L)
par(mfrow = c(2,1))
plot(params[[best]]$log_like, type='l', ylab = 'Likelihood', xlab = 'Number of iterations')
plot(params[[best]]$inf_cont, type='l', ylab = 'Log odds', xlab = 'Number of iterations')
fprintf('IC after = %.4f\n', information_content(pssm_matrixC(theta$theta_1, theta$theta_2)))
fprintf('LL after = %.4f\n\n', log_likelihood(Y, theta, w)$log_like)


toc()


# stopCluster(cl)
