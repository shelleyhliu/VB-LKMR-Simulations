#############################
# Shelley H. Liu 
# Data Simulation File for VB-LKMR
# Simulating gradual quadratic effects model
#############################

set.seed(1)
library(lars); library(mvtnorm); library(SuppDists); library(MCMCpack); library(grplasso); library(magic); library(kernlab); library(MASS); library(fields)

N 				= 100 # Number of subjects
T 				= 4 # Number of time points
P 				= N*T
M				  = 3 # Number of metals

TauG = N #Group size for tau. For the kernel case it should be N.

poly = polydot(degree = 2, offset = 1)

offDiagonal = function(x) {
		diagBelow = diag(x)
		i = 1
		while (i <= TauG) { #changed
			diagBelow=rbind(rep(0,length(x)	+i),cbind(diagBelow,rep(0,length(x) + i - 1)))
			i = i + 1
		}
		mat <- diagBelow + t(diagBelow) - diag(diag(diagBelow))
		return(mat)
}

#############################

cofactor		= 1e-7

ar 			= 0.8

time_mat <- matrix( ar, 4,4 ) #AR-1, so time points 4
diag(time_mat) <- 1
time_mat[1,3] = time_mat[3,1] = time_mat[2,4] = time_mat[4,2] = ar^2
time_mat[1,4] = time_mat[4,1] = ar^3

metal_mat <- matrix( c(1, 0.2, 0.3,0.2, 1, 0.5, 0.3, 0.5, 1), 3,3, byrow=TRUE) 

cmat = kronecker(time_mat, metal_mat)
obs <- matrix( mvrnorm(N, mu=rep(0, 3*4), Sigma=cmat), nrow = N)

Z = obs
Z 		= apply(Z, 2, scale)

Z.1 = Z[,1:3]
Z.2 = Z[,4:6]
Z.3 = Z[,7:9]
Z.4 = Z[,10:12]

X				= cbind(scale(matrix(rnorm(N, mean = 10, sd = 1))), scale(matrix(sample(1:2, N, replace = TRUE))))
beta			= c(1,1)

res.sd = 1 
Y = matrix(rnorm(n=N, mean=0, sd=res.sd)) + X %*% beta + 0.5*(Z.2[,1]^2 - Z.2[,2]^2 + 1/2*Z.2[,1]*Z.2[,2] + Z.2[,1] + Z.2[,2]) + 0.8*(Z.3[,1]^2 - Z.3[,2]^2 + 1/2*Z.3[,1]*Z.3[,2] + Z.3[,1] + Z.3[,2]) + 1*(Z.4[,1]^2 - Z.4[,2]^2 + 1/2*Z.4[,1]*Z.4[,2] + Z.4[,1] + Z.4[,2])
 
# Make W matrix

W 				= diag(N)
counter 		= 1
while(counter < T) {
	W 		= cbind(W, diag(N))
	counter = counter+1
}
