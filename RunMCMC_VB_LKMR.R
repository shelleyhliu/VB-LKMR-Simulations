#############################
# Shelley H. Liu 
# Run MCMC-LKMR
# Nov. 7, 2017
#############################

sig.shape 	= 5 
sig.scale 	= 5 

lambda1.sq 	= rgamma(1,shape = r1, rate = delta1) 
lambda2.sq 	= rgamma(1,shape = r2, rate = delta2) 
sig.sq0		= rinvgamma(1, shape = sig.shape, scale = sig.scale)
mean.bef 	= rep(0, P)
tau.sqf  	= rgamma(T, shape = (N + 1)/2, rate = lambda1.sq/2)
sig.sqf  	= rexp(T-1, rate = lambda2.sq/2)
c 			= length(beta)
conf 		= rep(1, c)
U 			= X

# 1. Block diagonal matrix 
list.G 		= list()
poly 		= polydot(degree=2, offset=1)
for (g in 1:T) {
	list.G[[g]] = solve(kernelMatrix(poly, Z[,(M*g-(M-1)):(M*g)]) + cofactor*diag(N)) / tau.sqf[g]
	}  
cov.bf1 	= do.call(adiag, list.G)

# 2. Diagonal

list.sig.sqf 		= c()
list.sig.sqf[1] 	= 1/sig.sqf[1]
list.sig.sqf[T] 	= 1/sig.sqf[(T-1)]
for (g in 2:(T - 1)) {
	list.sig.sqf[g] = 1/sig.sqf[g] + 1/sig.sqf[g-1]
}

list.mat.sig.sqf 	= list()
for (g in 1:T) {
	list.mat.sig.sqf[[g]] = list.sig.sqf[g] * diag(N)
}

cov.bf2 			= do.call(adiag, list.mat.sig.sqf)

#3. Off-diagonals
new.sig.sqf 		= rep(1/sig.sqf, each=TauG)
cov.bf3 			= offDiagonal((-1)*new.sig.sqf)
    
# 4. Beta covariance matrix
cov.bf 				= cov.bf1 + cov.bf2 + cov.bf3

# A few ways to invert cov.bf
SIG       			= solve(cov.bf + cofactor*diag(P))

beta.fp    			= rmvnorm(1,mean=mean.bef,sigma=SIG)     
    
# FOR POSTERIOR
sigsq0.post = lambda1f.post = lambda2f.post = NULL 
beta.f      = rbind( beta.fp,matrix(rep(NA,num.reps*P),ncol=P) )
tausqf.post = rbind( tau.sqf,matrix(rep(NA,num.reps*T),ncol=T) )
sigsqf.post = rbind( sig.sqf,matrix(rep(NA,num.reps*(T-1)),ncol=T-1) )
conf.post 	= rbind( conf, matrix(rep(NA,num.reps*c),nrow=num.reps) ) 
conf 		= matrix(conf, nrow = c)

######################################################    
# POSTERIOR 
cat(c("Job started at:",date()),fill=TRUE)
source("MCMCUpdates_VB_LKMR.R")
cat(c("Job finished at:",date()),fill=TRUE)

MCMC = list(Sigsq = sigsq0.post, Lam1 = lambda1f.post, Lam2 = lambda2f.post, Beta = beta.f, Tau = tausqf.post, Omega = sigsqf.post, Conf = conf.post)
#save(MCMC, file=paste0("MCMC", jobid, ".RData"))

