#############################
# Shelley H. Liu 
# Updates for MCMC-LKMR
# Nov. 7, 2017
#############################

for (iter in 1:num.reps)  {
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

	list.mat.sig.sqf	= list()
	for (g in 1:T) {
		list.mat.sig.sqf[[g]] = list.sig.sqf[g] * diag(N)
	}

	cov.bf2 			= do.call(adiag, list.mat.sig.sqf)

	#3. Off-diagonals
	new.sig.sqf 	= rep(1/sig.sqf, each=TauG)
	cov.bf3 			= offDiagonal((-1)*new.sig.sqf)
    
	# 4. Beta covariance matrix inverse
	cov.bf 			= cov.bf1 + cov.bf2 + cov.bf3

    cov.bef      	=  solve(1/sig.sq0*WTW+cov.bf+cofactor*diag(P))
    mean.bef     	= 1/sig.sq0 * cov.bef%*%t(W)%*%(Y - U%*%conf)
    
    beta.f[iter+1,] = rmvnorm(1,mean=mean.bef,sigma=cov.bef)

    # sig.sq0
    sh.sig      	= N/2 + sig.shape
    sc.sig      	= 1/2*t(Y-W%*%beta.f[iter+1,] - U%*%conf)%*%(Y-W%*%beta.f[iter+1,] - U%*%conf) + sig.scale  
    
    sig.sq0     	= rinvgamma(1, shape=sh.sig, scale=sc.sig)
    sigsq0.post 	= c(sigsq0.post, sig.sq0)
     
    gam <- c()
    for (j in 1:T){
    	term.beta 	= as.matrix(beta.f[iter+1, ((j-1)*N+1):(j*N) ])
    	term.gam 	= t(term.beta) %*% solve(kernelMatrix(poly, Z[,(M*j-(M-1)):(M*j)]) + cofactor*diag(N)) %*% term.beta 
    	gam[j]  	= rinvGauss(1, nu=sqrt(lambda1.sq/term.gam), lambda=lambda1.sq)          		
    	tau.sqf[j] 	= 1/gam[j] 
    }  	
    tausqf.post[iter+1,] = gam 
    
    et 					 = c()
    for (k in 1:(T-1)){
    	nu.k    		 = sqrt(lambda2.sq/sum((beta.f[iter+1,((k-1)*N+1):(k*N)]-beta.f[iter+1,(k*N+1):((k+1)*N)])^2))
    	et[k]  			 = rinvGauss(1, nu=nu.k, lambda=lambda2.sq)
    	sig.sqf[k] 		 = 1/et[k] 
      	}
    sigsqf.post[iter+1,] = et 
	    
    # lambda
    sh.lam1       = T*(N+1)/2 + r1 
    sc.lam1       = 1/2*sum(tau.sqf) + delta1 
    lambda1.sq    = rgamma(1, shape=sh.lam1, rate=sc.lam1)
    lambda1f.post = c(lambda1f.post, lambda1.sq)

    sh.lam2       = T - 1 + r2 
    sc.lam2       = 1/2*sum(sig.sqf) + delta2 
    lambda2.sq    = rgamma(1, shape=sh.lam2, rate=sc.lam2)
    lambda2f.post = c(lambda2f.post, lambda2.sq) 
    
    mean.conf 	= solve(t(U) %*% U) %*% t(U) %*% (Y - W%*%beta.f[iter+1,])
    sig.conf 	= sig.sq0 * solve(t(U) %*% U) 
    conf 		= rmvnorm(1,mean=mean.conf,sigma=sig.conf)
    conf.post[iter+1,] = conf  
    conf 		= matrix(conf, nrow=c)
}
