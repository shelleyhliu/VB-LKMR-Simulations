#############################
# Shelley H. Liu
# Mean Field Variational Bayes for Lagged Kernel Machine Regression
#############################

# Set flag for doing kernel approach:

doKernelApproach = TRUE

# Set data simulation

source("DataSimulationCode_VB_LKMR.R")

# Set tasks:

doMFVB 			= TRUE
doMCMC 			= TRUE
doPrelimPlots 	= FALSE
plotTRUE		= FALSE
doAnalysis		= FALSE

ptm 			= proc.time()

# Set constant matrices needed for MFVB:

XTX 			= crossprod(X)
WTW 			= crossprod(W)

# Set hyperparameters
r1 				= 60
r2 				= 45
delta1 		= 10
delta2 		= 10

if (doMFVB)
{
	ptm = proc.time()

   # Do iter_num iterations of MFVB updating:

   iter_num 			= 100

   E.q.recip.sigsq 		= 1
   E.q.beta 			= c(1,1)
   E.q.recip.tausq 		= rep(1, T)
   E.q.recip.omegasq 	= rep(1, T-1)
   E.q.lambda1sq 		= 1
   E.q.lambda2sq		= 1
   
   A.q.recip.tausq 		= rep(1, T)
   B.q.recip.tausq 		= rep(1, T)
   A.q.recip.omegasq 	= rep(1, T-1)
   B.q.recip.omegasq 	= rep(1, T-1)
   
   sig.shape = 5 
   sig.scale = 5
   a = sig.shape
   gamma = sig.scale
   
   logpyvec = c()
   logpyvec[1] = 0
   
   list.G 				= list()
   for (g in 1:T) {
   		if (!doKernelApproach) {  
       		list.G[[g]] 	= 1
   		} 
   		if (doKernelApproach) { 
			list.G[[g]]		= solve(kernelMatrix(poly, Z[,(M*g-(M-1)):(M*g)]) + cofactor*diag(N))
   		}  
   	}
 
	for (itnum in 1:iter_num) {
	
	list.GT = list()
	for (g in 1:T) {
		list.GT[[g]] = list.G[[g]] * E.q.recip.tausq[g]
	}

	cov.1 				= do.call(adiag, list.GT)
	
	# 2. Diagonal

	list.sig.sqf 		= c()
	list.sig.sqf[1] 	= E.q.recip.omegasq[1]
	list.sig.sqf[T] 	= E.q.recip.omegasq[T-1]
	for (g in 2:(T - 1)) {
		list.sig.sqf[g] = E.q.recip.omegasq[g] + E.q.recip.omegasq[g-1]
	}

	list.mat.sig.sqf 	= list()
	for (g in 1:T) {
			list.mat.sig.sqf[[g]] = list.sig.sqf[g] * diag(TauG) #changed
	}

	cov.2 				= do.call(adiag, list.mat.sig.sqf)

	#3. Off-diagonals
	cov.3 				= offDiagonal((-1)*(rep(E.q.recip.omegasq, each=TauG))) #changed

	    E.q.Sigma.h.inverse = cov.1 + cov.2 + cov.3 
   
	Cov.q.h 			= solve( E.q.recip.sigsq * WTW + E.q.Sigma.h.inverse + cofactor*diag(P)) #changed
	E.q.h				= E.q.recip.sigsq * Cov.q.h %*% t(W) %*% (Y - X %*% E.q.beta)
	E.q.beta 			= solve(XTX) %*% t(X) %*% (Y - W %*% E.q.h)


        B.q.recip.sigsq = as.numeric(crossprod(Y - X %*% E.q.beta - W %*% E.q.h)/2) + sig.scale #+ t(E.q.h) %*% E.q.Sigma.h.inverse %*% E.q.h/2) + sig.scale #changed
               
        if (doKernelApproach) {
        	A.q.recip.sigsq = N/2 + sig.shape 
        }

	E.q.recip.sigsq 	= A.q.recip.sigsq/B.q.recip.sigsq

  list.h 				= list()
	for (g in 1:T) {

                if (!doKernelApproach)
                {
		    list.h[[g]] 	= E.q.h[(1+(g-1)):(g)]

                } 
                if (doKernelApproach)
                {  
                   list.h[[g]] 		= E.q.h[(1+N*(g-1)):(N*g)]
                }
                   
		A.q.recip.tausq[g] 	= sqrt( E.q.lambda1sq / (t(list.h[[g]]) %*% list.G[[g]] %*% list.h[[g]]) ) #changed
		B.q.recip.tausq[g] 	= E.q.lambda1sq
		E.q.recip.tausq[g] 	= A.q.recip.tausq[g]
	}	
	for (g in 1:(T-1)) {
		A.q.recip.omegasq[g] = sqrt( E.q.lambda2sq / ( sum((list.h[[g+1]] - list.h[[g]])^2)) )
		B.q.recip.omegasq[g] = E.q.lambda2sq
		E.q.recip.omegasq[g] = A.q.recip.omegasq[g]
	}
		
	if (doKernelApproach) {
		A.q.lambda1sq = T*(N+1)/2 + r1
	}
    B.q.lambda1sq = 1/2*sum(sapply(E.q.recip.tausq, function(x) 1/x)) + delta1
    E.q.lambda1sq = A.q.lambda1sq/B.q.lambda1sq
   
    A.q.lambda2sq = (T - 1 + r2)
    B.q.lambda2sq = (1/2*sum(sapply(E.q.recip.omegasq, function(x) 1/x)) + delta2)
	E.q.lambda2sq = A.q.lambda2sq / B.q.lambda2sq 
	
   }

  VB.runtime	= proc.time()[3] - ptm[3]
}

if (doMCMC)
{
       if (doKernelApproach)
       {
        num.reps 	= 10000
	  		sel 		= 5000:num.reps
	  		ptm			= proc.time()
	  		source("RunMCMC_VB_LKMR.R")
	  		MCMC.runtime	= proc.time()[3] - ptm[3]
       }  
}
