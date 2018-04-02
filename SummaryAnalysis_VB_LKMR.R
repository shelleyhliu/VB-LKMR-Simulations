#############################
# Shelley H. Liu
# Analysis of output
# Nov. 7, 2017
#############################

Summary.MFVB 	= matrix(NA, nrow=T, ncol=4)

# Time point 1
ind 			= 1:N
h_hat.1			= E.q.h[ind]
h_true.1		= rep(0, N)

model.1 = lm(h_hat.1 ~ h_true.1)

Summary.MFVB[1,] = c(summary(model.1)$coef[1,1], 0, 0, sqrt( sum( (h_hat.1 - h_true.1)^2) / N))

# Time point 2
ind 			= (N+1):(2*N)
h_hat.2 		= E.q.h[ind]
h_true.2 		= 0.5*(Z.2[,1]^2 - Z.2[,2]^2 + 1/2*Z.2[,1]*Z.2[,2] + Z.2[,1] + Z.2[,2])

model.2 = lm(h_hat.2 ~ h_true.2)

Summary.MFVB[2,] = c(summary(model.2)$coef[1,1], summary(model.2)$coef[2,1], summary(model.2)$r.sq, sqrt( sum( (h_hat.2 - h_true.2)^2) / N))

# Time point 3
ind 			= (2*N+1):(3*N)
h_hat.3 		= E.q.h[ind]
h_true.3 		= 0.8*(Z.3[,1]^2 - Z.3[,2]^2 + 1/2*Z.3[,1]*Z.3[,2] + Z.3[,1] + Z.3[,2])

model.3 = lm(h_hat.3 ~ h_true.3)

Summary.MFVB[3,] = c(summary(model.3)$coef[1,1], summary(model.3)$coef[2,1], summary(model.3)$r.sq, sqrt( sum( (h_hat.3 - h_true.3)^2) / N))

# Time point 4
ind 			= (3*N+1):(4*N)
h_hat.4 		= E.q.h[ind]
h_true.4 		= 1*(Z.4[,1]^2 - Z.4[,2]^2 + 1/2*Z.4[,1]*Z.4[,2] + Z.4[,1] + Z.4[,2])

model.4 = lm(h_hat.4 ~ h_true.4)

Summary.MFVB[4,] = c(summary(model.4)$coef[1,1], summary(model.4)$coef[2,1], summary(model.4)$r.sq, sqrt( sum( (h_hat.4 - h_true.4)^2) / N))

print("VB: Intercept, Slope, R^2, RMSE")
print( round(Summary.MFVB, 2) )

###################################################

Summary.MCMC	= matrix(NA, nrow=T, ncol=4)

# Time point 1
ind 			= 1:N
h_hat.1			= apply(MCMC$Beta[sel,], 2, mean)[ind]

model.1 = lm(h_hat.1 ~ h_true.1)

Summary.MCMC[1,] = c(summary(model.1)$coef[1,1], 0, 0, sqrt( sum( (h_hat.1 - h_true.1)^2) / N))

# Time point 2
ind 			= (N+1):(2*N)
h_hat.2 		= apply(MCMC$Beta[sel,], 2, mean)[ind]

model.2 = lm(h_hat.2 ~ h_true.2)

Summary.MCMC[2,] = c(summary(model.2)$coef[1,1], summary(model.2)$coef[2,1], summary(model.2)$r.sq, sqrt( sum( (h_hat.2 - h_true.2)^2) / N))

# Time point 3
ind 			= (2*N+1):(3*N)
h_hat.3 		= apply(MCMC$Beta[sel,], 2, mean)[ind]

model.3 = lm(h_hat.3 ~ h_true.3)

Summary.MCMC[3,] = c(summary(model.3)$coef[1,1], summary(model.3)$coef[2,1], summary(model.3)$r.sq, sqrt( sum( (h_hat.3 - h_true.3)^2) / N))

# Time point 4
ind 			= (3*N+1):(4*N)
h_hat.4 		= apply(MCMC$Beta[sel,], 2, mean)[ind]

model.4 = lm(h_hat.4 ~ h_true.4)

Summary.MCMC[4,] = c(summary(model.4)$coef[1,1], summary(model.4)$coef[2,1], summary(model.4)$r.sq, sqrt( sum( (h_hat.4 - h_true.4)^2) / N))

print("MCMC: Intercept, Slope, R^2, RMSE")
print(round(Summary.MCMC, 2))

