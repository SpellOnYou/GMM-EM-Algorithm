# Newton-Raphson method for fitting Gamma distribution.

hessian = function(y){	
	n = length(y)
	hessian = function(params){
		h_aa = -1*n*trigamma(params[1])
		h_ab = n/params[2]
		h_bb = -1*n*params[1]/params[2]^2
		h = matrix(c(h_aa,h_ab,h_ab,h_bb),ncol=2)
		h
	}
}
score = function(y){
	n = length(y)
	score = function(params){
		s_a = n*log(params[2])+sum(log(y))-n*digamma(params[1])
		s_b = n*params[1]/params[2]-sum(y)
		s = c(s_a,s_b)
		s
	}
}

NRroot <- function(S, H, params, maxiter=100, tol=1e-22) {
	xs <- matrix(0, ncol=length(params), nrow=maxiter)
	for(i in 1:maxiter) {
		grad = solve(H(params))%*% S(params)
		cat("\nparams : ", params)
		params_1 <- params - grad #element wise minus?
		cat("\nparams_1 : ", params_1)
		if(sum(abs(params_1 - params)/ params_1) < tol) break
		xs[i,] <- params[1]
	}
	ans <- list(val=params, iter=i, trace=xs[1:(i-1),])
	ans
}

# Data initialization
set.seed(0821); alpha=2; beta=3
data <- round(rgamma(50,shape=alpha, scale=beta),2)
theta =c(alpha, beta); theta

# Making instance of score, hessian
my_score = score(data); my_hessian = hessian(data)
my_score(theta); my_hessian(theta)

#moments alpha, moments beta to initial value of parameter
moments_alpha <- mean(data)^2/var(data); moments_alpha

moments_beta <- var(data)/mean(data); moments_beta

params.0 = c(moments_alpha, moments_beta); params.0

# solve function test
solve(my_hessian(params.0))%*% my_score(params.0)

NRroot(my_score, my_hessian, params.0)


