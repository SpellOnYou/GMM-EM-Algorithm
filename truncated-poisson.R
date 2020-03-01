# log-likelihood of Truncated Poisson distribution

input = "7, 5, 6, 3, 5, 5, 4, 5, 4, 4, 4, 4, 4, 4, 3, 5, 4, 4, 4, 4, 5, 6, 6, 4, 7, 5, 5, 5, 4, 9, 7, 3, 3, 5, 5, 6, 4, 3, 4, 7, 4, 3, 5, 5, 5, 5, 4, 3, 5, 5, 3, 5, 4, 6, 2, 6, 6, 5, 3, 4, 4, 4, 8, 5, 7, 3, 8, 2, 4, 7, 9, 3, 6, 6, 6, 6, 4, 5, 4, 3, 6, 3, 4, 4, 6, 3, 7, 4, 5, 2, 6, 9, 4, 6, 12, 5, 4, 3, 10, 4, 8, 8, 8, 3, 4, 5"
data = strsplit(input, ", " )[[1]]
data = as.numeric(data); data

# For the fisher scoring, the Hessian is replaced by expected information. And information is given.

poi.pack <- function(x) {
	n <- length(x)
	xbar <- mean(x)
	log.lik <- function(lambda) n * xbar * log(lambda) - n*lambda - sum(lfactorial(x)) - n*log(1-exp(-lambda))
	neg.lik = function(lambda) {a = log.lik(lambda); -a} #to nlp
	score <- function(lambda) n * xbar / lambda - n - n*(exp(-lambda)) / (1-exp(-lamb
da))
	# hess <- function(p) # we don't need this function
	info <- function(lambda) + n / lambda - n*exp(-lambda)/((1-exp(-lambda))^2)
	list(ll=log.lik, score=score, info=info, neg = neg.lik)
}

fn = poi.pack(data)

# Newton-raphson and Fisher scoring

nr <- function(S, I, l0, maxiter = 100, tol = 1e-7) {
	ls <- numeric(maxiter)
	for(i in 1:maxiter) {
		l1 <- l0 + solve(I(l0)/S(l0)) # sign should be changed
		if(abs(l1 - l0) < tol ) break
		l0 <- l1
		ls[i] <- l0
	}
	if(i == maxiter) warning("maximum iteration!")
	list(value=l1, iter=i, trace=ls[1:(i-1)])
}

out = nr(fn$score, fn$info, mean(data)); out

# compare it with R-built in algorithm, nlm.

out2 = nlm(fn$neg, mean(data), hessian=TRUE); out2

cat("My algorithm : ", out$value, "nlm estimation : ",out2$estimate)

## My algorithm : 4.92627 nlm estimation : 4.926268

# compare ALL of the value by putting the value to the score function.

a = abs(fn$score(mean(data))) # mean(lambda)
b = abs(fn$score(out$value)) # our algorithm using fisher and NR
c = abs(fn$score(out2$estimate)) # embeded nlm function
b<c; b<a; c<a;

## The algorithm by myself (not usign library) outs figure which is most close to 0

# plot this result



a = abs(fn$score(mean(data))) # mean(lambda)
b = abs(fn$score(out$value)) # our algorithm using fisher and NR
c = abs(fn$score(out2$estimate)) # embeded nlm function
b<c; b<a; c<a;