##
## hmm.R

##
## 1. sample generation
##

## state names
snames <- c("Fair", "Bias")
K <- length(snames)

## output symbols
onames <- c("Head", "Tail")
R <- length(onames)

## initial state
p <- c(0.5, 0.5)
logp <- log(p)

## transition[state j, state k]
a <- matrix(c(0.7, 0.3,
              0.9, 0.1), ncol=2, byrow=TRUE)
loga <- log(a)

## emission[state k, observation m]
b <- matrix(c(0.2, 0.8,
              0.5, 0.5), ncol=2, byrow=TRUE)
logb <- log(b)

## 
## Generate(N, p, a, b)
##
gen <- function(size, p, a, b) {
    N <- size
    K <- nrow(b)
    M <- ncol(b)
    zs <- numeric(N)
    ys <- numeric(N)

    zs[1] <- sample(1:K, size=1, prob=p)
    ys[1] <- sample(1:M, size=1, prob=b[zs[1], ])
    for(n in 2:N) {
        zs[n] <- sample(1:K, size=1, prob=a[zs[n-1], ])
        ys[n] <- sample(1:M, size=1, prob=b[zs[n], ])    
    }

    list(states=zs, symbols=ys)
}


N <- 100
sam <-  gen(N, p, a, b)
z <- sam$state
y <- sam$symbol


##
## compute the complete joint prob : P(y,z|p,a,b)
##

loglik_joint <- function(y, z, p, a, b) {
      N <- length(y)
      logprobs <- numeric(N)
      logprobs[1] <- logp[z[1]] + logb[z[1], y[1]]
      for(n in 2:N) {
          logprobs[n] <- loga[z[n-1],z[n]] + logb[z[n], y[n]]
      }
      sum(logprobs)
}



loglik_joint(y, z, p, a, b)



##
## 2. forward-backwrad
##

p <- c(0.5, 0.5)
a <- matrix(c(0.7, 0.3,
              0.9, 0.1), ncol=2, byrow=TRUE)
b <- matrix(c(0.2, 0.8,
              0.5, 0.5), ncol=2, byrow=TRUE)
y <- c(1, 1, 2, 1, 1, 2, 1, 1, 1, 2)


## alpha[k, i] is stored at (s[k], i) in the trellis and expresses the total
## probability of ending up in state s[k] at time i
##
alpha <- matrix(0, nrow=N, ncol=K)
alpha[1,] <- p * b[,y[1]]
for(n in 2:N) {
    alpha[n,] <- colSums(alpha[n-1,] * a) * b[,y[n]]
}
log(sum(alpha[N,]))



forward <- function(y, p, a, b) {
    N <- length(y)
    K <- nrow(a)
    logp <- log(p)
    loga <- log(a)
    logb <- log(b)
    logalpha <- matrix(0, nrow=N, ncol=K)
    logalpha[1,] <- logp + logb[,y[1]]
    for(n in 2:N) {
        mm <- max(logalpha[n-1,])
        logalpha[n,] <- log(colSums(exp(logalpha[n-1,] - mm + loga))) + mm  + logb[,y[n]]
    }
    logalpha
}


beta <- matrix(0, nrow=N, ncol=K)
beta[N,] <- 1
for(n in (N-1):1) {
    beta[n,] <- colSums(t(a) * b[,y[n+1]] * beta[n+1,])
}
log(  sum(p * b[,y[1]] * beta[1,]))

backward <- function(y, p, a, b) {
      N <- length(y)
      K <- nrow(a)
      logp <- log(p)
      loga <- log(a)
      logb <- log(b)
      logbeta <- matrix(0, nrow=N, ncol=K)
      logbeta[N,] <- 0
      for(n in (N-1):1) {
          ##
          ## 
          ##
      }
      logbeta    
}

logsum <- function(x) {
    mx <- max(x)
    log(sum(exp(x - mx))) + mx
}

##
## compute the marginal prob : log P(y|p,a,b)
##
loglik_marginal <- function(y, p, a, b) {
    logalpha <- forward(y, p, a, b)
    N <- nrow(logalpha)
    logsum(logalpha[N,])
}

y <- c(1, 1, 2, 1, 1, 2, 1, 1, 1, 2)
loglik_marginal (y, p, a, b)

## or
logbeta <- backward(y, p, a, b)
logsum(logp + logb[,y[1]] + logbeta[1, ])

##
## 3. gamma
##
  logalpha <- forward(y, p, a, b)
  logbeta <- backward(y, p, a, b)

  temp <- logalpha + logbeta
  mx <- apply(temp, MARGIN=1, max)
  gamma <- prop.table(exp(temp - mx), margin=1)

##
## 4. viterbi
##
viterbi <- function(y, a, b, p) {
    N <- length(y)
    K <- nrow(a)
    loga <- log(a)
    logb <- log(b)
    logp <- log(p)
    
    logdelta <- matrix(0, nrow=N, ncol=K)
    psi <- matrix(0, nrow=N, ncol=K)

    ## Initialize
    logdelta[1,] <- logp + logb[,y[1]]
    psi[1,] <- 0

    ## Recursion
    for(n in 2:N) {
        ##for(k in 1:K) {
        ##  logdelta[n,k] <- max(logdelta[n-1,] + loga[,k]) + logb[k,y[n]]
        ##  psi[n,k] <- which.max(logdelta[n-1,] + loga[,k])
        ##}
        logdelta[n, ] <- apply(t(logdelta[n-1,] + loga) + logb[,y[n]], MARGIN=1, FUN=max)
        psi[n, ] <- apply(logdelta[n-1,] + loga, MARGIN=2, FUN=which.max)
    }

    ## Termination
    zstar <- numeric(N)
    zstar[N] <- which.max(logdelta[N,])
    logPzstar <- max(logdelta[N,])

    ## Backtracking
    for(n in (N-1):1) {
        zstar[n] <- psi[n+1,zstar[n+1]]
    }

    list(path=zstar, loglik=logPzstar)
}

viterbi(y, a, b, p)

##
## 5. implement Baum-Welch algorithm
##
## y : observed symbols
## p, a, b : initial values
##
baumwelch <- function(y, p, a, b, maxiter=100)  {
    ##
    ##
    ##
}





y <- c(1, 1, 2, 1, 1, 2, 1, 1, 1, 2)
K <- 2
R <- max(y)

a0 <- prop.table(matrix(runif(K*K), nrow=K, ncol=K), margin=1)
b0 <- prop.table(matrix(runif(K*R), nrow=K, ncol=R), margin=1)
p0 <- rep(1/K, K)
ans <- baumwelch(y, p0, a0, b0, maxiter=500L)
ans$loglik
plot(ans$logliks[-1])



