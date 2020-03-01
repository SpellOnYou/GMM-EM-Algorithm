##
## K-Component Gaussian Mixture
##
## See Mixture Models
##

## Model
K <- 3

p <- c(0.5, 0.4, 0.1)
m <- c(0, 4, 10)
s <- c(1, sqrt(2), 1)


## Generation
N <- 10000
z <- sample(1:3, N, replace=TRUE, prob=p)
y <- numeric(N)
for(i in 1:N) y[i] <- rnorm(1, m[z[i]], s[z[i]])

hist(y, breaks=sqrt(length(y)))

## complete MLE
prop.table(table(z))
z1 <- ifelse(z == 1, 1, 0)
z2 <- ifelse(z == 2, 1, 0)
z3 <- ifelse(z == 3, 1, 0)

m1 <- sum(y*z1) / sum(z1)
m2 <- sum(y*z2) / sum(z2)
m3 <- sum(y*z3) / sum(z3)


##----------------------------------------------------------------
## EM. latent z.
##

maxiter <- 500L
tol <- 1e-7

K <- 3
N <- length(y)

## initialize
p.old <- p.new <- rep(1/K, K)
m.old <- m.new <- sample(y, K)
s.old <- s.new <- rep(sd(y), K)

for(r in 1:maxiter) {
    ## E-step
    h <- matrix(0, ncol=K, nrow=N)
    for(j in 1:K) h[,j] <- p.old[j] * dnorm(y, m.old[j], s.old[j])
    h <- prop.table(h, margin=1)

    ## M-step
    for(j in 1:K) {
        nj <- sum(h[,j])
        p.new[j] <- nj / N
        m.new[j] <- sum(h[,j]*y) / nj
        s.new[j] <- sqrt( sum(h[,j]*(y-m.new[j])^2) / nj )
    }

    ## Stop Condition
    change <- sum(abs(1-p.old/p.new) + abs(1-m.old/m.new) + abs(1-s.old/s.new))
    if(change < tol) break
    p.old <- p.new
    m.old <- m.new
    s.old <- s.new
}
 


