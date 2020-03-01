##----
## EM clustering Sample Data
##
doctermfreq <- read.table("emcluster_sample.txt", header=TRUE)

docid <- doctermfreq$docid
termid <- doctermfreq$termid

##
## true class label data
##
##
docid_class <- unique(doctermfreq[, c("docid", "clabel")])
doc_class_mat <- xtabs(~ docid + clabel, docid_class)
term_class_mat <- xtabs(~ termid + clabel, data=doctermfreq)

##
## utility functions
##
colsums <- function(x) {
    if(is.matrix(x)) colSums(x)
    else x
}

## log(sum(probs)) = log(sum(exp(log(probs))))
##                 = log(sum(exp( log(probs) - max(log(probs)) ))) + max(log(probs))
logsumexp <- function(logps)  log(sum(exp(logps - max(logps)))) + max(logps)

##
## 
##
N <- max(docid)        # number of documents
M <- max(termid)       # number of terms
K <- 2                 # number of clusters
R <- 50                # number of repetitions of EM

loglik <- numeric(R)   # keep track of the log likelihood at each step

h <- matrix(0, nrow=N, ncol=K)               # P(class k | doc n)
logprior <- log(rep(1/K, K))                 # log P(class k)
cond1 <- matrix(runif(M*K), nrow=M, ncol=K)  # P(term m = 1 | class k)
logcond1 <- log(cond1)                       # log P(term m = 1 | class k)
logcond0 <- log(1-cond1)                     # log P(term m = 0 | class k)


for (r in 1:R) {
    ##
    ## E-step: compute h[i, k] = P(cluster k | document i)
    ##
    for(i in 1:N) {
        ## By Bayes Theorem,
        ##
        ## P(cluster k | doc i) = P(doc i | cluster k) * P(cluster k)  / P(doc i)
        ##
        ## At first, ignore the normalizing constant P(doc i)
        ##         P(doc i)  = sum of P(doc i | cluster k) * P(cluster k) for all k
        ##
        ## The priror probability P(cluster k) => use logprior[k]
        ##
        ## The conditional probability P(doc i | cluster k) => Naive Bayes Assumption:
        ##
        ## (1) compute only terms which occur in document i
        ##
        ## log P(doc i | cluster k) =  sum log P(term m = 1| cluster k) for only term m in document i.
        ##
        ## (2) compute all terms
        ##
        ## log P(doc i | cluster k) = sum log P(term m = 0/1 | cluster k) for all term m
        ##            = sum log P(term m = 1 | cluster k) for all term m in doc i
        ##            + sum log P(term m = 0 | cluster k) for all term m not in doc i
        ##
        ## Here we take the second method because there are only 4 terms.
        ##
        sumlogprob1 <- colsums(logcond1[termid[docid == i], ])    # terms in doc i
        sumlogprob0 <- colsums(logcond0[-termid[docid == i], ])   # terms not in doc i
        h[i, ] <- logprior + sumlogprob1 + sumlogprob0

    }
    ## logliks
    loglik[r] <- sum(apply(h, MARGIN=1, FUN=logsumexp))

    ## P(cluster k | doc i) 
    h <- t(apply(h, MARGIN=1, FUN=function (logps)  prop.table(exp(logps - max(logps)))))


    ##
    ## M-Step: update logprior and logcond
    ##
    num_cluster <- colSums(h)                  # expected number of cluster k
    logprior <- log(prop.table(num_cluster))

    for(w in 1:M) {
        ## compute conditional probability of term w given cluster k
        ##
        ## P(term w | cluster k)  = count(term w, cluster k) / count(cluster k)
        ##
        num_cluster_with_term <- colsums(h[docid[termid==w], ])    # expected number of cluster k with term w
        cond1[w, ] <- num_cluster_with_term / num_cluster
    }
    logcond1 <- log(cond1)
    logcond0 <- log(1 - cond1)
    
}

plot(loglik[2:R])
loglik[R]

exp(logprior)
cond1

##
## compare em clustering result with true class label
##
xtabs(~ apply(h, MARGIN=1, FUN=which.max) +
          apply(doc_class_mat, MARGIN=1, FUN=which.max))
