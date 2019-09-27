
fit.twonorm <- function(y, par, max_iter, var.equal=FALSE, trace=FALSE)
{
  # initialization
  loglik <- c();
  aic <- c(); bic <- c()
  res = c()
  par_trace=c(par)

  for (k in 1:max_iter) {
    # E step
    d1 = par[1]*dnorm(y, par[2], par[3])
    d2 = (1-par[1])*dnorm(y, par[4], par[5])
    h1 = d1/(d1+d2)
    h2 = 1-h1
    #responsibility trace
    res = rbind(res, h2)
    # M step
    par[1] = mean(h1)
    par[2] = sum(h1*y)/sum(h1)
    par[3] = sqrt(sum(h1*(y^2))/sum(h1) - (par[2])^2)
    par[4] = sum(h2*y)/sum(h2)      
    par[5] = sqrt(sum(h2*(y^2))/sum(h2) - (par[4])^2)

    par_trace = rbind(par_trace, par)      
    if(sum(abs(1-(par_trace[k]/par_trace[k+1]))) < 0.0000000001) break   

    like = par[1]*dnorm(y,par[2],par[3]) + (1-par[1])*dnorm(y,par[4],par[5])
    loglik = rbind(loglik, log(like))

    aic = c(aic, -2*sum(log(like)) + 2*length(par)) #Akaike information criterion
    bic = c(bic, -2*sum(log(like)) + length(par)*log(length(y)))
  }
  if (trace){out = list(parameters=c(par), loglik=loglik, params=par_trace, AIC=aic, BIC=bic, iterations = k)}
  else{out <- list(parameters=c(par), iterations= k)}
  out
}

y = c(-0.39, 0.12, 0.94, 1.67, 1.76, 2.44, 3.72, 4.28, 4.92, 5.53,
+        0.06, 0.48, 1.01, 1.68, 1.80, 3.25, 4.12, 4.60, 5.28, 6.22)

initial_values = c(0.5, 1, sqrt(var(y)), 6, sqrt(var(y)))

print(initial_values)

fit.twonorm(y, par=initial_values, max_iter=200, var.equal=FALSE, trace=FALSE)

hist(y)

result = fit.twonorm(y, par=initial_values, max_iter=200, var.equal=FALSE, trace=TRUE)$parameters
print(result)

my_theory = list(p=c(res[1], 1-res[1]), mu = c(res[2],res[4]), sigma=c(res[3], res[5]))

!install.packages("mixtools")

library(mixtools)

normalmixEM(y)

res_R = normalmixEM(y)
list(p = res_R$lambda, mu=res_R$mu, sigma=res_R$sigma)

my_theory

plot(res_R, which=2)
built-in function을 사용했습니다.
my.ecdf = ecdf(y); my.ecdf

plot(my.ecdf)

res = fit.twonorm(y, par=initial_values, max_iter=200, var.equal=FALSE, trace=TRUE)

my.ordered = sort(y); my.ordered

nrow(res$loglik)

plot(x=1:res$iterations, y=rowSums(res$loglik))

plot(x= 1:50, y=rowSums(res$loglik)[1:50])

plot(x= 1:30, y=rowSums(res$loglik)[1:30])

plot(x= 1:30, y=rowSums(res$loglik)[1:30])

plot(x= 1:20, y=rowSums(res$loglik)[1:20])

plot(x= 1:10, y=rowSums(res$loglik)[1:10])
