# Modified Rosenbrock target

d <- 40

dmr <- function(x, b = 1){ # Modifies Rosenbrock density 
  d <- length(x)
  ret <- 1
  for(i in 1:(d/2)){
    s <- 99*(i-1)/(d/2 - 1) + 1
    t1 <- exp(-0.5*(1 / s)*(x[i] - sqrt(2*s) * b)^2) 
    t2 <- exp(-0.5*(x[2*i] - (1/sqrt(2*s)) * (x[2*i-1]^2) / (1 + x[i]^2/(4*s)))^2)
    ret <- ret*t1*t2
  }
  return(ret)
}

drbf <- function(x){
  
}

grad_mr <- function(x){
  d <- length(x) - 1
  ret <- c(400*(x[2] - x[1]^2)*x[1] + 2*(x[1] - 1))
  for(i in 2:d-1){
    t1 <- 400*(x[i+1] - x[i]^2)*x[i] + 2*(x[i] - 1)
    t2 <- -200*(x[i] - x[i-1]^2)
    ret <- c(ret, t1+t2)
  }
  
  ret <- c(ret,-200*(x[d] - x[d-1]^2))
  return(ret)
}

bivariate_dmr <- function(x){
  t1 <- exp(-(100/20)*(x[2] - x[1]^2)^2)
  t2 <- exp(-(1/20)*(1-x[1]^2))
  return(c(t1,t2))
}

## Leap Frog method
LF_sampler_MR <- function(L = 1, eps = .1, d = 40, n = 1e4)
{
  qt <- matrix(nrow = n, ncol = d)
  qt[1,] <- 1 # starting here
  accept <- 1
  # vectorizing this outside to save time
  momentums <- matrix(rnorm(d*n), nrow = n, ncol = d)
  for(k in 2:n)
  {
    # new momentum
    p <- momentums[k,]
    q_prop <- qt[k-1,]
    # one half-Euler step
    p_prop <- p - eps/2 * grad_mr(q_prop)
    for(l in 1:L) # let the frog leap!
    {
      q_prop <- q_prop + eps * p_prop
      if(l != L) p_prop <- p_prop - eps*grad_mr(q_prop)
    }
    # one last half-Euler step
    p_prop <- p_prop - eps/2 * grad_mr(q_prop)
    # Accept-reject
    log.ratio <- (-t(q_prop)%*%q_prop - t(p_prop)%*%p_prop + t(qt[k-1,])%*%qt[k-1,] + t(p)%*%p)/2
    if(log(runif(1)) < as.numeric(log.ratio))
    {
      qt[k,] <- q_prop
      accept <- accept + 1
    } else{
      qt[k,] <- qt[k-1,]
    }
  }
  print(paste("Acceptance = ", accept/n))
  return(qt)
}
library(mcmcse)

chain1 <- LF_sampler_MR(L = 25, eps = .6)
ess(chain1)
chain2 <- LF_sampler_MR(L = 1, eps = 1)
ess(chain2)
chain3 <- LF_sampler_MR(L = 100, eps = .01)
ess(chain3)

heatmap()