library(pracma)
N <- 1e4
k <- 50
iid.samples <- rt(n = N, df = k)
plot(density(iid.samples))

# Solution : 
t_mh <- function(N = 1e3, k, h){
  
  out.vect <- numeric(length = N)
  acc.prob <- 0
  out.vect[1] <- 5
  
  for(t in 2:N){
    prop.curr <- rnorm(1,out.vect[t-1],sqrt(h))
    y.curr <- prop.curr^3
    
    alpha.num <- (dnorm(nthroot(out.vect[t-1],3),y.curr,sqrt(h))) * (dt(y.curr, k)) * (nthroot(out.vect[t-1],3)^(-2))
    alpha.den <- (dnorm(prop.curr,out.vect[t-1],sqrt(h))) * (dt(out.vect[t-1], k)) * (prop.curr^(-2))
    
    alpha <- alpha.num / alpha.den
    
    alpha <- min(1,alpha)
    
    U.curr <- runif(1)
    if(U.curr <= alpha){
      out.vect[t] <- prop.curr
      acc.prob <- acc.prob + 1
    }
    else{
      out.vect[t] <- out.vect[t-1]
    }
  }
  
  print(acc.prob/N)
  return(out.vect)
}

chain <- t_mh(N,k,10)
plot(density(iid.samples))
lines(density(chain), col = "red")
plot(density(chain))

plot.ts(iid.samples)
lines(1:N,chain,col = "red")