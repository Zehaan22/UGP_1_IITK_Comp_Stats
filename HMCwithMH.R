N <- 1e4
k <- 50
iid.samples <- rt(n = N, df = k)
plot(density(iid.samples))

# Solution : 
t_mh <- function(N = 1e3, k, h){
  
  out.vect <- numeric(length = N)
  acc.prob <- 0
  out.vect[1] <- 2
  
  for(t in 2:N){
    prop.curr <- rnorm(1,out.vect[t-1],sqrt(h))
    y.curr <- prop.curr^3
    
    alpha <- (dt(y.curr, k)*sqrt(out.vect[t-1]))/(dt(out.vect[t-1],k)*sqrt(y.curr)) # Since q is normal hence q(qbrt(x)|z)/q(qbrt(z)|x) = 1    
    
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

chain <- t_mh(N,k,20)
plot(density(iid.samples))
lines(density(chain), col = "red")
plot(density(chain))

plot.ts(iid.samples)
lines(1:N,chain,col = "red")