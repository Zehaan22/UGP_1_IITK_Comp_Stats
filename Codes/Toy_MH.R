# Implementing a toy problem using Metropolois-Hastings Algorithm
set.seed(42)
##############################
# Problem : We want to run the MH algorithm to sample from a t dist.

N <- 1e4
k <- 73
iid.samples <- rt(n = N, df = k)
plot(density(iid.samples))

# Density function is dt

##############################


##############################
# Solution : 
t_mh <- function(N = 1e3, k, h){
  
  out.vect <- numeric(length = N)
  acc.prob <- 0
  out.vect[1] <- 30
  
  for(t in 2:N){
    prop.curr <- rnorm(1,out.vect[t-1],sqrt(h))
    
    alpha <- dt(prop.curr, k)/dt(out.vect[t-1],k)
    
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

plot.ts(iid.samples)
lines(1:N,chain,col = "red")
##############################