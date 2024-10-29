## Standard Theoretical HMC
normalHMC <- function(s = 1, n = 1e4)
{
  qt <- numeric(length = n)
  qt[1] <- 0 # starting here
  # don't need a starting value of p
  for(k in 2:n)
  {
    # new momentum
    p <- rnorm(1)
    #initial conditions for new H
    r2 <- qt[k-1]**2 + p**2
    # choosing +- with probability 1/2
    r <- sample(c(sqrt(r2), -sqrt(r2)), size = 1)
    a <- acos(qt[k-1]/r)
    # simulating Hamiltonian forward s time units
    qt[k] <- r*cos(a + s)
    p <- -r * sin(a + s)
  }
  return(qt)
}

chain1 <- normalHMC(s = .1)
chain2 <- normalHMC(s = 1)
chain3 <- normalHMC(s = 5)

# Comparing the chains through plots
plot(density(chain1), col = "blue", main = "Density from HMC sampling")
lines(density(chain2), col = "green")
lines(density(chain3), col = "red")
lines(density(rnorm(1e4)), col = "black")
legend("right", c(
  "Truth",
  "S = .1",
  "S = 1",
  "S = 5"
),
fill = c(
  "black",
  "blue",
  "green",
  "red"
))

plot.ts(chain1, col = "blue", ylab = "", main = "Time Series Plot for the HMC Samples")
lines(chain2, col = rgb(0,1, 0, alpha = 0.6))
lines(chain3, col = rgb(1,0, 0, alpha = 0.4))

## Leap Frog method

normalLF_HMC <- function(L = 10, eps = .1, n = 1e4)
{
  qt <- numeric(length = n)
  qt[1] <- 0 # starting here
  accept <- 1
  # vectorizing this outside to save time
  momentums <- rnorm(n)
  for(k in 2:n)
  {
    # new momentum
    p <- momentums[k]
    q_prop <- qt[k-1]
    # one half-Euler step
    p_prop <- p - eps/2 * q_prop
    for(l in 1:L) # let the frog leap!
    {
      q_prop <- q_prop + eps * p_prop
      if(l != L) p_prop <- p_prop - eps*q_prop
    }
    # one last half-Euler step
    p_prop <- p_prop - eps/2 * q_prop
    # Accept-reject
    log.ratio <- (-q_prop**2 - p_prop**2 + qt[k-1]**2 + p**2)/2
    if(log(runif(1)) < log.ratio)
    {
      qt[k] <- q_prop
      accept <- accept + 1
    } else{
      qt[k] <- qt[k-1]
    }
  }
  print(paste("Acceptance = ", accept/n))
  return(qt)
}

#Keeping L * epsilon = 1
chain1 <- normalLF_HMC(L = 10, eps = .1)
chain2 <- normalLF_HMC(L = 1, eps = 1)
chain3 <- normalLF_HMC(L = 100, eps = .01)

plot(density(chain1), col = "blue", main = "Density from HMC sampling")
lines(density(chain2), col = "green")
lines(density(chain3), col = "red")
lines(density(rnorm(1e4)), col = "black")
legend("right", c(
  "Truth",
  "L = 10, eps = .1",
  "L = 1, eps = 1",
  "L = 100, eps = .01"
),
fill = c(
  "black",
  "blue",
  "green",
  "red"
))

plot.ts(chain1, col = "blue", ylab = "", main = "Time Series Plot for the HMC Samples")
lines(chain2, col = rgb(0,1, 0, alpha = 0.6))
lines(chain3, col = rgb(1,0, 0, alpha = 0.4))
