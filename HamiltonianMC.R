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

plot(density(chain1), col = "blue")
lines(density(chain2), col = "green")
lines(density(chain3), col = "red")

plot.ts(chain1, col = "blue")
plot.ts(chain2, col = "green")
plot.ts(chain3, col = "red")