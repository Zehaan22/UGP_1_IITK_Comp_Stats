source("Toy_NUTS.R")


# For Theoretical HMC
chain1 <- normalHMC(s = .1)
chain2 <- normalHMC(s = 1)
chain3 <- normalHMC(s = 5)

# Comparing the chains through plots
plot(density(chain1), col = "blue")
lines(density(chain2), col = "green")
lines(density(chain3), col = "red")

plot.ts(chain1, col = "blue")
plot.ts(chain2, col = "green")
plot.ts(chain3, col = "red")

# For LeapFrog
#Keeping L * epsilon = 1
chain1 <- normalLF_HMC(L = 10, eps = .1)
chain2 <- normalLF_HMC(L = 1, eps = 1)
chain3 <- normalLF_HMC(L = 100, eps = .01)

# Comparing the chains through plots
plot(density(chain1), col = "blue")
lines(density(chain2), col = "green")
lines(density(chain3), col = "red")

plot.ts(chain1, col = "blue")
plot.ts(chain2, col = "green")
plot.ts(chain3, col = "red")

#Trying out different values of epsilon
chain1 <- nut_sampler(e = .1)
chain2 <- nut_sampler(e = 1)
chain3 <- nut_sampler(e = .01)

# Comparing the chains through plots
plot(density(chain1), col = "blue")
lines(density(chain2), col = "green")
lines(density(chain3), col = "red")

plot.ts(chain1, col = "blue")
plot.ts(chain2, col = "green")
plot.ts(chain3, col = "red")