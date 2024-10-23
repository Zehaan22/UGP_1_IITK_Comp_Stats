source("Toy_NUTS.R")
library(mcmcse)

# For Theoretical HMC
print("Theoretical HMC:")
chain1 <- normalHMC(s = .1)
ess(chain1)
chain2 <- normalHMC(s = 1)
ess(chain2)
chain3 <- normalHMC(s = 5)
ess(chain3)
print("###########################")

# Comparing the chains through plots
par(mfrow = c(1,1))
plot(density(chain1), col = "blue", main = "Theoretical HMC")
lines(density(chain2), col = "green")
lines(density(chain3), col = "red")

par(mfrow = c(3,1))
plot.ts(chain1, col = "blue", main = "Theoretical HMC | s = 0.1")
plot.ts(chain2, col = "green", main = "Theoretical HMC | s = 1")
plot.ts(chain3, col = "red", main = "Theoretical HMC | s = 5")

acf(chain1, main = "Theoretical HMC | s = 0.1")
acf(chain2, main = "Theoretical HMC | s = 1")
acf(chain3, main = "Theoretical HMC | s = 5")

# For LeapFrog
#Keeping L * epsilon = 1
print("Leap Frog:")
chain1 <- normalLF_HMC(L = 10, eps = .1)
ess(chain1)
chain2 <- normalLF_HMC(L = 1, eps = 1)
ess(chain2)
chain3 <- normalLF_HMC(L = 100, eps = .01)
ess(chain3)
print("###########################")

# Comparing the chains through plots
par(mfrow = c(1,1))
plot(density(chain1), col = "blue", main = "Leap Frog")
lines(density(chain2), col = "green")
lines(density(chain3), col = "red")

par(mfrow = c(3,1))
plot.ts(chain1, col = "blue", main = "Leap Frog | L = 10 e = 0.1")
plot.ts(chain2, col = "green", main = "Leap Frog | L = 1, e = 1")
plot.ts(chain3, col = "red", main = "Leap Frog | L = 100, e = .01")

acf(chain1, main = "Leap Frog | L = 10 e = 0.1")
acf(chain2, main = "Leap Frog | L = 1, e = 1")
acf(chain3, main = "Leap Frog | L = 100, e = .01")


#Trying out different values of epsilon in NUTS
print("NUTS Algorithm:")
chain1 <- nut_sampler(e = .1)
ess(chain1)
chain2 <- nut_sampler(e = 1)
ess(chain2)
chain3 <- nut_sampler(e = .01)
ess(chain3)
print("###########################")

# Comparing the chains through plots
par(mfrow = c(1,1))
plot(density(chain1), col = "blue", main = "NUTS Algorithm")
lines(density(chain2), col = "green")
lines(density(chain3), col = "red")

par(mfrow = c(3,1))
plot.ts(chain1, col = "blue", main = "NUTS Algorithm | e = 0.1")
plot.ts(chain2, col = "green", main = "NUTS Algorithm | e = 1")
plot.ts(chain3, col = "red", main = "NUTS Algorithm | e = 0.01")

acf(chain1, main = "NUTS Algorithm | e = 0.1")
acf(chain2, main = "NUTS Algorithm | e = 1")
acf(chain3, main = "NUTS Algorithm | e = 0.01")