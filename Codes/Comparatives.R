setwd("~/Research/UGP_1_IITK_Comp_Stats/Codes")
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

vals <- numeric(100)
for(i in 1:100){
  chain <- normalHMC(s = i)
  vals[i] <- ess(chain)
}
par(mfrow = c(1,1))
plot.ts(vals/1e4, main = "Cyclic nature of s", ylab = "ess/N")

# Comparing the chains through plots
par(mfrow = c(1,1))
plot(density(chain1), col = "blue", main = "Theoretical HMC")
lines(density(chain2), col = "green")
lines(density(chain3), col = "red")
lines(density(rnorm(1e4)), col = "black")
legend("right", c(
  "Truth",
  "s = .1",
  "s = 1",
  "s = 5"
),
fill = c(
  "black",
  "blue",
  "green",
  "red"
))



par(mfrow = c(3,1))
plot.ts(chain1, col = "blue", main = "Theoretical HMC | s = 0.1", cex.main = 2)
plot.ts(chain2, col = "green", main = "Theoretical HMC | s = 1", cex.main = 2)
plot.ts(chain3, col = "red", main = "Theoretical HMC | s = 5", cex.main = 2)

par(cex.main = 2)
acf(chain1, main = "Theoretical HMC | s = 0.1", cex.main = 3)
acf(chain2, main = "Theoretical HMC | s = 1",  cex.main = 3)
acf(chain3, main = "Theoretical HMC | s = 5", cex.main = 3)

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

ess.heat.map.plt <- matrix(0,nrow = 50,ncol = 50)
for(l in 1:50){
  for(eps in 1:50){
    chain <- normalLF_HMC(L = l, eps = eps/10)
    ess.heat.map.plt[l,eps] <- ess(chain)
  }
  print(l)
}

library(ggplot2)

# Dummy data
L <- 1:50
eps <- (1:50)/10
data <- expand.grid(L=L, eps=eps)
data$Ess <- c(ess.heat.map.plt)/1e4

# Heatmap 
ggplot(data, aes(L, eps, fill= Ess)) + 
  geom_tile()


# Create example data        
data <- ess.heat.map.plt

# Column names    
colnames(data) <- (1:50)/10
rownames(data) <- 1:50

# Remove dendrogram
# Manual color range
my_colors <- colorRampPalette(c("yellow", "blue"))

# Heatmap with margins around the plot
heatmap(data, col = my_colors(100), main = "Heatmap for ESS/N", 
        xlab = "Eps", ylab = "L",Colv = NA, Rowv = NA)

# Comparing the chains through plots
par(mfrow = c(1,1))
plot(density(chain1), col = "blue", main = "Leap Frog")
lines(density(chain2), col = "green")
lines(density(chain3), col = "red")
lines(density(rnorm(1e4)), col = "black")
legend("right", c(
  "Truth",
  "eps = .1",
  "eps = 1",
  "eps = .01"
),
fill = c(
  "black",
  "blue",
  "green",
  "red"
))

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
lines(density(rnorm(1e4)), col = "black")
legend("right", c(
  "Truth",
  "eps = .1",
  "eps = 1",
  "eps = .01"
),
fill = c(
  "black",
  "blue",
  "green",
  "red"
))

par(mfrow = c(3,1))
plot.ts(chain1, col = "blue", main = "NUTS Algorithm | e = 0.1")
plot.ts(chain2, col = "green", main = "NUTS Algorithm | e = 1")
plot.ts(chain3, col = "red", main = "NUTS Algorithm | e = 0.01")

acf(chain1, main = "NUTS Algorithm | e = 0.1")
acf(chain2, main = "NUTS Algorithm | e = 1")
acf(chain3, main = "NUTS Algorithm | e = 0.01")
