# Trying to sample from a multivariate norm

grad_skw_norm <- function(x, s.2 = 1, a = 1){
  d <- length(x)
  ret.grad <- numeric(d)
  dnorms <- dnorm(a*x/sqrt(s.2))
  pnorms <- pnorm(a*x/sqrt(s.2))
  ret.grad <- x/(s.2) - dnorms/pnorms
  return(ret.grad)
}


LF_smplr_skw_norm <- function(L = 10, eps = .01, d = 40, n = 1e4)
{
  qt <- matrix(nrow = n, ncol = d)
  qt[1,] <- 0
  accept <- 1
  momentums <- matrix(rnorm(d*n), nrow = n, ncol = d)
  for(k in 2:n)
  {
    # new momentum
    p <- momentums[k,]
    q_prop <- qt[k-1,]
    # one half-Euler step
    p_prop <- p - eps/2 * grad_skw_norm(q_prop)
    for(l in 1:L) 
    {
      q_prop <- q_prop + eps * p_prop
      if(l != L) p_prop <- p_prop - eps*grad_skw_norm(q_prop)
    }
    # one last half-Euler step
    p_prop <- p_prop - eps/2 * grad_skw_norm(q_prop)
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
  return(list(qt,accept/n))
}

# Just to check if everything works
chain <- LF_smplr_skw_norm()

# Running the heatmap code

grid.length.eps <- 40
grid.length.Ls <- 15
n <- 1e5

EPSs <- seq(0.01,0.25, length = grid.length.eps)
Ls <- 1:grid.length.Ls
ess.heat.map.plt <- matrix(0, nrow = grid.length.Ls, ncol = grid.length.eps)
ess.heat.map.plt.avg <- matrix(0, nrow = grid.length.Ls, ncol = grid.length.eps)

for(l in 1:grid.length.Ls){
  for(eps in 1:grid.length.eps){
    chain <- LF_smplr_skw_norm(L = l, eps = EPSs[eps], n = n)
    if(chain[[2]] <= 2e-4){
      ess.heat.map.plt[l,eps] <- 0
      ess.heat.map.plt.avg[l,eps] <- 0
      break
    }else {
      foo <- ess(chain[[1]])
      ess.heat.map.plt[l,eps] <- min(foo)
      ess.heat.map.plt.avg[l,eps] <- mean(foo)
    }
  }
  print(l)
}

write.csv(ess.heat.map.plt, file = "Min_Heat_map_data_skwd.csv")
write.csv(ess.heat.map.plt.avg, file = "Avg_Heat_map_data_skwd.csv")

library(ggplot2)

foo <- read.csv("Min_Heat_map_data_MVT_norm.csv")
foo <- foo[,2:41]
L <- Ls
eps <- EPSs
data <- expand.grid(L=L, eps=eps)
data$Ess <- c(ess.heat.map.plt)/n

# Heatmap 
ggplot(data, aes(L, eps, fill= Ess)) + 
  geom_tile()+
  guides(fill=guide_legend(title="min(ESS) / N"))


