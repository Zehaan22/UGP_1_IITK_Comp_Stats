# Trying to sample from a multivariate norm
grad_mvt_norm <- function(x, s.2 = 1, a=3){
  d <- length(x)
  ret.grad <- numeric(d)
  dnorms <- dnorm(a*x/sqrt(s.2))
  pnorms <- pnorm(a*x/sqrt(s.2))
  ret.grad <- x/(s.2) - dnorms/pnorms
  return(ret.grad)
}


# Defining the auxiliary functions
leapFrog <- function(p,q,e){
  p_ <- p - (e/2)*grad_mvt_norm(q)
  q_ <- q + e*p_
  p_ <- p_ - (e/2)*grad_mvt_norm(q_)
  return(list(p_,q_))
}

travelRoad <- function(p,q,u,v,j,e, d = 40, del_max = 1e-6){
  if (j == 0){ #Base Case
    rets <- leapFrog(p,q,v*e)
    p_ <- rets[[1]]
    q_ <- rets[[2]]
    if(all(u <= exp(-(q^2+p^2)/2))){
      set_C <- array(0,dim = c(1,2,d))
      set_C[1,1,] <- rets[[1]]
      set_C[1,2,] <- rets[[2]]
    }
    else{
      set_C <- array(0,dim = c(1,2,d))
    }
    s_ <- all(-(q_^2)/2 - (p^2)/2 > log(u) - del_max)
    return(list(p_,q_,p_,q_,set_C,s_))
  }
  else{ # Recursive Case
    rets <- travelRoad(p,q,u,v,j-1,e)
    p_min <- rets[[1]]
    q_min <- rets[[2]]
    p_pls <- rets[[3]]
    q_pls <- rets[[4]]
    Set_C_1 <- rets[[5]]
    s_1 <- rets[[6]]
    if (v == -1){
      rets_min <- travelRoad(p_min,q_min,u,v,j-1,e)
      p_min <- rets_min[[1]]
      q_min <- rets_min[[2]]
      Set_C_2 <- rets_min[[5]]
      s_2 <- rets_min[[6]]
    }
    else{
      rets_min <- travelRoad(p_pls,q_pls,u,v,j-1,e)
      p_pls <- rets_min[[3]]
      q_pls <- rets_min[[4]]
      Set_C_2 <- rets_min[[5]]
      s_2 <- rets_min[[6]]
    }
    s_ <- s_1*s_2*all((q_pls-q_min)*p_min>=0)*all((q_pls-q_min)*p_pls>=0)
    Set_C_1 <- abind(Set_C_1, Set_C_2,along = 1)
    return(list(p_min,q_min,p_pls,q_pls,Set_C_1,s_))
  }
}

nut_sampler <- function(e = 0.1, d=40, n=1e4){
  q = matrix(0,nrow = n, ncol = d)
  momentums <- matrix(rnorm(d*n),nrow = n, ncol = d)
  jumps = numeric(n)
  q_0 = q[1,]
  for(i in 1:n){
    p <- momentums[i,]
    maxes <- exp(-((p^2)+(q_0^2))/2)
    u <- numeric(d)
    for(k in 1:d){
      u[k] <- runif(1,0,maxes[k])
    }
    
    q_min <- q_0
    q_pls <- q_0
    p_min <- p
    p_pls <- p
    j = 0
    
    Set_C <- array(0,dim = c(1,2,d))
    Set_C[1,1,] <- q_0
    Set_C[1,2,] <- p
    s = 1
    while(s == 1){
      v = sample(c(-1,1),1)
      if(v == -1){
        rets <- travelRoad(p_min,q_min,u,v,j,e)
        p_min <- rets[[1]]
        q_min <- rets[[2]]
        Set_C_1 <- rets[[5]]
        s_1 <- rets[[6]]
      }
      else{
        rets <- travelRoad(p_pls,q_pls,u,v,j,e)
        p_pls <- rets[[3]]
        q_pls <- rets[[4]]
        Set_C_1 <- rets[[5]]
        s_1 <- rets[[6]]
      }
      if(s_1 == 1){
        Set_C <- abind(Set_C,Set_C_1,along = 1)
      }
      s <- s_1*(all((q_pls-q_min)*p_min>=0))*(all((q_pls-q_min)*p_pls>=0))
      j <- j + 1
    }
    n_0 <- sample(1:dim(Set_C)[1],1)
    q_0 <- Set_C[n_0,1,]
    # q_0 <- q_0[2]
    q[i,] <- q_0
    jumps[i] <- j
  }
  print(summary(jumps))
  return(list(q,jumps))
}

# just to check if all works
chain <- nut_sampler()
grid.length.eps <- 40
EPSs <- seq(0.01,0.25, length = grid.length.eps)

means <- numeric(grid.length.eps)
maxes <- numeric(grid.length.eps)
mins <-  numeric(grid.length.eps)
#ess.mins <- numeric(grid.length.eps)
for(e in 1:grid.length.eps){
  chain <- nut_sampler(e = EPSs[e])
  means[e] <- mean(chain[[2]])
  maxes[e] <- max(chain[[2]])
  mins[e] <- min(chain[[2]])
  #ess.mins[e] <- min(ess(chain[[1]]))
}

dat <- read.csv("Data/Min_Heat_map_data_skwd.csv")
dat <- dat[1:15,2:41]
fin.ans <- matrix(0,nrow = grid.length.eps, ncol = 5)
for(i in 1:grid.length.eps){
  fin.ans[i,1] <- EPSs[i]
  fin.ans[i,2] <- which.max(dat[,i])
  fin.ans[i,3] <- mins[i]
  fin.ans[i,4] <- means[i]
  fin.ans[i,5] <- maxes[i]
}
write.csv(fin.ans, "Data/Fin_skwd.csv")