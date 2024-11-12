## Code for the complete nuts algorithm

# Defining the auxiliary functions
leapFrog <- function(p,q,e){
  p_ <- p - (e/2)*q
  q_ <- q + e*p_
  p_ <- p_ - (e/2)*q_
  return(list(p_,q_))
}

travelRoad <- function(p,q,u,v,j,e, del_max = 1e-6){
  if (j == 0){ #Base Case
    rets <- leapFrog(p,q,v*e)
    p_ <- rets[[1]]
    q_ <- rets[[2]]
    if(u <= exp(-(q^2+p^2)/2)){
      set_C <- c(rets[[1]],rets[[2]])
    }
    else{
      set_C <- NaN
    }
    s_ <- -(q_^2)/2 - (p^2)/2 > log(u) - del_max
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
    s_ <- s_1*s_2*((q_pls-q_min)*p_min>=0)*((q_pls-q_min)*p_pls>=0)
    Set_C_1 <- rbind(Set_C_1, Set_C_2)
    return(list(p_min,q_min,p_pls,q_pls,Set_C_1,s_))
  }
}

### New addition
find_e <- function(q){
  e <- 1
  p <- rnorm(1)
  new_st <- leapFrog(p,q,e)
  st_ratio <- (dnorm(new_st[[1]])*dnorm(new_st[[2]]))/(dnorm(p)*dnorm(q))
  if (st_ratio > 0.5){
    a <- 1
  }else{
    a <- -1
  }
  while(st_ratio^a > 2^(-a)){
    e <- (2^a)*e
    new_st <- leapFrog(new_st[[1]],new_st[[2]],e)
    st_ratio <- (dnorm(new_st[[1]])*dnorm(new_st[[2]]))/(dnorm(p)*dnorm(q))
  }
  return(e)
}

nut_sampler <- function(e,n=1e4){
  q = numeric(n)
  jumps = numeric(n)
  q_0 = 0
  for(i in 1:n){
    
    p <- rnorm(1)
    u <- runif(1,0,exp(-((p^2)+(q_0^2))/2))
    q_min <- q_0
    q_pls <- q_0
    p_min <- p
    p_pls <- p
    j = 0
    
    Set_C <- matrix(c(q_0,p),ncol = 2)
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
        Set_C <- rbind(Set_C,Set_C_1)
      }
      s <- s_1*((q_pls-q_min)*p_min>=0)*((q_pls-q_min)*p_pls>=0)
      j <- j + 1
    }
    n_0 <- sample(1:dim(Set_C)[1],1)
    q_0 <- Set_C[n_0,]
    q_0 <- q_0[2]
    q[i] <- q_0
    jumps[i] <- j
  }
  print(summary(jumps))
  return(q)
}

## Running an example

nut_sampler(find_e(0))
