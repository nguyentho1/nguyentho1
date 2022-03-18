#Step 1: GIPr(phi,lambda_GIPr)
#Step 2: normal distribution F_N(x|omega,sigma)
#Step 3: Markov chain 
#Step 4: ARL 
################################################
# GGIPr(phi,lambda_GIPr) & other equations in order to construct GIPr(phi,lambda_GIPr) algorithm 

#construct function g0(r, phi), g1(r, phi), g2(r, phi)

#compute g0(r, phi) for p.m.f_GIPr(x|phi, lambda)

####   g0(r, phi)
g0 <- function(x,phi){
  if(phi < 1){phi*(1 - phi^(x+1))/(1 - phi)
  }
  else{x + 1
  }
}

####   g1(x, phi) for E(X) 

g1 <- function(x,phi){
  
  if(phi < 1){(phi^2)*(x*phi^(x+1) - (x+1)*phi^x + 1)/(1 - phi)^2
  }
  else{0.5*x*(x+1)
  }
}

## computing expected value (X) = mu_GIPr from g1
EX_GIPr <- function(x,r,phi,lambda_GIPr){(g1(x,phi) + (r + 1 - g0(x,phi))*lambda_GIPr)/(r + 1)}

####   g2(x, phi) for E(X^2)

g2 <- function(x,phi){
  if(phi < 1){phi*(phi*(phi+1) - phi^(x+1)*((phi*x)^2 - phi*(2*x^2 + 2*x - 1) + (x+1)^2))/(1-phi)^3
  }
  else{x*(2*x^2 + 3*x + 1)/6
  }
}
## computing expected value (X^2) from g2
EX2_GIPr <- function(x,r,phi,lambda_GIPr){(g2(x,phi) + (r + 1 - g0(x,phi))*lambda_GIPr*(lambda_GIPr + 1))/(r+1)}

# compute var(X) = E(X^2) - (E(X))^2
Var_GIPr <- EX2_GIPr - EX_GIPr^2
#############################################

#Step 2: normal distribution F_N(x|omega,sigma), f_N(x|omega,sigma)
#F_N(x|omega,sigma) #pnorm(x, mean, sd)
#we have omega ~ mean; sigma ~ sd 
#f_N(x|omega,sigma) ~ dnorm(x, mean, sd)
#############################################

#Step 3: Compute Var*(X) = Var(X) + sigma^2 
#denote Var*(X) ~ Var(X) "continuousify"
VarX_continuousify = Var_GIPr + sigma^2 

# Compute UCL, LCL 
UCL_continuousify <- EX_GIPr + K*sqrt((lambda_EWMA*VarX_continuousify)/(2 - lambda_EWMA))
### we divide length intervals from 0 to UCL_continuousify into c subintervals with 2*delta width
## c <-    (filling a certain quantity later)
delta <- UCL_continuousify/(2*c) 

#the midpoint of the c_th subinterval is equal to (2*c - 1)/delta 
H_j <- (2*c - 1)/delta 
#step 4: Compute ARL1, SDRL1 based on transition probability matrix of Markov chain 

###ARL_GIPr_EWMA_continuousify   
W <- matrix(0, ncol = p, nrow = p)
p <- c + 1

for (i in 1:p){
  t <- (lambda_continuousify*Hj - Hj)(lambda_continuousify)
  l <- 0 
  while(TRUE){ 
    t0 <- max(0, (lambda_continuosify*x + (1 - lambda_continuousify)*H(t-1)
  
  
  
  
  
  }
















