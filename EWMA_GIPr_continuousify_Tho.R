#Step 1: GIPr(phi,lambda)
#Step 2: normal distribution F_N(x|omega,sigma)
#Step 3: Markov chain 
#Step 4: ARL 
################################################
# GGIPr(phi,lambda) & other equations in order to construct GIPr(phi,lambda) algorithm 

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

####    g2(x, phi) for E(X^2)

g2 <- function(x,phi){
  if(phi < 1){phi*(phi*(phi+1) - phi^(x+1)*((phi*x)^2 - phi*(2*x^2 + 2*x - 1) + (x+1)^2))/(1-phi)^3
  }
  else{x*(2*x^2 + 3*x + 1)/6
  }
}


# compute var(X) = E(X^2) - (E(X))^2








#############################################

#Step 2: normal distribution F_N(x|omega,sigma), f_N(x|omega,sigma)
#F_N(x|omega,sigma) #pnorm(x, mean, sd)
#we have omega ~ mean; sigma ~ sd 
#f_N(x|omega,sigma) ~ dnorm(x, mean, sd)
#############################################

#Step 3: Compute Var*(X) = Var(X) + sigma^2 
# Compute UCL, LCL 

#step 4: Markov chain 

#Step 5: Compute ARL1, SDRL1,...


















