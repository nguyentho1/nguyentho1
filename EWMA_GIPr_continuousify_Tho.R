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

############################################
# compute var(X) = E(X^2) - (E(X))^2
#Var_GIPr <- EX2_GIPr - EX_GIPr^2
Var_GIPr <- EX2_GIPr(x = x,r = r,phi = phi, lambda_GIPr = lambda_GIPr) - EX_GIPr(x = x,r = r,phi = phi, lambda_GIPr = lambda_GIPr)^2

#############################################

#Step 2: normal distribution F_N(x|omega,sigma), f_N(x|omega,sigma)
#F_N(x|omega,sigma) ~ pnorm(x, mean, sd)
#we have omega ~ mean; sigma ~ sd 
#f_N(x|omega,sigma) ~ dnorm(x, mean, sd)
#############################################

#Step 3: Compute Var*(X) = Var(X) + sigma^2 
#denote Var*(X) ~ Var(X) "continuousify"
VarX_star = Var_GIPr + sigma^2 

# Compute UCL
UCL_star <- EX_GIPr(x = x,r = r,phi = phi, lambda_GIPr = lambda_GIPr) + K*sqrt((lambda_EWMA*VarX_star)/(2 - lambda_EWMA))
### we divide length intervals from 0 to UCL_star into m subintervals with 2*Delta width
Delta <- UCL_star/(2*m) 

#the midpoint of the k_th subinterval is equal to (2*k - 1)/delta 
H_k <- (2*k - 1)/Delta 
#######################################################################################
#step 4: Compute ARL1, SDRL1 based on transition probability matrix of Markov chain 
#p.m.f GIPr (omega,theta) in which omega = x, thetaX = c(phi, lambda_GIPr)

ff_GIPr<-function(omega,theta){
  ifelse(x>=(r+1),(((r+1-g0(r,phi))*dpois(x,lambda_GIPr,log=FALSE)))/(r+1),
         ((phi^(x+1)+(r+1-g0(r,phi))*dpois(x,lambda_GIPr,log=FALSE)))/(r+1))
}

F_Xstar <- sum(omega)*ff_GIPr(omega = x,theta = thetaX)*pnorm(x,omega,sigma)
###########################################################################
#Denote:
  r = 2
  phi = 0.5 
  lambda_GIPr = 5
  x = 3
  H = 10.8  #based on sample size of the previous paper
  m = 100*H
  tn = m + 1
  K = 3
  lambda_EWMA = 0.2
  thetaT = c(phi, lambda_GIPr)
sigma = 0.15
# Creat matrix Q and transition probability matrix 
tn = m + 1 
Q <- matrix(0, nrow = tn, ncol = tn)
# Creat Qk0 and Qkj generic elements 
F_Xstar_k0 <- function(x,theta, sigma){sum(omega)*(ff_GIPr(omega,theta)*pnorm(x,omega,sigma))}

F_Xstar_kj <- function(x1,x2,theta,sigma){
  sum(omega)*(ff_GIPr(omega,theta)*pnorm(x1,omega,sigma) 
              - ff_GIPr(omega,theta) * pnorm(x2,omega,sigma))}
##################################################
for (k in 1:tn){Hk = (2*k - 1)*Delta
  for( j in 1:tn){
  if( j == 1 ){
    x0 <- (-1 + lambda_EWMA)*Hk/lambda_EWMA + 1
    Q[k,j] <- F_Xstar_k0(x = x0, theta = thetaT, sigma = sigma)
  } 
  else{
  x1 <- (2*j*Delta + Hk*Delta - Hk)/lambda_EWMA + 1;
  x2 <- (2*j*Delta - 2*Delta + Hk*Delta - Hk)/lambda_EWMA + 1;
  }
  F_Xstar_kj <- function(x1,x2,theta,sigma){
    sum(omega)*(ff_GIPr(omega,theta)*pnorm(x1,omega,sigma) 
                - ff_GIPr(omega,theta)*pnorm(x2,omega,sigma))}
  Q[k,j] = Q[k,j] + F_Xstar_kj(x1,x2,thetaT,sigma);
  }
 }
#################################################
  I <- diag(1,tn)
  One <- array(1,tn)
  ARL <- solve(I-Q, One)
  print(ARL[1])
  nu1 <- ARL[1]

#######################################
#SDRL1

W<-solve(I-Q)
avec<-rep(0,tn)
y1<-rep(1,tn)
nu2<-2*avec%*%W%*%W%*%(Q)%*%y1
SDRL<-sqrt(nu2-nu1^2)
print(c(r,phi,omega,lambda_EWMA,nu1,SDRL))
#########################################################################
  

















