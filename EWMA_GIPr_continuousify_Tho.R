#Step 1: GIPr(r,phi,lambda_GIPr)
#Step 2: normal distribution F_N(x|omega,sigma): I use function pnorm(q,mean,sd) for 
#Step 3: Markov chain: I compute UCL, denote m subintervals and define transient states in Q(m + 1, m + 1) matrix with loop 
#Step 4: ARL and SDRL 
######################################################################################################################
# GGIPr(phi,lambda_GIPr) & other equations in order to construct GIPr(r,phi,lambda_GIPr) algorithm                   #                                                                                                                 #
#Create function g0(r, phi), g1(r, phi), g2(r, phi): Like in your previous paper 02 parameters GIPr                  #
#Compute E(X), Var(X)                                                                                                #
######################################################################################################################
# Denote 
# GIPr
r = 2
phi = 0.5 
lambda_GIPr = 5.0

# EWMA + normal distribution 
m = 200                              # m is the number of subgroups [0,UCL]
n = m + 1                            # 0 -> m => Q matrix
lambda_EWMA = 0.2 
K = 3
omega = 2.5                         # Do I need to create a data.frame to calculate omega or sigma?
sigma = 0.3 

#################################### step 1: GIPr(phi,lambda_GIPr)
###  g0(r, phi)
g0 <- function(x,phi){
  if(phi < 1){phi*(1 - phi^(x+1))/(1 - phi)
  }
  else{x + 1
  }
}
###  g1(x, phi) for E(X)
g1 <- function(x,phi){
  if(phi < 1){(phi^2)*(x*phi^(x+1) - (x+1)*phi^x + 1)/(1 - phi)^2
  }
  else{0.5*x*(x+1)
  }
}
### E(X) 
EX_GIPr <- function(x,r,phi,lambda_GIPr){(g1(x,phi) + (r + 1 - g0(x,phi))*lambda_GIPr)/(r + 1)}
###  g2(x, phi) for E(X^2)
g2 <- function(x,phi){
  if(phi < 1){phi*(phi*(phi+1) - phi^(x+1)*((phi*x)^2 - phi*(2*x^2 + 2*x - 1) + (x+1)^2))/(1-phi)^3
  }
  else{x*(2*x^2 + 3*x + 1)/6
  }
}
### (X^2) 
EX2_GIPr <- function(x,r,phi,lambda_GIPr){(g2(x,phi) + (r + 1 - g0(x,phi))*lambda_GIPr*(lambda_GIPr + 1))/(r+1)}

###   Var(X) = E(X^2) - (E(X))^2 
#Var_GIPr <- EX2_GIPr - EX_GIPr^2
Var_GIPr <- EX2_GIPr(x,r,phi, lambda_GIPr) - EX_GIPr(x,r,phi, lambda_GIPr)^2

##############################   p.m.f GIPr ####################################
ff_GIPr<-function(x,r,phi,lambda_GIPr){
  ifelse(x>=(r+1),(((r+1-g0(r,phi))*dpois(x,lambda_GIPr,log=FALSE)))/(r+1),
         ((phi^(x+1)+(r+1-g0(r,phi))*dpois(x,lambda_GIPr,log=FALSE)))/(r+1))
}
##################################  Step 2: # Mixture of normal distribution
###  normal distribution F_N(x|omega,sigma), f_N(x|omega,sigma)

###  from equation (2) mixture of normal distribution 
Fx_cdf <- sum(x)*ff_GIPr(x,r,phi, lambda_GIPr)*pnorm(x,omega,sigma)    # Do I understand and code correctly this function?

### EWMA    
#Compute Var*(X) = Var(X) + sigma^2                         # (equation 4)
Var_EWMA <- Var_GIPr + sigma^2 

###  UCL   
UCL <- EX_GIPr(x,r,phi,lambda_GIPr) + K*sqrt((lambda_EWMA*Var_EWMA)/(2 - lambda_EWMA))

###    Zi*  
jscore <- function(x,k,lambda_EWMA){round((1 - lambda_EWMA)*k + lambda_EWMA*x)}             #(equation 5)

### we divide length intervals from 0 to UCL_star into m subintervals with 2*Delta width
Delta <- UCL/(2*m)                                                  # m is given 200 

###  the midpoint of the k_th subinterval is equal to (2*k - 1)/Delta 
Hk <- (2*k - 1)/Delta 
###########################################################################################################################
# Define transient probabilities matrix 
Q <- matrix(0, nrow = n, ncol = n)
for (k in 1: n){Hk = (2*k - 1)*Delta; 
for( j in 1: n){Hj = (2*j - 1)*Delta
if(j == 1 ){j0 <- (-1 + lambda_EWMA)*Hk/lambda_EWMA
  FX.k0 <- sum(j0)*ff_GIPr(j0,r,theta)*pnorm(j0, omega, sigma) 
  Q[k,j] <- FX.k0
} 
else{j <- max(0, jscore(x,k,lambda_EWMA))
j1 <- (2*j*Delta + Hk*Delta - Hk)/lambda_EWMA
j2 <- (2*j*Delta - 2*Delta + Hk*Delta - Hk)/lambda_EWMA
F1 <- sum(j1)*ff_GIPr(j1,r,phi,lambda_GIPr)*pnorm(j1,omega,sigma)
F2 <- sum(j2)*ff_GIPr(j2,r,phi,lambda_GIPr)*pnorm(j2,omega,sigma)
FX.kj <- F1 - F2 
Q[k,j] = FX.kj
}
}
}
#########################################################################################################################################
### Compute ARL, SDRL
  I <- diag(1,n)
  One <- array(1,n)
  ARL <- solve(I-Q, One)
  print(ARL[1])
  v1 <- ARL[1]

###   SDRL1

W<-solve(I-Q)
avec<-rep(0,n)
y1<-rep(1,n)
v2<-2*avec%*%W%*%W%*%(Q)%*%y1
SDRL<-sqrt(v2-v1^2 + v1)
print(c(r,phi,lambda_GIPr,lambda_EWMA,sigma,nu1,SDRL))
#########################################################################
  

















