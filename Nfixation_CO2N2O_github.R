#install.packages("deSolve")
library(deSolve)


#######################################################################################
#Model parameters:

I_low  <- 0.5
I_boreal_Dentener<-2.926198
I_boreal_Vet<-6.041241
I_boreal_S4<-3.838037
I_temperate_Dentener<-7.970901
I_temperate_Vet<-4.433534
I_temperate_S4<-10.92945
I_tropical_Dentener<-6.217917
I_tropical_Vet<-4.062762
I_tropical_S4<-8.953253

v0 <- 0.2
vF <- 0.2
wF<-107.9255
w0<-215.851
u0<-0.04295999
uF<-0.04295999
psi <- 1/121
m <- 1/5
k <- 3.472
phi <- 0.001

beta0_boreal <- 4.056749
betaF_boreal <- 3.651074
gamma0_boreal <- 0.001253138
gammaF_boreal <- 0.003132845
beta0_temperate <- 4.12859
betaF_temperate <- 4.545
gamma0_temperate <- 0.0006564928
gammaF_temperate <- 0.00235
beta0_tropical <- 5.05
betaF_tropical <- 5.555
gamma0_tropical <- 0.00094
gammaF_tropical <- 0.00235

EF_boreal <- 0.04031742
EF_temperate <- 0.04031742
EF_tropical <- 1.528

Fmax <- 0.01
Fmin_obligate<-Fmax
z_obligate<-1
Fmin_facultative<-0
z_facultative<-1
Fmin_down<-Fmax/2
z_down<-1

#######################################################################################
#Initial conditions: 

Binit <- 1
Einit <- 0

#######################################################################################
#Time parameters

t0 <- 0
dt <- 0.1
tf <- 100
tv <- seq(t0,tf,dt)
tsteps <- length(tv)

#######################################################################################
#Model functions:

N_model_ode <- function(t,x,parms){
  with(as.list(c(parms,x)),{
    BF <- x[1]
    B0 <- x[2]
    L <- x[3]
    A <- x[4]
    E <- x[5]
    
    F <- max(Fmin,min(betaF/(wF*(1+gammaF*(BF+B0)))-z*vF*A,Fmax))
    gF <- min(wF*(vF*A + F),betaF/(1+gammaF*(BF+B0)))
    g0 <- min(w0*v0*A,beta0/(1+gamma0*(BF+B0)))
    
    dBFdt <- BF*(gF-uF)
    dB0dt <- B0*(g0-u0)
    dLdt <- uF*BF/wF + u0*B0/w0 - L*(m + phi)
    dAdt <- I - BF*(gF-wF*F)/wF - B0*g0/w0 - k*A +  m*L - EF*A
    dEdt <- EF*A - psi*E
    
    list(c(dBFdt,dB0dt,dLdt,dAdt,dEdt))
  })
}

N_model_ode_nonfixer <- function(t,x,parms){
  with(as.list(c(parms,x)),{
    B0 <- x[1]
    L <- x[2]
    A <- x[3]
    E <- x[4]
    
    g0 <- min(w0*v0*A,beta0/(1+gamma0*(B0)))
    
    dB0dt <- B0*(g0-u0)
    dLdt <- u0*B0/w0 - L*(m + phi)
    dAdt <- I - B0*g0/w0 - k*A +  m*L - EF*A
    dEdt <- EF*A - psi*E
    
    list(c(dB0dt,dLdt,dAdt,dEdt))
  })
}

#######################################################################################
#Ecosystems with N-fixing trees:

#Assign values to I_var, z_var, Fmin_var, beta0_var, betaF_var, gamma0_var, gammaF_var and EF_var:
I_var<-I_tropical_Dentener
z_var<-z_facultative
Fmin_var<-Fmin_facultative
beta0_var<-beta0_tropical
betaF_var<-betaF_tropical
gamma0_var<-gamma0_tropical
gammaF_var<-gammaF_tropical
EF_var<-EF_tropical

parms_Nfixer <- list(I=I_var,z=z_var,Fmin=Fmin_var,beta0=beta0_var,betaF=betaF_var,
                     gamma0=gamma0_var,gammaF=gammaF_var,EF=EF_var)
Ainit <- (I_var*w0*gamma0_var*m+I_var*gamma0_var*w0*phi-phi*beta0_var+phi*u0)/
  (w0*gamma0_var*(k+EF_var)*(m+phi))*0.5
Linit <- (beta0_var-u0)/(w0*gamma0_var*(m+phi))*0.75

op_Nfixer <- lsoda(c(Binit,Binit,Linit,Ainit,Einit),tv,N_model_ode,parms_Nfixer)
  
BF_Nfixer <- op_Nfixer[,2]
B0_Nfixer <- op_Nfixer[,3]
L_Nfixer <- op_Nfixer[,4]
A_Nfixer <- op_Nfixer[,5]
E_Nfixer <- op_Nfixer[,6]
  
CO2_Nfixer <- -(BF_Nfixer[length(tv)]+B0_Nfixer[length(tv)]-Binit*2)/10000/tf*44/12
N2O_Nfixer <- (E_Nfixer[length(tv)]-Einit)/10000/tf*44/28*298
netCO2N2O_Nfixer <- CO2_Nfixer+N2O_Nfixer

#######################################################################################
#Ecosystems without N-fixing trees:

#Assign values to I_var, beta0_var, gamma0_var and EF_var:
I_var<-I_tropical_Dentener
beta0_var<-beta0_tropical
gamma0_var<-gamma0_tropical
EF_var<-EF_tropical

parms_nonfixer <- list(I=I_var,beta0=beta0_var,gamma0=gamma0_var,EF=EF_var)
Ainit <- (I_var*w0*gamma0_var*m+I_var*gamma0_var*w0*phi-phi*beta0_var+phi*u0)/
  (w0*gamma0_var*(k+EF_var)*(m+phi))*0.5
Linit <- (beta0_var-u0)/(w0*gamma0_var*(m+phi))*0.75

op_nonfixer <- lsoda(c(Binit,Linit,Ainit,Einit),tv,N_model_ode_nonfixer,parms_nonfixer)
  
B0_nonfixer <- op_nonfixer[,2]
L_nonfixer <- op_nonfixer[,3]
A_nonfixer <- op_nonfixer[,4]
E_nonfixer <- op_nonfixer[,5]
  
CO2_nonfixer <- -(B0_nonfixer[length(tv)]-Binit)/10000/tf*44/12
N2O_nonfixer <- (E_nonfixer[length(tv)]-Einit)/10000/tf*44/28*298
netCO2N2O_nonfixer <- CO2_nonfixer+N2O_nonfixer
