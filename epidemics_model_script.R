# This script shows the MTE functions and underlying assumptions for the temperature-dependent model parameters and temperature-independent parameters
# We then show how these parameters are used in the R0 calculation for the Daphnia magna - Ordospora colligata experimental system in Kirk et al 2020 Proc B
# Instead of calculating R0, the parameters could also be used in one's preferred simulator to simulate dynamics deterministically or stochastically #

require(deSolve) #loading the ODE solver
require(FME)
require(plotly)

minTemp  <- 273.15  # The minimum temperature 
maxTemp  <- 313.15  # The maximum temperature
resol <- 0.01  # Resolution
temp <- seq(from = minTemp, to = maxTemp, by = resol)


# Initialize vectors for each of the parameters
mu <- vector(mode="numeric", length(temp)) # natural mortality rate in Weibull survival
beta <- vector(mode="numeric", length(temp)) # ageing rate in Weibull survival
alpha <- vector(mode="numeric", length(temp)) # per-parasite virulence
chi <- vector(mode="numeric", length(temp)) # contact rate
lambda <- vector(mode="numeric", length(temp)) # shedding rate
sigma <- vector(mode="numeric", length(temp)) # probability of infection
theta <- vector(mode="numeric", length(temp)) # degradation rate
omega <- vector(mode="numeric", length(temp)) # shedding rate after death
equil.abundance<- vector(mode="numeric", length(temp)) # 


# Define MTE parameters 
# Note when parameterizing, Kirk et al. 2018 used 15C as reference T0, and Kirk et al 2019 used 20C as reference T0
chi0   =  1.304617e+02; sigma0= 5.684252e-09; omega0= 131; beta0 <- 2.281060732; mu0 <- 0.007507025 # Normalization constant (i.e.. value at reference temperature)                           
E_chi  =  2.90155e-01; E_sigma = 1.106441; E_omega= -0.295; E_beta <- 0.390532437; E_mu <- 0.801482243   # Activation energy                        
E_sigmaL = 6.766889; E_omegaL= 357; E_muL <- E_mu*5    # Inactivation energy at lower threshold   
E_chiH =  4.129647; E_sigmaH = 5.66131; E_omegaH= 6.32; E_betaH <- 1.472373114; E_muH <- 4.239515185  #  Inactivatio' energy at upper threshold      
T0_chi = 293.15; T0_sigma = 293.15; T0_omega= 288.15; T0_beta <- 288.15 ; T0_mu <- 288.15  # Reference temperatures (taken as 20C for contact and infection, 15 for mu, beta, and omega)                        
T_sigmaL = 273.15+1.414539e+01;  T_omegaL= 273.15+9.9; T_muL <-  8.957624193 + 273.15# Lower thresholds                          
T_chiH = 273.15+3.513025e+01  ; T_sigmaH = 273.15+2.343956e+01  ; T_omegaH= 273.15+23.9 ; T_betaH <- 28.8419577 +273.15; T_muH <-  30.23264714 +273.15 #Upper thresholds                        
cs = 2.574784e-01 ; is = -1.181799; # allometric scaling values for contact rate and infection rate (all other parameters don't use body size here)


k <- 8.62*10^(-5) # Boltzmann's constant


# Assumptions here #
# 1) Within-host parasite abundance as a function of equilibrium abundance
abundance.prop = 0.182

# 2) How long the clusters take to burst 
burst.time = 7

# 3) How many spores are in each cluster
spores.per.cluster = 24

# 4) How many of the spores are released into the environment from each burst cluster 
burst.size = 0.5

# 5) Daphnia size
size = 2700
mass <-  (0.009*(size*0.001)^2.63)*0.001



# Calculate parameters across temperature #

# Note that sigma in this loop is the infection RATE, will calculate infection probability in step after #
for(t in seq(temp) ){
  chi[t] <- ((((chi0*(mass^cs)*exp((-E_chi/k)*(1/temp[t]-1/T0_chi))/(1+exp(E_chiH/k*(-1/temp[t]+1/T_chiH)))))/3.5E7)*1440); # values at the end are to standardize to the volume and time of the experiment
  mu[t] <-     mu0 * exp((-E_mu/k)  * (1/temp[t]-1/T0_mu))* (1+ exp(E_muL/k * (1/temp[t]-1/T_muL))  +  exp(E_muH/k * (-1/temp[t]+1/T_muH)));
  beta[t] <- beta0*exp((-E_beta/k) * ((1/temp[t])-(1/T0_beta)))*((1+exp((E_betaH/k)*((-1/temp[t])+(1/T_betaH))))^(-1))
  sigma[t] <-  sigma0*(mass^is)*exp((-E_sigma/k) *  (1/temp[t]-1/T0_sigma)) * (1+  exp(E_sigmaL/k  * (1/temp[t]-1/T_sigmaL)) +  exp(E_sigmaH/k * (-1/temp[t]+1/T_sigmaH)))^(-1); 
  equil.abundance[t] <-  omega0  * exp((-E_omega/k) *  (1/temp[t]-1/T0_omega)) * (1+  exp(E_omegaL/k  * (1/temp[t]-1/T_omegaL)) +  exp(E_omegaH/k * (-1/temp[t]+1/T_omegaH)))^(-1)
  omega[t] <- (equil.abundance[t]*abundance.prop)*spores.per.cluster
  alpha[t] <- (equil.abundance[t]*abundance.prop)*5.12E-6; # parasite intensity multiplied by per parasite virulence
  lambda[t] <- (equil.abundance[t]*abundance.prop)*(spores.per.cluster*burst.size)/burst.time
}



# Now to calculate infection probablity, need to calculate contact in microliters/min to calculate gut residence time #
standardized.contact <- c()

# contact here is in microliters / min
for(i in 1:length(temp)){
  standardized.contact[i] <- ((chi0*(mass^cs)*exp((-E_chi/k)*(1/temp[i]-1/T0_chi))/(1+exp(E_chiH/k*(-1/temp[i]+1/T_chiH)))))
}

# Calculate ingestion rate (cells/hour), and we have 10 algal cells per microliter 
ing.cells.rate <- ((standardized.contact)*10)*60

# From Kooijman 1993 and Evers and Kooijman: volume of gut in cells = gv (9900) * length^3 (mm^3) / ingested cells rate, this gives gut.res time in hours
# Note that the concentration of algae here was much lower than in Kirk et al. 2019, therefore ingeston of cells is lower and gut residence time is longer
gut.res.time <- (9900*((size/1000)^3))/ing.cells.rate

#probability of infection 
sigma.prob <- 1 - exp(-sigma*gut.res.time)




#### Parameters ###
# Add in temperature independent parameters, and then subset continuous temperature functions from above for our 8 temperatures from 10 - 13.5 

# Input from stocks #
inputS = 3 + 0.5350318
inputI = 0.4649682

# Max recruitment rate #
recruit = 1.33

# Carrying capacity #
kappa = 170

# Harvesting rate #
h = 0.0235

# Degradation of dead infecteds rate #
theta = 0.1
  
# Environmental parasite mortality rate (induced via changing the medium)
gamma <-  0.02857143



# In Kirk et al. 2018, mortality rate was able to change over time (defined by a Weibull survival model) whereas here we use a constant hazard rate (i.e., exponential survival distribution)
# To translate between the two, for each temperature we will simulate survival using the mu (mortality scale parameter) and beta (mortality shape parameter)
# And then use this to find average lifespan and then average mortality rate #

hazard <- matrix(nrow=4001,ncol=2)
hazard[,1] <- temp

for(i in c(1:4001)){
  mu.sim <- mu[i]
  beta.sim <- beta[i]
  
  U<-1 ; #Starting conditions
  N0<-c(U)
  
  TT<-seq(0.1,400,0.1)   #Set amount of time steps to run in simulation
  step <- 0.1
  
  parms<-c(beta.sim,mu.sim)
  
  survival<-function(t,y,p){
    beta<-p[1]; mu<-p[2];
    U<-y[1];
    
    dU <- -beta*mu^(beta)*(t^(beta-1))*U
    
    list(c(dU))
  }
  
  #run lsoda	
  predictions <-lsoda(N0,TT,survival,parms)	
  hazard[i,2] <- 1/predictions[which(predictions[,2] <= 0.5)[[1]],1]
  
}

# Set mu.t to be equal to the hazard, which is now mu for an exponential survival distribution #
mu.t <- hazard[,2]



# Show temperature-dependent relationships that will then be used in R0 calculation #
quartz()
par(mfrow=c(2,2))
par(mar=c(1,5,1,1))
par(oma=c(4,0,0,0))
plot(mu.t, type="l", lwd=3, ylim=c(0,0.5), xlab=NA, xaxt="n",ylab="Mortality Rate (1/day)", cex.lab=1)
axis(side=1, at=c(0,500,1000,1500,2000,2500,3000,3500,4000),labels=c(0,5,10,15,20,25,30,35,40),tick=TRUE)

plot(chi, type="l", lwd=3, xlab=NA,xaxt="n",ylab="Contact Rate (1/day)",cex.lab=1)
axis(side=1, at=c(0,500,1000,1500,2000,2500,3000,3500,4000),labels=c(0,5,10,15,20,25,30,35,40),tick=TRUE)


plot(sigma.prob, type="l", lwd=3, xaxt="n",ylab="Probability of Infection",cex.lab=1)
axis(side=1, at=c(0,500,1000,1500,2000,2500,3000,3500,4000),labels=c(0,5,10,15,20,25,30,35,40),tick=TRUE)

plot(equil.abundance*abundance.prop, xaxt="n",type="l", lwd=3, ylab="Within-Host Infection Intensity",cex.lab=1)
axis(side=1, at=c(0,500,1000,1500,2000,2500,3000,3500,4000),labels=c(0,5,10,15,20,25,30,35,40),tick=TRUE)
mtext(expression(paste('Temperature ('~degree~'C)')),side=1, cex=1.5, line = 3, outer=TRUE)





### Now that all the parameters are defined, we can calculate R0 across temperature at experimental density of 170 hosts #
R0 <- c()
for(t in 1:length(temp)){  
  R0[t] <- (((lambda[t] / (mu.t[t] + alpha[t]+h)) + (omega[t]*((mu.t[t]+alpha[t])/(mu.t[t]+alpha[t]+h)))) * (chi[t]*sigma.prob[t]*170/gamma))
}

# This crosses R0 = 1 at 11.97 degrees C #
quartz()
par(mfrow=c(1,1))
plot(R0 ~ c(temp-273.15), type="l", lwd=3, ylab="R0", xlab="Temperature (C)")
abline(h=1, col="red", lty=2) # red dashed line shows R0 = 1 


### We can also calculate R0 at different temperatures and different host densities #
R0.matrix.thin <- matrix(nrow=2501, ncol=201) # Look across 25 degrees and 0-200 hosts #

for(t in 501:3001){  
  for(i in 1:201){
    density = seq(0,200,1)[i]
    R0.matrix.thin[t-500,i] <- (((lambda[t] / (mu.t[t] + alpha[t]+h)) + (omega[t]*((mu.t[t]+alpha[t])/(mu.t[t]+alpha[t]+h)))) * (chi[t]*sigma.prob[t]*density/gamma))
  }
}

temp.thin <- temp[501:3001]
temps.thin <- temp.thin-273.15
density <- seq(0,200,1)
mycols = c("navyblue","green","yellow","orange","firebrick4")

quartz()
filled.contour(x=temps.thin, y=density, z=R0.matrix.thin, xlab="Temperature", 
               ylab="Density of susceptible large adult females", 
               main="R0",levels=c(0,1,6,11,16,21),
               col=mycols)

