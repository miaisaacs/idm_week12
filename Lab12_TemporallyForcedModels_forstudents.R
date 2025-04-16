## Lab: Temporally forced models

library(deSolve)

########################################################################
## Part 0 (pre-lab): Play with the cosine/sine function in R
########################################################################
## first check out the sinusoidal forcing
omega=2*pi/365 # the frequency
times=seq(1,6*365,by=1)  # time steps
y=cos(omega*times); # cos() is the cosine function in R; use sin() for sine
# change to sin if you'd like
# y=sin(omega*times); 

# plot the sinusodal function
par(cex=1.2,mar=c(3,3,1,1),mgp=c(2,.5,0))
plot(times,y,type='l',lwd=2,xlab='Time (year)',ylab='cos(wt)')
abline(h=0,col='grey30')  # add a horizontal line at 0
abline(h=c(-1,1),lty=2,col='grey30')  # add horizontal lines at -1 and 1 (lower and upper bounds)
abline(v=seq(1,365*6,by=365),lty=2,col='grey30') # vertical lines to show the year divisions




########################################################################
## Part 1: Compare models with vs. without temporal forcing
########################################################################
## SEIR model without seasonal forcing
SEIRdem = function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS = mu*N - beta * S * I / N - mu * S
    dE = beta * S * I / N - alpha * E - mu * E
    dI = alpha * E - mu * I - gamma * I
    list(c(dS, dE, dI))
  })
}


## SEIR model with sinusoidal seasonal forcing (SEIRsine)
## CODE ON YOUR OWN

SEIRsine <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    beta = beta0 * (1 + beta1 * cos(omega * time))
    dS = mu*N - beta * S * I / N - mu * S
    dE = beta * S * I / N - alpha * E - mu * E
    dI = alpha * E - mu * I - gamma * I
    
    list(c(dS, dE, dI))
  })
}



# Parameters and inital conditions
N=1; # population size, in fraction, i.e. 100%
E0=N*1e-3; # initial exposed
I0=N*1e-3; # initial infectious
S0=N*6e-2; # initial susceptible
R0=17; # R0
alpha=1/8; # 1/latent period, per day
gamma=1/5; # 1/infection period, per day
mu=1/50/365; # mortality rate, per day
omega=2*pi/365; # frequency of seasonal forcing
beta1 = 0.02 # amplitude of seasonal forcing


state=c(S=S0, E=E0,I=I0);

parms0=c(beta=R0*gamma,alpha=alpha,gamma=gamma,mu=mu); # parameters for the SEIRdem model
parms1=c(beta0=R0*gamma,beta1=beta1,omega=omega,alpha=alpha,gamma=gamma,mu=mu); # for the SEIRsine model
times=seq(1,1000*365,by=1) # need to run 1000 yrs to have the 3rd sim stable!
## NOTE IT MAY TAKE A FEW MINUTES TO RUN EACH SIMULATION FOR 1000 YEARS

# put together the code to run the model 
# and plot the last 10 years here on your own:

# EXAMPLE CODE:
if(F){
  
  sim0=ode(y=state,times=times,func=SEIRdem,parms=parms0)
  sim1=ode(y=state,times=times,func=SEIRsine,parms=parms1)
  
  
  
  ## compare the two
  ## get the last 10 years using the tail function
  tsim0=tail(sim0,365*10); # the last 10 yrs for the SEIRdem
  tsim1=tail(sim1,365*10); # the last 10 yrs for the SEIRsine

}


sim0=ode(y=state,times=times,func=SEIRdem,parms=parms0)
sim1=ode(y=state,times=times,func=SEIRsine,parms=parms1)

tsim0=tail(sim0,365*10); # the last 10 yrs for the SEIRdem
tsim1=tail(sim1,365*10); # the last 10 yrs for the SEIRsine


par(mfrow=c(1,1),cex=1,mar=c(3,3,1,1),mgp=c(1.8,.5,0))
plot(tsim1[,'time'],tsim1[,'I'],ylab='%I',xlab='Time (days)',type='l',lwd=2)
lines(tsim0[,'time'],tsim0[,'I'],lty=2)
legend('topright',c('With seasonal forcing', 'No seasonal forcing'),
       lty=c(1,2),cex=.9,bty='n')





########################################################################
## Part 2: Set up and compare different seasonal forcing
########################################################################
# School term time forcing models, without the correction
SEIRterm <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # real-tim beta based on school term-time forcing (uncorrected) 
    beta= beta0 * (1 + b.term * Term[time])
    
    dS = mu*N - beta * S * I / N - mu * S
    dE = beta * S * I / N - alpha * E - mu * E
    dI = alpha * E - mu * I - gamma * I
    
    list(c(dS, dE, dI))
  })
}

# School term time forcing models, with the correction
SEIRterm.corrected <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # real-tim beta based on school term-time forcing (corrected) 
    correction.factor=(1/365*((1+b.term)*273+(1-b.term)*92))
    beta= beta0/correction.factor * (1 + b.term * Term[time])
    
    dS = mu*N - beta * S * I / N - mu * S
    dE = beta * S * I / N - alpha * E - mu * E
    dI = alpha * E - mu * I - gamma * I
    
    list(c(dS, dE, dI))
  })
}

correction.factor <- (1/365) * ((1 + b.term) * 273 + (1 - b.term) * 92)

beta <- function(day) {
  beta0 / correction.factor * (1 + b.term * Term[day])
}

beta(1)
beta(12)
beta(50)


## SETTING UP THE SCHOOL TERM-TIME FUNCTION:
# term-time forcing
holidays=c(1:6, 100:115, 200:251, 300:307, 356:365, 0);  # the days in a year that are holiays
times=seq(1,100*365); # in day
Term=rep(1,length(times));  # initialize a vector to store the Term
ind=(1:length(Term) %% 365) %in% holidays  # find those days that are school holidays
Term[ind]=-1; # set them to -1



# Parameters and initial conditions:
N=1; # Note the state variables (N, S, E, I, etc) here are fractions, e.g. here 1=100%
S0=6e-2*N; E0=I0=1e-5*N
state=c(S=S0, E=E0,I=I0);

beta0=3; b.term=.5;
alpha=1/8; gamma=1/5; mu=1/50/365; 


# Run the uncorrected school term-time forcing
parmsTerm=c(beta0=beta0,b.term=b.term,alpha=alpha,gamma=gamma,mu=mu);
simTerm=ode(y=state,times=times,func=SEIRterm,parms=parmsTerm)

# Run the corrected school term-time forcing
parmsTerm=c(beta0=beta0,b.term=b.term,alpha=alpha,gamma=gamma,mu=mu);
simTermCorrected=ode(y=state,times=times,func=SEIRterm.corrected,parms=parmsTerm)

# Run the Sinusiodal function:
omega=2*pi/365; 
parmsSine=c(beta0=beta0,beta1=.3,alpha=alpha,gamma=gamma,mu=mu);
simSine=ode(y=state,times=times,func=SEIRsine,parms=parmsSine)


# plot the first and last 10 yrs together for the SEIRterm (no correction)
par(mfrow=c(2,1),mar=c(3,3,1,1),mgp=c(1.5,.5,0),cex=.9)
plot(head(simTerm[,c('time','I')],365*10),type='l',xlab='Time (days)',ylab='%I',main='First 10 years')  # first 10 yrs
plot(tail(simTerm[,c('time','I')],365*10),type='l',xlab='Time (days)',ylab='%I',main='Last 10 years')  # last 10 yrs

# plot results from the two school term-time simulations and compare
# to do so, first we put the results together using the cbind function:
resI=cbind(simTerm[,c('time','I')],simTermCorrected[,'I']); colnames(resI)=c('time','I.term','I.termCorrected')
# now plot:
par(mfrow=c(2,1),mar=c(3,3,1,1),mgp=c(1.5,.5,0),cex=.9)
# plot first 10 years [Note: the model is still spinning up]
matplot(head(resI[,'time'],365*10),head(resI[,c('I.term','I.termCorrected')],365*10),type='l',lty=1,xlab='Time (days)',ylab='%I',col=c('blue','red'),main='First 10 years')  # first 10 yrs
legend('topright',c('Term-time: Uncorrected','Term-time: Corrected'),col=c('blue','red'),cex=.8,lty=1,bty='n')
# plot the last 10 years - more stablized
matplot(tail(resI[,'time'],365*10),tail(resI[,c('I.term','I.termCorrected')],365*10),type='l',lty=1,xlab='Time (days)',ylab='%I',col=c('blue','red'),main='Last 10 years')  # first 10 yrs
legend('topright',c('Term-time: Uncorrected','Term-time: Corrected'),col=c('blue','red'),cex=.8,lty=1,bty='n')


# plot results from all three simulations and compare
# to do so, first we put the three results together using the cbind function:
resI.all3=cbind(simTerm[,c('time','I')],simTermCorrected[,'I'],simSine[,'I']); colnames(resI.all3)=c('time','I.term','I.termCorrected','I.sine')
par(mfrow=c(2,1),mar=c(3,3,1,1),mgp=c(1.5,.5,0),cex=.9)
# plot first 10 years [Note: the model is still spinning up]
matplot(head(resI.all3[,'time'],365*10),head(resI.all3[,c('I.term','I.termCorrected','I.sine')],365*10),type='l',lty=1,xlab='Time (days)',ylab='%I',col=c('blue','red','orange'),main='First 10 years')  # first 10 yrs
legend('topright',c('Term-time: Uncorrected','Term-time: Corrected','Sinusoidal'),col=c('blue','red','orange'),cex=.8,lty=1,bty='n')
# plot the last 10 years - more stablized
matplot(tail(resI.all3[,'time'],365*10),tail(resI.all3[,c('I.term','I.termCorrected','I.sine')],365*10),type='l',lty=1,xlab='Time (days)',ylab='%I',col=c('blue','red','orange'),main='Last 10 years')  # first 10 yrs
legend('topright',c('Term-time: Uncorrected','Term-time: Corrected','Sinusoidal'),col=c('blue','red','orange'),cex=.8,lty=1,bty='n')





########################################################################
## Part 3: Test how the epidemic dynamics alter with varying seasonal forcing
########################################################################
# Use the app (ShinyApp_Seasonality.R)

