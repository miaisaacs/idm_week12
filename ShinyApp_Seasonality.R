## R shiny app to Simulate Disease Seasonality Using Temporally Forced Models
## and Test how seasonality changes with R0 and amplitude of seasonality
## Model: SEIR; Parameters based on childhood infections
## 4/15/18, by Wan Yang

library(shiny)
library(deSolve)


SEIRsin <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    beta = beta0 * (1 + beta1 * sin(omega * time +pi/4))
    
    dS = mu*N - beta * S * I / N - mu * S
    dE = beta * S * I / N - alpha * E - mu * E
    dI = alpha * E - mu * I - gamma * I
    
    list(c(dS, dE, dI))
  })
}
N=1; # 100%

ui <- fluidPage(h4('Simulating Disease Seasonality Using Temporally Forced Models\nTest how seasonality changes with R0 and amplitude of seasonality'),
                h5('Note: it may take a while to run... wait... wait... wait...'),
                # outputs
                sidebarLayout(
                  sidebarPanel(width=4,
                               sliderInput(inputId = 'R0',h6('R0'),value=17,min=6,max=20,step=1),
                               sliderInput(inputId = 'beta11',h6('Amplitude (β1) for Run 1'),value=0.02,min=0,max=0.25,step=.01),
                               sliderInput(inputId = 'beta12',h6('Amplitude (β1) for Run 2'),value=0.1,min=0,max=0.25,step=.01),
                               sliderInput(inputId = 'beta13',h6('Amplitude (β1) for Run 3'),value=0.23,min=0,max=0.25,step=.01),
                               h6('Note: other parameters based on childhood infections: \nlatent period = 8 days; infectious period = 5 days')
                               ),
                              
                              mainPanel(
                                plotOutput(outputId = 'plots',width = '100%', height = "550px")
                                )
                              )
                
                )

server <- function(input, output){
  
  output$plots=renderPlot({
    
    R0=input$R0         # R0 for all three runs
    beta11=input$beta11 # beta1 for run1
    beta12=input$beta12 # beta1 for run2
    beta13=input$beta13 # beta1 for run3
    
    # other parameters and initial conditions
    alpha=365/8; gamma=365/5; mu=1/50; omega=2*pi/1; # note time is in year
    beta0=R0*gamma
    N=1; # 100%
    E0=I0=N*1e-3; S0=N*6e-2; 
    state=c(S=S0, E=E0,I=I0);
    
    tm_step=7; # run it by week
    times=seq(0,1000,by=tm_step/365) # need to run for ~1000 yrs to have the 3rd sim stable!
    
    parms1=c(beta0=beta0,beta1=beta11,omega=omega,alpha=alpha,gamma=gamma,mu=mu);
    parms2=c(beta0=beta0,beta1=beta12,omega=omega,alpha=alpha,gamma=gamma,mu=mu);
    parms3=c(beta0=beta0,beta1=beta13,omega=omega,alpha=alpha,gamma=gamma,mu=mu);
    
    sim1=ode(y=state,times=times,func=SEIRsin,parms=parms1)
    sim2=ode(y=state,times=times,func=SEIRsin,parms=parms2)
    sim3=ode(y=state,times=times,func=SEIRsin,parms=parms3)
    
    ## plot the last 10 years
    tsim1=tail(sim1,round(365/tm_step*10,0)); # the last 10 yrs for the SEIRdem
    tsim2=tail(sim2,round(365/tm_step*10,0)); # the last 10 yrs for the SEIRsin
    tsim3=tail(sim3,round(365/tm_step*10,0)); # the last 10 yrs for the SEIRsin
    
  
    # plot results
    par(mfrow=c(3,1),cex=1,mar=c(3,3,1,1),oma=c(0,0,1,0),mgp=c(1.8,.5,0))
    plot(tsim1[,'time'],log(tsim1[,'I']),ylab='Log(fraction infective)',xlab='Time (year)',
         type='l',lwd=2, ylim=c(-15,-5))
    mtext(bquote(beta[1]==.(beta11)),side=3,line=-1.5,outer=F,adj=.95,cex=1.2)
    mtext(bquote(R[0]==.(R0)),side=3,line=0,outer=F,adj=.5,cex=1.5)
    abline(v=990:1000,col='grey',lty=2)
    plot(tsim2[,'time'],log(tsim2[,'I']),ylab='Log(fraction infective)',xlab='Time (year)',
         type='l',lwd=2, ylim=c(-15,-5))
    mtext(bquote(beta[1]==.(beta12)),side=3,line=-1.5,outer=F,adj=.95,cex=1.2)
    abline(v=990:1000,col='grey',lty=2)
    plot(tsim3[,'time'],log(tsim3[,'I']),ylab='Log(fraction infective)',xlab='Time (year)',
         type='l',lwd=2, ylim=c(-15,-5))
    mtext(bquote(beta[1]==.(beta13)),side=3,line=-1.5,outer=F,adj=.95,cex=1.2)
    abline(v=990:1000,col='grey',lty=2)
    
  })
  
}

shinyApp(ui=ui, server = server)