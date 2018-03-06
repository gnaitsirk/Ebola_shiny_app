
library(shiny)
library(ggplot2)
library(dplyr)
library(deSolve)
library(tidyr)


ui <- fluidPage(
  titlePanel("Ebola Outbreak"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("Input.beta.I", "transmission rate in the community", 
                  min=0, max=0.2, value=0.084),
      sliderInput("Input.beta.H", "transmission rate in the hospital", 
                  min=0, max=.5, value=0.1134),
      sliderInput("Input.beta.F", "transmission rate at funerals", 
                  min=0, max=2, value=.929, step=.2)
    ),
    mainPanel(
      plotOutput("coolplot"),
      wellPanel(
        tags$body(h3("Stochastic Compartmental Model for Ebola Epidemic, Democratic Republic of Congo, 1995"), 
                  p(" S: Susceptible, E: Exposed, I: Infectious, H: Hospital, F:Funeral, R: Recovered, Inci: Incidence") 
            
    )
  )
)))



server <- function(input, output) {
  
  EBOLA=function(t,state,parameters){
    with(as.list(c(state,parameters)),{
      # rate of change
      dS=-S/N*(beta.I*I+beta.H*H+beta.F*FF);
      dE=S/N*(beta.I*I+beta.H*H + beta.F*FF) - alpha*E
      dI=alpha*E-I*(gamma.h*theta1+gamma.d*(1-theta1)*delta1+gamma.i*(1-theta1)*(1-delta1));
      dH=gamma.h*theta1*I-gamma.dh*delta2*H - gamma.ih*(1-delta2)*H
      dFF= gamma.dh*delta2*H + gamma.d*(1-theta1)*delta1*I - gamma.f*FF
      dR= gamma.f*FF + gamma.i*(1-theta1)*(1-delta1)*I + gamma.ih*(1-delta2)*H
      dInci=alpha*E
      # return the rate of change
      list(c(dS,dE,dI,dH,dFF,dR,dInci))
    }) # end with(as.list...)
  }
  
  
  ####################################################################################
  ## simulate the DRC 1995 outbreak
  ####################################################################################
  ## parameters for the DRC 1995 outbreak (Table 3 in Legrand et al. 2007)
  N=2e5; E0=H0=FF0=R0=0; I0=3; Inci0=3; S0=N-I0;
  state=c(S=S0,E=E0,I=I0,H=H0,FF=FF0,R=R0,Inci=Inci0);
  # no intervention
  alpha=1/7; # incubation period: 7 days
  gamma.h=1/5; # from onset to hopspitalization: 5 days
  gamma.d=1/9.6; # from onset to death: 9.6 days
  gamma.i=1/10; # from onset to end of infectiousness for survivors: 10 days
  gamma.f=1/2; # from death to traditional burial 2 days
  gamma.ih=1/(10-5); # from hospitalization to end of infectiousness for survivors
  gamma.dh=1/(9.6-5); # from hospitalization to death
  theta1=.67; # proportion infectious in the hospital
  delta1=.8; # CFR for unhospitalized
  delta2=.8; # CFR for hospitalize
  p=.85; ## reduced by a factor of p as we are running it deterministically
  #beta.I=.588/7; # transmission rate in the community
  #beta.H=.794/7; # transmission rate in the hospital
  #beta.F=7.653/7 * p # transmission rate at funerals; 
  times=0:280
  filtered <- reactive({
    
    
    parmsNoCtrl=c(alpha=alpha, # incubation period: 7 days
                  gamma.h=gamma.h, # from onset to hopspitalization: 5 days
                  gamma.d=gamma.d, # from onset to death: 9.6 days
                  gamma.i=gamma.i, # from onset to end of infectiousness for survivors: 10 days
                  gamma.f=gamma.f, # from death to traditional burial 2 days
                  gamma.ih=gamma.ih, # from hospitalization to end of infectiousness for survivors
                  gamma.dh=gamma.dh, # from hospitalization to death
                  theta1=theta1, # proportion infectious in the hospital
                  delta1=delta1, # CFR for unhospitalized
                  delta2=delta2, # CFR for hospitalize
                  beta.I=input$Input.beta.I, # transmission rate in the community
                  beta.H=input$Input.beta.H, # transmission rate in the hospital
                  beta.F=input$Input.beta.F # transmission rate at funerals, 
    );
    
    sim=ode(y=state,times=times,func=EBOLA,parms=parmsNoCtrl);
    sim<-as.data.frame(sim)
    sim %>% 
      select(time, S, E, I, H, FF, R, Inci) %>% 
      gather(Compartments, Cases, -time)
    
    
  })
  
  output$coolplot <- renderPlot({
    ggplot(filtered(), aes(x=time, y=Cases, 
                           group=Compartments, color=Compartments)) + 
      xlab("Time (Days)") + ylab("Number of People") + geom_line(alpha=.8) + theme_minimal() + 
      theme(legend.position="bottom") + theme(legend.title=element_blank())
    
    
  })
  
  
}

shinyApp(ui = ui, server = server)
