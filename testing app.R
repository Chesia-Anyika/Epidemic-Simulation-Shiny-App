#############################################################################
# STA3040-Endsem 
# Chesia Anyika
# ID - 665567
############################################################################

#############################################################################
#
# Necessary Libraries and Defining states
#
############################################################################
library(shiny)
library(ggplot2)
library(dplyr)
library(sna)
library(knitr)
library(reshape2)
library(igraph)
#library(animation)

STATES <- 7

STATENAMES <-  c("Unexposed",
                 "Asymptomatic & contagious",
                 "Symptomatic and contagious",
                 "Symptomatic and not contagious",
                 "Post-COVID immune",
                 "Naturally immune",
                 "Death")


STATELABELS <-  c("Unexposed","Asymptomatic\n & contagious",
                  "Symptomatic \n& contagious",
                  "Symptomatic \n& not contagious",
                  "Post-COVID immune",
                  "Naturally immune",
                  "Death")


#############################################################################
#
# Define the agent, and transition matrix
#
############################################################################

makeAgent <- function(psychstate,biostate,age=30)
{
  
  return (list(psychstate=psychstate,
               biostate=biostate,
               age=age,
               nextbiostate=NA,
               biostatecountdown=NA))
}


# * 1. Unexposed
# * 2. Asymptomatic but infected/contagious
# * 3. Symptomatic and contagious
# * 4. Symptomatic and not contagious
# * 5. Post-COVID Immune
# * 6. Naturally immune (will not contract)
# * 7. Death


bioTransition <- matrix(0,STATES,STATES)
bioMin <- matrix(1,STATES)      #state time minimum
bioMax <- matrix(1,STATES)      #state time maximum



bioMin[2] <- 3            #infected but asymptomatic for 3 to 15 days
bioMax[2] <- 15          
bioTransition[2,3] <- .15  #transition to infected with symptoms
bioTransition[2,5] <- .85  #transition to no longer contagious/cured


bioMin[3] <-    3             #symptoms + contagion
bioMax[3] <- 8                #symptoms + contagion max
bioTransition[3,4] <-  .95    #transition to no longer contagious
bioTransition[3,7] <-  .05    #transition to death state 


bioMin[4] <- 1          #symptoms but no longer contagious
bioMax[4] <- 7
bioTransition[4,5] <- 1  #Transition to 'immune' cured state.


setAgentState<- function(agent, biostate)
{
  agent$biostate <- biostate
  if(sum(bioTransition[biostate,])>0) # this state transitions to something else.
  {
    ##which state do we go to?
    agent$biostatecountdown <- sample(x=seq(bioMin[biostate],bioMax[biostate]),1) #how long will we state in this state?
    agent$nextbiostate <- sample(1:STATES, prob=bioTransition[agent$biostate,],size=1)
    
  } else{
    agent$biostatecountdown <- NA
    agent$nextbiostate <- NA   ##just so we can tell if the agent is finished.
  }
  return(agent) 
}

transitionAgent<- function(agent)
{
  return(setAgentState(agent,agent$nextbiostate))
}

updateAgent<- function(agent)
{
  if(!is.na(agent$biostatecountdown))
  {
    agent$biostatecountdown <- agent$biostatecountdown -1
    if(agent$biostatecountdown <=0)  ##new state
    {
      agent <- transitionAgent(agent)
      
    }
  }
  return(agent)
}

#############################################################################
#
# Rudimentary biotransition graph
#
############################################################################


# Convert the bioTransition matrix to a graph object
bioTransitionGraph <- function() {
  g <- graph_from_adjacency_matrix(bioTransition, mode = "directed", weighted = TRUE)
  E(g)$label <- E(g)$weight
  V(g)$label <- STATENAMES
  return(g)
}


#############################################################################
#
# Defining networks
#
############################################################################

makeNetwork<- function(numAgents,numsets=3,steps=1,power=1)
{
  ord <- sample(numAgents)
  tmp<-as_adjacency_matrix(sample_pa(numAgents,power=power),sparse=F)
  tmp <- tmp + t(tmp)
  tmp <- tmp[ord,ord]
  if(numsets>1)
  {
    for(i in 2:numsets)
    {
      ord <- sample(numAgents)
      sn2 <-as_adjacency_matrix(sample_pa(numAgents,power=power),sparse=F)[ord,ord]
      tmp <- tmp + sn2 + t(sn2)
      
    }
  }
  if(steps>1)
  {
    for(i in 2:steps)
    {
      tmp <- tmp + tmp %*% tmp  ## iterate to neighbors
    }
  }
  (tmp>0)+0
}

## This tries to make a network with families, schools, workplaces..
## 
makeDemographicNetwork <- function(numAgents=1000, 
                                   householdsize=2.5,
                                   numSchools=4,
                                   numWorkplaces=25,
                                   agePyramid = c(20,8,10,12,14,13,12,11))
  
{
  
  numChildren <- floor(sum(agePyramid[1:2])/sum(agePyramid) * numAgents)
  numAdults <- numAgents-numChildren
  
  ##Adult ages age:
  
  adultAge<- sample(c(3:8)*10-5,size=numAdults,replace=T,prob=agePyramid[-(1:2)])
  childrenAge <- sample(c(5:15),size=numChildren,replace=T,prob=c(1:11))#agePyramid[1:2])
  ages <- c(adultAge,childrenAge)
  ##Each household
  numHouseholds <-floor( numAgents/householdsize)
  householdXY <- matrix(runif(numHouseholds*2)*10,ncol=2)
  
  #Now, let's assign people to households.  Adults can be assigned to anywhere;
  ##children need to be assigned to a household that already exists.
  
  ##assign 'head-of-household' ##first the adults, then the children
  myHousehold <- rep(NA,numAgents)
  myHousehold[1:numHouseholds] <- 1:numHouseholds
  myHousehold[(numHouseholds+1):numAgents] <-  sample(numHouseholds,replace=T,size=numAgents-numHouseholds)
  
  
  ##now, all households are assigned.  Let's assign workplaces and schools.  
  ##for simplicity, every adult has a workplace, and every child has a school, even though children
  ## or adults may work from home/too young for school, or be retired,
  ##wokplaces 1 is hosptial
  ## workplace 1+1;numschools are considered schools
  
  ##Workplace size is roughly zipfs-law
  workPlaces <-sample(1:numWorkplaces,size=numAdults,prob=1/(1+1:numWorkplaces),replace=T)
  schoolPlaces <- sample(1+(1:numSchools),size=numChildren,replace=T)
  work <- c(workPlaces,rep(0,numChildren))
  school <- c(rep(0,numAdults),schoolPlaces)
  pool <- list()
  for(i in 1:numAgents)
  {
    tmpAgent <- makeAgent(psychstate=1,
                          biostate=1,
                          age=ages[i])
    
    tmpAgent$household = myHousehold[i]
    tmpAgent$work = work[i]  ##work or school
    tmpAgent$school=school[i]
    tmpAgent$head <- (i<numHouseholds)
    tmpAgent$xy <- householdXY[myHousehold[i] ]
    pool[[i]] <-tmpAgent
  }
  
  
  ##now, let's lay out the social network
  schoolmat <- (outer(school,school,"==") * outer(school,school,"*"))>0
  workmat <- (outer(work,work,"==") * outer(work,work,"*"))>0
  housemat <- ( outer(myHousehold,myHousehold,"=="))>0
  agentXY <- householdXY[myHousehold,] +  matrix(.15* rnorm(length(myHousehold)*2),ncol=2)
  neighborhood <- as.matrix(dist(agentXY)) < runif(nrow(agentXY)^2)*2
  #  neighborhood <- outer(agentXY,agentXY, function(a,b){(a-b)^2 < runif(1)*.5})
  network <- schoolmat+ workmat + housemat + neighborhood
  return (list(xy=agentXY, pool=pool,network=network,
               neighborhood=neighborhood,
               worknet=workmat,
               schoolnet=schoolmat,
               housenet=housemat))
}

#############################################################################
#
# Plotting the networks
#
############################################################################
mygplot <- function(coord, network,states,main="",edgecol="grey40",add=F)
{
  if(is.null(coord))
  {
    coord  <- gplot.layout.fruchtermanreingold(network,layout.par=list(niter=500))
  }
  
  newmin <- mean(coord[,2]) - (-min(coord[,2]) + mean(coord[,2])) * 1.4
  palette=c("white","yellow","red","green","darkgreen","blue","black")
  if(add==F)
    plot(coord,col="black",bty="n",pch=16,cex=2.7,xaxt="n",yaxt="n",main=main,
         xlab="",ylab="",axes=F,
         ylim=c(newmin,max(coord[,2])),type="n")
  
  for(i in 1:nrow(network))
  {
    if(sum(network[i,])>0)
    {
      segments(coord[i,1],
               coord[i,2],
               coord[network[i,]>0,1,drop=F],
               coord[network[i,]>0,2,drop=F],col=edgecol)
    }
  }
  points(coord,pch=16,cex=2.3,col= palette[states])
  
  points(coord,pch=1,cex=2.3,col="black")
  legend(mean(coord[,1]),min(coord[,2]),bty='n',y.intersp=.7,cex=.8,
         STATENAMES, pch=16,col=palette)
  
  return (coord)
}

net <- makeDemographicNetwork(500, numSchools=4,
                              numWorkplaces=50)


#############################################################################
#
# USER INTERFACE
#
############################################################################
# Define UI
ui <- fluidPage(
  titlePanel("Modeling Psychological Impacts on Epidemic Spread"),
  
  # Adding Background Explanation
  tags$div(class="background-explanation",
           tags$p("This is a Shiny App Interface for codes borrowed from Lamia Alam & Shane Mueller on their research work on Modeling psychological impacts on epidemic spread while accounting for the demographic networks."),
           tags$p("They try to incorporate opinion dynamics, network modeling, and a bit of the game theory to understand epidemic spread.")
  ),
  
  tabsetPanel(
    tabPanel("COVID-19 State Transitions", 
             HTML("The model defines an agent, and 7 states for the agent's progression through the COVID-19 disease. <br><br>
       These states are namely:"),
             tags$ul(
               tags$li("Unexposed"),
               tags$li("Asymptomatic but infected/contagious"),
               tags$li("Symptomatic and contagious"),
               tags$li("Symptomatic and not contagious"),
               tags$li("Post-COVID Immune"),
               tags$li("Naturally immune (will not contract)"),
               tags$li("Death")
             ),
             tags$p("The below graph represents the biological transitions between various states in the disease's spread."),
             sidebarLayout(
               sidebarPanel(
                 helpText("Visualization of COVID-19 state transitions using a directed graph."),
                 actionButton("submit1", "Generate Graph") # Changed submit to submit1 to avoid ID conflict
               ),
               mainPanel(
                 plotOutput("bioTransitionPlot", width = "100%", height = "800px")
               )
             )),
    tabPanel("Demographic Networks",
             # Adding explanation text
             HTML("Next, a network using preferential treatment to represent the social world is made. Additionally,agents are split into age groups as follows:"),
             tags$ul(
               tags$li("Children - agents aged from 5 years to 15 years"),
               tags$li("Adults - agents aged from 25 years to 75 years")
             ),
             HTML("4 distinct networks are made and added together to make a final fifth. <br><br> These models are namely:"),
             tags$ul(
               tags$li("The Household Network - incorporates both Children and Adults"),
               tags$li("The School Network - incorporates only Children"),
               tags$li("The Work Network - incorporates only Adults"),
               tags$li("The Neighbourhood Network - incorporates both children and Adults"),
               tags$li("The Combined Network - a combination of all 4 Networks")
             ),
             tags$p("The diverse networks are made to represent how we each have different types of relationships and interactions-- you may be the president of your drama club and so are a central member, but you are also in a pottery class where you don't talk to very many people (although your teacher does)"),
             sidebarLayout(
               sidebarPanel(
                 selectInput("networkType", "Choose a Network Type:",
                             choices = c("Household", "Schools", "Workplaces", "Neighborhood", "Combined"))
               ),
               mainPanel(
                 plotOutput("networkPlot", width = "100%", height = "800px")
               )
             )),
    tabPanel("Day by Day Spread of Infection Simulation",
             # Adding explanation text
             tags$p("Lastly, taking into account the preferential networks created, this visualisation shows a simulation of the spread of the COVID-19 Pandemic over a 100 day period"),
             sidebarLayout(
               sidebarPanel(
                 sliderInput("numDays", "Number of Days", min = 1, max = 100, value = 10),
                 actionButton("runSimulation", "Run Simulation")
               ),
               mainPanel(
                 plotOutput("epidemicPlot")
               )
             )
    )
  )
)


#############################################################################
#
# USER INTERFACE
#
############################################################################

# Define combined server logic
server <- function(input, output) {
  
  # App 1: COVID-19 State Transitions Plot
  graphData <- eventReactive(input$submit1, {  # Note the changed ID to submit1
    bioTransitionGraph()  # Assumes bioTransitionGraph() is defined elsewhere
  }, ignoreNULL = FALSE)
  
  output$bioTransitionPlot <- renderPlot({
    g <- graphData()
    plot(g, edge.label=E(g)$label, vertex.size=20, vertex.label.color="black",
         vertex.label.cex=0.8, edge.arrow.size=0.5)
  })
  
  # Demographic Networks Plot
  output$networkPlot <- renderPlot({
    networkType <- input$networkType
    if(networkType == "Household") {
      mygplot(coord=net$xy, net$housenet, rep(1, nrow(net$network)), main="Households", edgecol="blue")
    } else if(networkType == "Schools") {
      mygplot(coord=net$xy, net$schoolnet+0, rep(1, nrow(net$network)), main="Schools", edgecol="orange")
    } else if(networkType == "Workplaces") {
      mygplot(coord=net$xy, net$worknet+0, rep(1, nrow(net$network)), main="Workplaces", edgecol="darkgreen")
    } else if(networkType == "Neighborhood") {
      mygplot(coord=net$xy, net$neighborhood, rep(1, nrow(net$network)), main="Neighborhood", edgecol="brown")
    } else if(networkType == "Combined") {
      mygplot(coord=net$xy, net$network, rep(1, nrow(net$network)), main="Combined network")
      mygplot(coord=net$xy, net$housenet, rep(1, nrow(net$network)), main="Household network", edgecol="blue", add=TRUE)
      mygplot(coord=net$xy, net$schoolnet+0, rep(1, nrow(net$network)), main="School network", edgecol="orange", add=TRUE)
      mygplot(coord=net$xy, net$worknet+0, rep(1, nrow(net$network)), main="Workplace network", edgecol="darkgreen", add=TRUE)
      mygplot(coord=net$xy, net$neighborhood, rep(1, nrow(net$network)), main="Neighborhood network", edgecol="brown", add=TRUE)
    }
  })
  
  # App 2: Epidemic Simulation Plot
  output$epidemicPlot <- renderPlot({
    # Your simulation code here
    numAgents <- 1000
    numDays <- input$numDays
    naturalImmunity <- 0
    net <- makeDemographicNetwork(numAgents)
    
    pool <- net$pool
    cc <- net$xy
    socialnetwork <- net$network
    
    
    numInteractions <- rep(8, numDays)  ##how many interactions per day per agent on average?
    contagionProb <- rep(0.03, numDays)  ##normal contagion probability after contact
    sampleFromNetwork <- rep(1.0, numDays)  ##how likely you are to stick with 'your' network
    
    plotNetwork <- TRUE
    
    if(plotNetwork) {
      cc <- mygplot(coord = cc, socialnetwork, rep(1, nrow(socialnetwork)), main = "Initial state")
    }
    
    disthistory <- matrix(NA, ncol = 7, nrow = numDays)
    
    pool <- list()
    for(i in 1:numAgents) {
      pool[[i]] <- makeAgent(psychstate = 1,
                             biostate = sample(c(1, 6), p = c(1 - naturalImmunity, naturalImmunity), 1))
    }
    
    ##infect patient 0
    numInfected <- 5
    for(i in sample(numAgents, numInfected)) {
      pool[[i]] <- setAgentState(pool[[i]], 2) ##infect this person
    }
    
    for(day in 1:numDays) {
      
      # who are you going to talk to today.
      sneezers <- rep(1:numAgents, each = numInteractions[day])
      sneezedons <- rep(NA, length(sneezers))
      
      for(i in 1:length(sneezers)) {
        if(runif(1) < (1 - sampleFromNetwork[day])) {
          sneezedons[i] <- sample(numAgents, 1)
        } else {
          sneezedons[i] <- sample(1:numAgents, prob = socialnetwork[sneezers[i], ], 1)
        }
      }
      
      for(i in 1:length(sneezers)) {
        agent1 <- pool[[ sneezers[i] ]]
        agent2 <- pool[[ sneezedons[i] ]]
        
        ##this constitutes the rules of infection.
        if((agent1$biostate == 2 || agent1$biostate == 3 ) & agent2$biostate == 1 & runif(1) < contagionProb[day]) {
          pool[[ sneezedons[i] ]] <- setAgentState(agent2, 2)##infect!
        }
        
      }
      ##increment each agent 1-day.
      for(i in 1:numAgents) {
        pool[[i]] <- updateAgent(pool[[i]])
      }
      
      states <- sapply(pool, FUN = function(x) {x$biostate})
      distrib <- table(factor(states, levels = 1:7))
      disthistory[day,] <- distrib
      if(plotNetwork) {
        mygplot(cc, socialnetwork, states, main = paste("Day", day))
      }
      
    }
  })
 
}

# Run the application 
shinyApp(ui = ui, server = server)
