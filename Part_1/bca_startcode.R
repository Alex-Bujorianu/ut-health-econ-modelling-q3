## Advanced Health Economic Modeling
## Assignment Part I: Health Economic Model
##
## Good luck!

## Section 1: Initialization ----

# Clear the workspace
rm(list=ls()); gc();

# Load the required packages, make sure the required packages are installed. For the installation of the packages use line 13 (uncomment shift+Ctrl+C, or remove #)

# install.packages(c("simmer", "simmer.plot", "fitdistrplus")) # Install packages

library(simmer);
library(simmer.plot);
library(fitdistrplus);

# Set the working directory
#setwd("R/");

# Load functions for extracting monitored attributes
source("Part_1/getSingleAttribute.R", echo=T);
source("Part_1/getMultipleAttributes.R", echo=T);

## Section 2: Data analysis ----

# Load the dataset
load("Data/trial_dataset.RData");

#Load the test results from the distribution fitting file and step 1.4
source("Part_2/Dx-distribution-fitting.R", echo=F);
source("Part_1/step-1.4.R", echo=F);

# Define parameters

func.sex <- function() {
  prob_male <- 0.33; #33% of patients are male
  sex <- ifelse(prob_male >= runif(1), 1, 0);
  return(sex)
}

func.condition <- function() {
  prob_poor <- 0.21; #21% of patients are in poor clinical condtion
  condition <- ifelse(prob_poor >= runif(1), 1, 0);
  return(condition)
}

func.age <- function(sex) {
  #Male patients have lognorm
  if (sex==1) {
    mean_age <- mean_male_ages; # mean age of male patients is about 53
    sd_age <- sd_male_ages; # SD of male patients
    location <- log(mean_age^2 / sqrt(sd_age^2 + mean_age^2))
    shape <- sqrt(log(1 + (sd_age^2 / mean_age^2)))
    age <- rlnorm(n=1, meanlog=location,sdlog=shape)
    return(age)
  }
  else {
    age <- rweibull(1, shape = dist_weibull_female$estimate[1], scale = dist_weibull_female$estimate[2])
    return(age)
  }
}

#Baseline function for response. About 45% of patients are responders.
func.tx1.response <- function() {
  prob_tx1_response <- 0.488
  response <- ifelse(prob_tx1_response <= runif(1), 1, 0) #1 is response, 0 is non-response
  return(response)
}
func.tx2.response <- function() {
  prob_tx2_response <- 1;
  response <- ifelse(prob_tx2_response >= runif(1), 1, 0); #1 is response, 0 is non-response
  return(response)
}

#Function to find cycle 1 cost
func.tx1cost<- function(Tx1.cycles, Tx1.time, Tx1.Complications) {
  Tx1_cyclecost <- 504;
  Tx1_daycost <- 8;
  Minor_cost <- 381;
  Major_cost <- 11032;
  if (Tx1.Complications == 1){
    cycle1cost <- Tx1_cyclecost * Tx1.cycles + Tx1_daycost * Tx1.time + 1*Minor_cost;
  }
  else if (Tx1.Complications == 2){
    cycle1cost <- Tx1_cyclecost * Tx1.cycles + Tx1_daycost * Tx1.time + 1*Major_cost;
  }
  else{
    cycle1cost <- Tx1_cyclecost * Tx1.cycles + Tx1_daycost * Tx1.time;
  }
  return(cycle1cost)
}

#Function to find cycle 2 cost
func.tx2cost<- function(Tx2.cycles, Tx2.time, Tx2.Complications) {
  Tx2_cyclecost <- 4450;
  Tx2_daycost <- 14;
  Minor_cost <- 381;
  Major_cost <- 11032;
  if (Tx2.Complications == 1){
    cycle2cost <- Tx2_cyclecost * Tx2.cycles + Tx2_daycost * Tx2.time + 1*Minor_cost;
  }
  else if (Tx2.Complications == 2){
    cycle2cost <- Tx2_cyclecost * Tx2.cycles + Tx2_daycost * Tx2.time + 1*Major_cost;
  }
  else{
    cycle2cost <- Tx2_cyclecost * Tx2.cycles + Tx2_daycost * Tx2.time;
  }
  return(cycle2cost)
}

#Function to find total cost
func.cost<- function(Tx1.Cycles, Tx1.time, Tx1.Complications, Tx2.Cycles, Tx2.time, Tx2.Complications) {
  total_cost <- func.tx1cost(Tx1.Cycles, Tx1.time, Tx1.Complications)+ func.tx2cost(Tx2.Cycles, Tx2.time, Tx2.Complications);
  return(total_cost)
}

#Function to find costs for the exp model taking into account the costs of the novel diagnostic tests
func.exp.costs <- function(Tx1.Cycles, Tx1.time, Tx1.Complications, Tx2.Cycles, Tx2.time, Tx2.Complications) {
  total_cost <- func.tx1cost(Tx1.Cycles, Tx1.time, Tx1.Complications)+ func.tx2cost(Tx2.Cycles, Tx2.time, Tx2.Complications);
  #Note that we do tests at the baseline, so even with 0 cycles we will do at least 1 test
  total_cost <- total_cost + ((Tx1.Cycles+1) * (278+256+194)) + ((Tx2.Cycles+1) * (278+256+194))
  return(total_cost)
}


#Utility of a patient depending on where they are in a cycle
func.utility <- function(position){
  if  (position<5){ #in tx1
    utility =0.55
    if (position==2){ #major complication
      utility	=utility -0.1
    }
    if (position==3){ #minor complication
      utility	=utility -0.05
    }
    if (position==4){#in follow up
      utility =1.1*utility
    }
  }
  else {
    utility=0.5
      if (position==6){ #major complication
        utility	=utility -0.1
      }
      if (position==7){ #minor complication
        utility	=utility -0.05
      }
      if (position==8){#in follow up
        utility =1.1*utility
      }
  }
  return(utility)
}

#This function does NOT return a total. It only returns the qalys during a specific cycle.
#Therefore, it must be called up during each cycle and the result added to the patient's cumulative total.
qaly=0
func.qaly <- function(position, Tx1.time, Tx2.time, followup1.time){
  if (position<4){
    qaly=qaly+(func.utility(position) * Tx1.time /365)
  }
  else if (position==4){
    qaly=qaly+(func.utility(position)*followup1.time/365) #add the utility during followup, which is 10% higher
  }
  else if (position>4 && position<8){
    qaly=qaly+(func.utility(position)*Tx2.time/365)
  }
  #Patient has completed full cycle.
  else {
    qaly=qaly+(func.utility(position)*palliative.time()/365)
  }
  return(qaly)
}

## Section 3: Supportive functions ----

#This is a helper function for Tx1.Event
#I had to refactor
Tx1.continue <- function(response) {
  if (response==1) {
    test_results <- SimDesign::rmvnorm(n = 1, mean = v_means, sigma = m_cov)
    first_test <- test_results[1]
    second_test <- test_results[2]
    third_test <- test_results[3]
    if (first_test > test1_boundary || second_test > test2_boundary || third_test > test3_boundary) {
      return(FALSE) #patient is not responding to treatment, do not overtreat, go to next treatment pathway
    }
    else {
      return(TRUE)
    }
  }
  else {
    test_results_nonresponder <- SimDesign::rmvnorm(n = 1, mean = non_v_means, sigma = non_m_cov)
    first_test_non <- test_results_nonresponder[1]
    second_test_non <- test_results_nonresponder[2]
    third_test_non <- test_results_nonresponder[3]
    if (first_test_non > test1_boundary || second_test_non > test2_boundary || third_test_non > test3_boundary) {
      return(FALSE) #patient is not responding to treatment, do not overtreat, go to next treatment pathway
    }
    else {
      return(TRUE)
    }
  }
}

# Function for determining the event to happen in Tx1 for the exp model
Tx1.Event.alt <- function(cycles, response) {
  #about 10% chance of minor complication, 4% major, 3% death
  minor_comp <- ifelse(runif(1) < 0.1, 1, 0); #10% chance of minor complication, 4% major, 3% death
  major_comp <- ifelse(runif(1) < 0.04, 1, 0);
  death <- ifelse(runif(1) < 0.03, 1, 0);
  #First check to see if the patient is alive, then check for complications
  if (death == 1) {
    return(2)
  }
  else if (cycles > 4){ #if the patient has survived the cycles, they are taken out of the simulation
    return(5)
  }
  #Stop the treatment if the patient is not responding, do not overtreat
  #We do these tests BEFORE the treatment cycle. So complications should not occur. 
  else if (cycles > 0 && Tx1.continue(response)==FALSE) {
    return(5)
  }
  else if (minor_comp == 1) {
    return(4) 
  }
  else if (major_comp == 1) {
    return(3)
  }
  else {
    return(1)
  }                                                                         
}

# Function for determining the event to happen in Tx1
Tx1.Event <- function(cycles) { #careful! This function now takes cycles as a parameter
  minor_comp <- ifelse(runif(1) < 0.1, 1, 0); #10% chance of minor complication, 4% major, 3% death
  major_comp <- ifelse(runif(1) < 0.04, 1, 0);
  death <- ifelse(runif(1) < 0.03, 1, 0);
  if (cycles>4){ #if the patient has survived the cycles, they are taken out of the simulation
    return(5)
  } 
  else if (minor_comp == 1) {
    return(4) 
  }
  else if (major_comp == 1) {
    return(3)
  }
  else if (death == 1) {
    return(2)
  }
  else {
    return(1)
  }                                                                                                  # A return value equal to 0 skips the branch and continues to the next activity.
} 

# Function for determining the event to happen in Tx2
Tx2.Event <- function(cycles) {
  minor_comp <- ifelse(runif(1) < 0.1, 1, 0); #10% chance of minor complication
  major_comp <- ifelse(runif(1) < 0.04, 1, 0);
  death <- ifelse(runif(1) < 0.03, 1, 0);
  if (cycles>4){ #if the patient has survived the cycles, they are taken out of the simulation
    return(5)
  } 
  else if (minor_comp == 1) {
    return(4) #function returns and stops execution
  }
  else if (major_comp == 1) {
    return(3)
  }
  else if (death == 1) {
    return(2)
  }
  else {
    return(1)
  }                                                                                                  # A return value equal to 0 skips the branch and continues to the next activity.
}

#Helper function for Tx2.Event.alt
Tx2.continue <- function(Tx1.cycles, response) {
  #Responders
  if (response==1) {
    decision <- switch(Tx1.cycles, #This switch evaluates based on the ORDER. If Tx1.cycles = 4 it evaluates the 4th statement.
                       'first-cycle'=ifelse(0.52 >= runif(1), TRUE, FALSE), #if responding, continue
                       'second-cycle'=ifelse(0.54 >= runif(1), TRUE, FALSE),
                       'third-cycle'=ifelse(0.56 >= runif(1), TRUE, FALSE),
                       'fourth-cycle'=ifelse(0.58 >= runif(1), TRUE, FALSE),
                       'fifth-cycle'=ifelse(0.6 >= runif(1), TRUE, FALSE)
    )
    return(decision)
  }
  #Non responders
  else {
    decision <- switch(Tx1.cycles, #This switch evaluates based on the ORDER. If Tx1.cycles = 4 it evaluates the 4th statement.
                       'first-cycle'=ifelse(0.9 >= runif(1), TRUE, FALSE), #if responding, continue
                       'second-cycle'=ifelse(0.78 >= runif(1), TRUE, FALSE),
                       'third-cycle'=ifelse(0.66 >= runif(1), TRUE, FALSE),
                       'fourth-cycle'=ifelse(0.54 >= runif(1), TRUE, FALSE),
                       'fifth-cycle'=ifelse(0.42 >= runif(1), TRUE, FALSE)
    )
    return(decision)
  }
}
#Event function for Tx2 in exp sim
Tx2.Event.alt <- function(cycles, Tx1.cycles, response) {
  minor_comp <- ifelse(runif(1) < 0.1, 1, 0); #10% chance of minor complication
  major_comp <- ifelse(runif(1) < 0.04, 1, 0);
  death <- ifelse(runif(1) < 0.03, 1, 0);
  #Events are mutually exclusive. Test for the events with bigger consequences first.
  #First check for death
  if (death == 1) {
    return(2)
  }
  else if (cycles>4){ #if the patient has survived the cycles, they are taken out of the simulation
    return(5)
  }
  #CHECK IF THE PATIENT HAS COMPLETED THE TREATMENT CYCLE FIRST
  #Otherwise, this simulation will run in an infinite loop.
  #Check for overtreatment. If Tx2.continue gives us FALSE, put the patient out of the treatment cycle
  else if (cycles > 0 && Tx2.continue(Tx1.cycles, response) == FALSE) {
    return(5)
  } 
  else if (major_comp == 1) {
    return(3)
  }
  else if (minor_comp == 1) {
    return(4) 
  }
  else {
    return(1)
  }          
}


# Function for determining the event to happen in Follow up 1
followup1.event <- function() {
  #1: patient lives; 2: patient dies
  prob_death <- ifelse(runif(1) < beta_death_followup, 2, 1);
  return(prob_death)
}

# Functions for determining the time-to-events

# Function for defining the time spent on a cycle of Tx1 
Tx1.time <- function(Tx1.Event) {
  #If patient has no issues or only minor complications, the full cycle duration occurs
  if (Tx1.Event == 1 || Tx1.Event == 4) {
  return(30);
  }
  else if (Tx1.Event == 2) {
    time <- runif(n=1, min=2,max=28) #incorporating uncertainties of time to death
    return(time)
  }
  else if (Tx1.Event == 3) {
    time <- rgamma(n=1, shape=47.61169, rate=8.09397) #incorporating uncertainties of time to major complicatons
    return(time)
  }
}

#Tx2 is identical to tx1
Tx2.time <- function(Tx2.Event) {
  #If patient has no issues or only minor complications, the full cycle duration occurs
  if (Tx2.Event == 1 || Tx2.Event == 4) {
    return(30);
  }
  else if (Tx2.Event == 2) {
    return(15); #deaths tend to occur on day 15
  }
  else if (Tx2.Event == 3) {
    return(6);
  }
}

# Function for defining the time spent in Follow up 1 
followup1.time <- function(followup1.event) {
  if (followup1.event == 1) {
    return(63); #If patient survives this period lasts 63 days
  }
  else {
    return(42); #deaths tend to occur on day 42
  }
}

# Function for defining the time spent in palliative care 
palliative.time <- function() {
  return(100);
}


## Section 4: Discrete event simulation model ----

# Define the model structure for the current practice, i.e. best standard care (BSC)
bsc.model <- trajectory()%>%
  
  # Initialization: do not forget to initialise these cycle attributes or the patients will all die.
  set_attribute(key="Tx1.Cycles", value = 0) %>%
  set_attribute(key="Tx2.Cycles", value = 0) %>%
  set_attribute(key="Tx1.Complications", value=0) %>% #0 for no complications, 1 for minor, 2 for major
  set_attribute(key="Tx2.Complications", value=0) %>%
  set_attribute(key="Tx1.Time", value=0) %>%
  set_attribute(key="Tx2.Time", value=0) %>%
  set_attribute(key="position", value=0) %>%
  set_attribute(key="condition", mod="+", value=function() func.condition()) %>%
  set_attribute(key="response", mod="+", value=function() func.tx1.response()) %>%
  set_attribute(key="sex", value=function() func.sex()) %>%
  set_attribute(key="age", value=function() func.age(get_attribute(bsc.sim, "sex"))) %>%
  set_attribute(key="qalys", value=0) %>%
  set_attribute(key="cost", value=0) %>% #keep track of each patient's cost
  set_attribute(key="Alive", value=1) %>%                                                                          # define an attribute to check whether the patient is alive
  
  # First-line treatment
  set_attribute(key="Tx1.Event", value=function() Tx1.Event(cycles = get_attribute(bsc.sim, "Tx1.Cycles"))) %>%                                                 # select the event to happen in this treatment cycle          
  branch(option=function() get_attribute(bsc.sim, "Tx1.Event"), continue=c(T, F, F, T, T), #don't forget 5th is false
         
         # Event 1: Full cycle
         trajectory() %>%
           set_attribute(key="Tx1.Time", value=function() Tx1.time(get_attribute(bsc.sim, "Tx1.Event"))) %>%       # determine how long the cycle will last
           seize(resource="Tx1", amount=1) %>%                                                                     # occupy a place in first-line treatment
           timeout_from_attribute(key="Tx1.Time") %>%                                                              # stay in first-line treatment for the determined time
           release(resource="Tx1", amount=1) %>%                                                                   # leave first-line treatment
           set_attribute(keys = "Tx1.Cycles", mod = "+", value = 1)%>%
           set_attribute(key="position", value=1) %>%
           set_attribute(key="qalys",  mod = "+", value=function() func.qaly(get_attribute(bsc.sim, "position"), 
                                                                  get_attribute(bsc.sim, "Tx1.Time"), 0, 0)) %>%
           rollback(amount=9, times=Inf),                                                                                             # go back for another cycle (Hint: look at plot trajectory)
         
         # Event 2: Death
         trajectory() %>%
           set_attribute(key="Tx1.Time", value=function() Tx1.time(get_attribute(bsc.sim, "Tx1.Event"))) %>%       # determine how long the cycle will last
           seize(resource="Tx1", amount=1) %>%                                                                     # occupy a place in first-line treatment
           timeout_from_attribute(key="Tx1.Time") %>%                                                              # stay in first-line treatment for the determined time
           release(resource="Tx1", amount=1) %>%
           log_("Patient has died during first line treatment") %>% # leave first-line treatment
           set_attribute(keys = "Tx1.Cycles", mod = "+", value = 1) %>%
           set_attribute(key="qalys",  mod = "+", value=function() func.qaly(get_attribute(bsc.sim, "position"), 
                                                                             get_attribute(bsc.sim, "Tx1.Time"), 0, 0)) %>%
           #Calculate the costs when a patient dies
           set_attribute(key="cost", value=function() func.cost(
             get_attribute(bsc.sim, "Tx1.Cycles"),
             get_attribute(bsc.sim, "Tx1.Time"),
             get_attribute(bsc.sim, "Tx1.Complications"),
             get_attribute(bsc.sim, "Tx2.Cycles"),
             get_attribute(bsc.sim, "Tx2.Time"),
             get_attribute(bsc.sim, "Tx2.Complications")
           )) %>%
           set_attribute(key="Alive", value=0),                                                                     # update that the patient has died
           
         # Event 3: Major Complications
         trajectory() %>%
           set_attribute(key="Tx1.Time", value=function() Tx1.time(get_attribute(bsc.sim, "Tx1.Event"))) %>%       # determine how long the cycle will last
           seize(resource="Tx1", amount=1) %>%                                                                     # occupy a place in first-line treatment
           timeout_from_attribute(key="Tx1.Time") %>%                                                              # stay in first-line treatment for the determined time
           release(resource="Tx1", amount=1) %>%                                                                   # leave first-line treatment
           set_attribute(keys = "Tx1.Cycles", mod = "+", value = 1)%>%
           set_attribute(keys="Tx1.Complications", value=2) %>% #Now we know the patient had a major comp in tx1
           set_attribute(key="position", value=2) %>%
           set_attribute(key="qalys",  mod = "+", value=function() func.qaly(get_attribute(bsc.sim, "position"), 
                                                                             get_attribute(bsc.sim, "Tx1.Time"), 0, 0)) %>%
           rollback(amount=10, times=Inf),                                                                          # go back for another cycle (Hint: look at plot trajectory)
         
         # Event 4: Minor Complications
         trajectory() %>%
           set_attribute(key="Tx1.Time", value=function() Tx1.time(get_attribute(bsc.sim, "Tx1.Event"))) %>%       # determine how long the cycle will last
           seize(resource="Tx1", amount=1) %>%                                                                     # occupy a place in first-line treatment
           timeout_from_attribute(key="Tx1.Time") %>%                                                              # stay in first-line treatment for the determined time
           release(resource="Tx1", amount=1) %>%                                                                   # leave first-line treatment
           set_attribute(keys = "Tx1.Cycles", mod = "+", value = 1)%>%
           set_attribute(keys="Tx1.Complications", value=1) %>% #1 for minor complications
           set_attribute(key="position", value=3) %>%
           set_attribute(key="qalys",  mod = "+", value=function() func.qaly(get_attribute(bsc.sim, "position"), 
                                                                             get_attribute(bsc.sim, "Tx1.Time"), 0, 0)) %>%
           rollback(amount=10, times=Inf), 
         #Fifth trajectory: the patient has survived all treatment cycles and is out of first line treatment
         trajectory()%>%
           timeout(10)
  )%>% 
  
  #First follow up period
  set_attribute(key="followup1.event", value=function() followup1.event()) %>%
  branch(option=function() get_attribute(bsc.sim, "followup1.event"), continue=c(T, F),
  
  trajectory() %>% #First option: they survive the follow up
    set_attribute(key="Fu1.Time", value=function() followup1.time(get_attribute(bsc.sim, "followup1.event"))) %>%       
    seize(resource="Fu1", amount=1) %>%                                                                     
    timeout_from_attribute(key="Fu1.Time") %>%                                                              
    release(resource="Fu1", amount=1) %>%                                                                   
    # rollback(amount=6, times=Inf),
    set_attribute(key="position", value=4) %>%
    set_attribute(key="qalys",  mod = "+", value=function() func.qaly(get_attribute(bsc.sim, "position"), 
                                                                      get_attribute(bsc.sim, "Tx1.Time"), 0, 10)) %>%
    timeout(10),

  trajectory() %>% #Second option: they die during the follow up
    set_attribute(key="Fu1.Time", value=function() followup1.time(get_attribute(bsc.sim, "followup1.event"))) %>%       
    seize(resource="Fu1", amount=1) %>%                                                                     
    timeout_from_attribute(key="Fu1.Time") %>%                                                              
    release(resource="Fu1", amount=1) %>%
    log_("Patient has died during follow up") %>%
    timeout(10)%>%
    set_attribute(key="qalys",  mod = "+", value=function() func.qaly(get_attribute(bsc.sim, "position"), 
                                                                      get_attribute(bsc.sim, "Tx1.Time"), 0, 10)) %>%
    set_attribute(key="cost", value=function() func.cost(
      get_attribute(bsc.sim, "Tx1.Cycles"),
      get_attribute(bsc.sim, "Tx1.Time"),
      get_attribute(bsc.sim, "Tx1.Complications"),
      get_attribute(bsc.sim, "Tx2.Cycles"),
      get_attribute(bsc.sim, "Tx2.Time"),
      get_attribute(bsc.sim, "Tx2.Complications")
    )) %>%
    set_attribute(key="Alive", value=0)
  ) %>%
         #Second line treatment
           set_attribute(key="Tx2.Event", value=function() Tx2.Event(cycles = get_attribute(bsc.sim, "Tx2.Cycles"))) %>%                                                 # select the event to happen in this treatment cycle          
           branch(option=function() get_attribute(bsc.sim, "Tx2.Event"), continue=c(T, F, F, T, T), #the 5th has to be false
                  # Event 1: Full cycle
                  trajectory() %>%
                    set_attribute(key="Tx2.Time", value=function() Tx2.time(get_attribute(bsc.sim, "Tx2.Event"))) %>%       # determine how long the cycle will last
                    seize(resource="Tx2", amount=1) %>%                                                                     # occupy a place in first-line treatment
                    timeout_from_attribute(key="Tx2.Time") %>%                                                              # stay in first-line treatment for the determined time
                    release(resource="Tx2", amount=1) %>%                                                                   # leave first-line treatment
                    set_attribute(keys = "Tx2.Cycles", mod = "+", value = 1)%>%
                    set_attribute(key="position", value=5) %>%
                    set_attribute(key="qalys",  mod = "+", value=function() func.qaly(get_attribute(bsc.sim, "position"), 
                                                                                      get_attribute(bsc.sim, "Tx1.Time"), 
                                                                                      get_attribute(bsc.sim, "Tx2.Time"), 10)) %>%
                    rollback(amount=9, times=Inf),                                                                       # go back for another cycle (Hint: look at plot trajectory)
                 
                  # Event 2: Death
                  trajectory() %>%
                    set_attribute(key="Tx2.Time", value=function() Tx2.time(get_attribute(bsc.sim, "Tx2.Event"))) %>%       # determine how long the cycle will last
                    seize(resource="Tx2", amount=1) %>%                                                                     # occupy a place in first-line treatment
                    timeout_from_attribute(key="Tx2.Time") %>%                                                              # stay in first-line treatment for the determined time
                    release(resource="Tx2", amount=1) %>%
                    log_("Patient has died during second line treatment") %>%
                    set_attribute(keys = "Tx2.Cycles", mod = "+", value = 1)%>% # leave second-line treatment
                    set_attribute(key="qalys",  mod = "+", value=function() func.qaly(get_attribute(bsc.sim, "position"), 
                                                                                      get_attribute(bsc.sim, "Tx1.Time"), 
                                                                                      get_attribute(bsc.sim, "Tx2.Time"), 10)) %>%
                    set_attribute(key="cost", value=function() func.cost(
                      get_attribute(bsc.sim, "Tx1.Cycles"),
                      get_attribute(bsc.sim, "Tx1.Time"),
                      get_attribute(bsc.sim, "Tx1.Complications"),
                      get_attribute(bsc.sim, "Tx2.Cycles"),
                      get_attribute(bsc.sim, "Tx2.Time"),
                      get_attribute(bsc.sim, "Tx2.Complications")
                    )) %>%
                    set_attribute(key="Alive", value=0),                                                                   # update that the patient has died
                  
                  # Event 3: Major Complications
                  trajectory() %>%
                    set_attribute(key="Tx2.Time", value=function() Tx1.time(get_attribute(bsc.sim, "Tx2.Event"))) %>%       # determine how long the cycle will last
                    seize(resource="Tx2", amount=1) %>%                                                                     # occupy a place in first-line treatment
                    timeout_from_attribute(key="Tx2.Time") %>%                                                              # stay in first-line treatment for the determined time
                    release(resource="Tx2", amount=1) %>%                                                                   # leave first-line treatment
                    set_attribute(keys = "Tx2.Cycles", mod = "+", value = 1)%>%
                    set_attribute(keys="Tx2.Complications", value=2) %>%
                    set_attribute(key="position", value=6) %>%
                    set_attribute(key="qalys",  mod = "+", value=function() func.qaly(get_attribute(bsc.sim, "position"), 
                                                                                      get_attribute(bsc.sim, "Tx1.Time"), 
                                                                                      get_attribute(bsc.sim, "Tx2.Time"), 10)) %>%
                    rollback(amount=10, times=Inf),                                                                          # go back for another cycle (Hint: look at plot trajectory)
                  
                  # Event 4: Minor Complications
                  trajectory() %>%
                    set_attribute(key="Tx2.Time", value=function() Tx1.time(get_attribute(bsc.sim, "Tx2.Event"))) %>%       # determine how long the cycle will last
                    seize(resource="Tx2", amount=1) %>%                                                                     # occupy a place in first-line treatment
                    timeout_from_attribute(key="Tx2.Time") %>%                                                              # stay in first-line treatment for the determined time
                    release(resource="Tx2", amount=1) %>%                                                                   # leave first-line treatment
                    set_attribute(keys = "Tx2.Cycles", mod = "+", value = 1)%>%
                    set_attribute(keys="Tx2.Complications", value = 1) %>%
                    set_attribute(key="position", value=7) %>%
                    set_attribute(key="qalys",  mod = "+", value=function() func.qaly(get_attribute(bsc.sim, "position"), 
                                                                                      get_attribute(bsc.sim, "Tx1.Time"), 
                                                                                      get_attribute(bsc.sim, "Tx2.Time"), 10)) %>%
                    rollback(amount=10, times=Inf),
                  #Fifth trajectory: the patient has survived all treatment cycles
                  trajectory()%>%
                    timeout(10)
         
  ) %>%
  set_attribute(key="palliative.time", value=function() palliative.time()) %>%
  branch(option=function() 1, continue=c(F),
  trajectory() %>%
    seize(resource="Fu2", amount=1) %>%
    timeout_from_attribute(key="palliative.time") %>%
    log_("Patient has survived the all treatment cycles") %>%
    set_attribute(key="position", value=8) %>%
    set_attribute(key="qalys",  mod = "+", value=function() func.qaly(get_attribute(bsc.sim, "position"), 
                                                                      get_attribute(bsc.sim, "Tx1.Time"), 
                                                                      get_attribute(bsc.sim, "Tx2.Time"), 10)) %>%
    set_attribute(key="cost", value=function() func.cost(
      get_attribute(bsc.sim, "Tx1.Cycles"),
      get_attribute(bsc.sim, "Tx1.Time"),
      get_attribute(bsc.sim, "Tx1.Complications"),
      get_attribute(bsc.sim, "Tx2.Cycles"),
      get_attribute(bsc.sim, "Tx2.Time"),
      get_attribute(bsc.sim, "Tx2.Complications")
    )) %>%
    release(resource="Fu2", amount=1)
  )


## EXP Model

exp.model <- trajectory()%>%
  
  # Initialization: do not forget to initialise these cycle attributes or the patients will all die.
  set_attribute(key="Tx1.Cycles", value = 0) %>%
  set_attribute(key="Tx2.Cycles", value = 0) %>%
  #Do not forget to initialise the times!
  set_attribute(key="Tx1.Time", value=0) %>%
  set_attribute(key="Tx2.Time", value=0) %>%
  set_attribute(key="Tx1.Complications", value=0) %>% #0 for no complications, 1 for minor, 2 for major
  set_attribute(key="Tx2.Complications", value=0) %>%
  set_attribute(key="position", value=0) %>%
  set_attribute(key="condition", mod="+", value=function() func.condition()) %>%
  set_attribute(key="response", mod="+", value=function() func.tx1.response()) %>%
  set_attribute(key="sex", value=function() func.sex()) %>%
  set_attribute(key="age", value=function() func.age(get_attribute(exp.sim, "sex"))) %>%
  set_attribute(key="qalys", value=0) %>%
  set_attribute(key="Alive", value=1) %>%                                                                          # define an attribute to check whether the patient is alive
  set_attribute(key="cost", value=0) %>%
  
  # First-line treatment
  set_attribute(key="Tx1.Event", value=function() Tx1.Event.alt(cycles = get_attribute(exp.sim, "Tx1.Cycles"),
                                                                response = get_attribute(exp.sim, "response")) ) %>%                                                 # select the event to happen in this treatment cycle          
  branch(option=function() get_attribute(exp.sim, "Tx1.Event"), continue=c(T, F, F, T, T),
         
         # Event 1: Full cycle
         trajectory() %>%
           set_attribute(key="Tx1.Time", value=function() Tx1.time(get_attribute(exp.sim, "Tx1.Event"))) %>%       # determine how long the cycle will last
           seize(resource="Tx1", amount=1) %>%                                                                     # occupy a place in first-line treatment
           timeout_from_attribute(key="Tx1.Time") %>%                                                              # stay in first-line treatment for the determined time
           release(resource="Tx1", amount=1) %>%                                                                   # leave first-line treatment
           set_attribute(keys = "Tx1.Cycles", mod = "+", value = 1)%>%
           set_attribute(key="position", value=1) %>%
           set_attribute(key="qalys",  mod = "+", value=function() func.qaly(get_attribute(exp.sim, "position"), 
                                                                             get_attribute(exp.sim, "Tx1.Time"), 0, 0)) %>%
           rollback(amount=9, times=Inf),                                                                                             # go back for another cycle (Hint: look at plot trajectory)
         
         # Event 2: Death
         trajectory() %>%
           set_attribute(key="Tx1.Time", value=function() Tx1.time(get_attribute(exp.sim, "Tx1.Event"))) %>%       # determine how long the cycle will last
           seize(resource="Tx1", amount=1) %>%                                                                     # occupy a place in first-line treatment
           timeout_from_attribute(key="Tx1.Time") %>%                                                              # stay in first-line treatment for the determined time
           release(resource="Tx1", amount=1) %>%
           log_("Patient has died during first line treatment") %>% # leave first-line treatment
           set_attribute(keys = "Tx1.Cycles", mod = "+", value = 1) %>%
           set_attribute(key="qalys",  mod = "+", value=function() func.qaly(get_attribute(exp.sim, "position"), 
                                                                             get_attribute(exp.sim, "Tx1.Time"), 0, 0)) %>%
           set_attribute(key="Alive", value=0)%>%                                                                  # update that the patient has died
           set_attribute(key="cost", mod="+", value=function()func.exp.costs(get_attribute(exp.sim, "Tx1.Cycles"),
                                                                             get_attribute(exp.sim, "Tx1.Time"),
                                                                             get_attribute(exp.sim, "Tx1.Complications"),
                                                                             get_attribute(exp.sim, "Tx2.Cycles"),
                                                                             get_attribute(exp.sim, "Tx2.Time"),
                                                                             get_attribute(exp.sim, "Tx2.Complications"))),
         
         # Event 3: Major Complications
         trajectory() %>%
           set_attribute(key="Tx1.Time", value=function() Tx1.time(get_attribute(exp.sim, "Tx1.Event"))) %>%       # determine how long the cycle will last
           seize(resource="Tx1", amount=1) %>%                                                                     # occupy a place in first-line treatment
           timeout_from_attribute(key="Tx1.Time") %>%                                                              # stay in first-line treatment for the determined time
           release(resource="Tx1", amount=1) %>%                                                                   # leave first-line treatment
           set_attribute(keys = "Tx1.Cycles", mod = "+", value = 1)%>%
           set_attribute(keys="Tx1.Complications", value=2) %>% #Now we know the patient had a major comp in tx1
           set_attribute(key="position", value=2) %>%
           set_attribute(key="qalys",  mod = "+", value=function() func.qaly(get_attribute(exp.sim, "position"), 
                                                                             get_attribute(exp.sim, "Tx1.Time"), 0, 0)) %>%
           rollback(amount=10, times=Inf),                                                                          # go back for another cycle (Hint: look at plot trajectory)
         
         # Event 4: Minor Complications
         trajectory() %>%
           set_attribute(key="Tx1.Time", value=function() Tx1.time(get_attribute(exp.sim, "Tx1.Event"))) %>%       # determine how long the cycle will last
           seize(resource="Tx1", amount=1) %>%                                                                     # occupy a place in first-line treatment
           timeout_from_attribute(key="Tx1.Time") %>%                                                              # stay in first-line treatment for the determined time
           release(resource="Tx1", amount=1) %>%                                                                   # leave first-line treatment
           set_attribute(keys = "Tx1.Cycles", mod = "+", value = 1)%>%
           set_attribute(keys="Tx1.Complications", value=1) %>% #1 for minor complications
           set_attribute(key="position", value=3) %>%
           set_attribute(key="qalys",  mod = "+", value=function() func.qaly(get_attribute(exp.sim, "position"), 
                                                                             get_attribute(exp.sim, "Tx1.Time"), 0, 0)) %>%
           rollback(amount=10, times=Inf), 
         #Fifth trajectory: the patient has survived all treatment cycles and is out of first line treatment
         trajectory()%>%
           timeout(10)
  )%>% 
  
  #First follow up period
  set_attribute(key="followup1.event", value=function() followup1.event()) %>%
  branch(option=function() get_attribute(exp.sim, "followup1.event"), continue=c(T, F),
         
         trajectory() %>% #First option: they survive the follow up
           set_attribute(key="Fu1.Time", value=function() followup1.time(get_attribute(exp.sim, "followup1.event"))) %>%       
           seize(resource="Fu1", amount=1) %>%                                                                     
           timeout_from_attribute(key="Fu1.Time") %>%                                                              
           release(resource="Fu1", amount=1) %>%                                                                   
           # rollback(amount=6, times=Inf),
           set_attribute(key="position", value=4) %>%
           set_attribute(key="qalys",  mod = "+", value=function() func.qaly(get_attribute(exp.sim, "position"), 
                                                                             get_attribute(exp.sim, "Tx1.Time"), 0, 10)) %>%
           timeout(10),
         
         trajectory() %>% #Second option: they die during the follow up
           set_attribute(key="Fu1.Time", value=function() followup1.time(get_attribute(exp.sim, "followup1.event"))) %>%       
           seize(resource="Fu1", amount=1) %>%                                                                     
           timeout_from_attribute(key="Fu1.Time") %>%                                                              
           release(resource="Fu1", amount=1) %>%
           log_("Patient has died during follow up") %>%
           timeout(10)%>%
           set_attribute(key="qalys",  mod = "+", value=function() func.qaly(get_attribute(exp.sim, "position"), 
                                                                             get_attribute(exp.sim, "Tx1.Time"), 0, 10)) %>%
           set_attribute(key="Alive", value=0)%>%                                                                  # update that the patient has died
           set_attribute(key="cost", mod="+", value=function()func.exp.costs(get_attribute(exp.sim, "Tx1.Cycles"),
                                                                             get_attribute(exp.sim, "Tx1.Time"),
                                                                             get_attribute(exp.sim, "Tx1.Complications"),
                                                                             get_attribute(exp.sim, "Tx2.Cycles"),
                                                                             get_attribute(exp.sim, "Tx2.Time"),
                                                                             get_attribute(exp.sim, "Tx2.Complications")))
         
  ) %>%
  #Second line treatment
  set_attribute(key="Tx2.Event", value=function() Tx2.Event.alt(cycles = get_attribute(exp.sim, "Tx2.Cycles"),
                                                                Tx1.cycles = get_attribute(exp.sim, "Tx1.Cycles"),
                                                                response = get_attribute(exp.sim, "response")))%>%                                                 # select the event to happen in this treatment cycle          
  branch(option=function() get_attribute(exp.sim, "Tx2.Event"), continue=c(T, F, F, T, T), 
         # Event 1: Full cycle
         trajectory() %>%
           set_attribute(key="Tx2.Time", value=function() Tx2.time(get_attribute(exp.sim, "Tx2.Event"))) %>%       # determine how long the cycle will last
           seize(resource="Tx2", amount=1) %>%                                                                     # occupy a place in first-line treatment
           timeout_from_attribute(key="Tx2.Time") %>%                                                              # stay in first-line treatment for the determined time
           release(resource="Tx2", amount=1) %>%                                                                   # leave first-line treatment
           set_attribute(keys = "Tx2.Cycles", mod = "+", value = 1)%>%
           set_attribute(key="position", value=5) %>%
           set_attribute(key="qalys",  mod = "+", value=function() func.qaly(get_attribute(exp.sim, "position"), 
                                                                             get_attribute(exp.sim, "Tx1.Time"), 
                                                                             get_attribute(exp.sim, "Tx2.Time"), 10)) %>%
           rollback(amount=9, times=Inf),                                                                       # go back for another cycle (Hint: look at plot trajectory)
         
         # Event 2: Death
         trajectory() %>%
           set_attribute(key="Tx2.Time", value=function() Tx2.time(get_attribute(exp.sim, "Tx2.Event"))) %>%       # determine how long the cycle will last
           seize(resource="Tx2", amount=1) %>%                                                                     # occupy a place in first-line treatment
           timeout_from_attribute(key="Tx2.Time") %>%                                                              # stay in first-line treatment for the determined time
           release(resource="Tx2", amount=1) %>%
           log_("Patient has died during second line treatment") %>%
           set_attribute(keys = "Tx2.Cycles", mod = "+", value = 1)%>% # leave second-line treatment
           set_attribute(key="qalys",  mod = "+", value=function() func.qaly(get_attribute(exp.sim, "position"), 
                                                                             get_attribute(exp.sim, "Tx1.Time"), 
                                                                             get_attribute(exp.sim, "Tx2.Time"), 10)) %>%
           set_attribute(key="Alive", value=0)%>%                                                                  # update that the patient has died
           set_attribute(key="cost", mod="+", value=function()func.exp.costs(get_attribute(exp.sim, "Tx1.Cycles"),
                                                                             get_attribute(exp.sim, "Tx1.Time"),
                                                                             get_attribute(exp.sim, "Tx1.Complications"),
                                                                             get_attribute(exp.sim, "Tx2.Cycles"),
                                                                             get_attribute(exp.sim, "Tx2.Time"),
                                                                             get_attribute(exp.sim, "Tx2.Complications"))),
         
         # Event 3: Major Complications
         trajectory() %>%
           set_attribute(key="Tx2.Time", value=function() Tx1.time(get_attribute(exp.sim, "Tx2.Event"))) %>%       # determine how long the cycle will last
           seize(resource="Tx2", amount=1) %>%                                                                     # occupy a place in first-line treatment
           timeout_from_attribute(key="Tx2.Time") %>%                                                              # stay in first-line treatment for the determined time
           release(resource="Tx2", amount=1) %>%                                                                   # leave first-line treatment
           set_attribute(keys = "Tx2.Cycles", mod = "+", value = 1)%>%
           set_attribute(keys="Tx2.Complications", value=2) %>%
           set_attribute(key="position", value=6) %>%
           set_attribute(key="qalys",  mod = "+", value=function() func.qaly(get_attribute(exp.sim, "position"), 
                                                                             get_attribute(exp.sim, "Tx1.Time"), 
                                                                             get_attribute(exp.sim, "Tx2.Time"), 10)) %>%
           rollback(amount=10, times=Inf),                                                                          # go back for another cycle (Hint: look at plot trajectory)
         
         # Event 4: Minor Complications
         trajectory() %>%
           set_attribute(key="Tx2.Time", value=function() Tx1.time(get_attribute(exp.sim, "Tx2.Event"))) %>%       # determine how long the cycle will last
           seize(resource="Tx2", amount=1) %>%                                                                     # occupy a place in first-line treatment
           timeout_from_attribute(key="Tx2.Time") %>%                                                              # stay in first-line treatment for the determined time
           release(resource="Tx2", amount=1) %>%                                                                   # leave first-line treatment
           set_attribute(keys = "Tx2.Cycles", mod = "+", value = 1)%>%
           set_attribute(keys="Tx2.Complications", value = 1) %>%
           set_attribute(key="position", value=7) %>%
           set_attribute(key="qalys",  mod = "+", value=function() func.qaly(get_attribute(exp.sim, "position"), 
                                                                             get_attribute(exp.sim, "Tx1.Time"), 
                                                                             get_attribute(exp.sim, "Tx2.Time"), 10)) %>%
           rollback(amount=10, times=Inf),
         #Fifth trajectory: the patient has survived all treatment cycles
         trajectory()%>%
           timeout(10)
         
  ) %>%
  set_attribute(key="palliative.time", value=function() palliative.time()) %>%
  branch(option=function() 1, continue=c(F),
         trajectory() %>%
           seize(resource="Fu2", amount=1) %>%
           timeout_from_attribute(key="palliative.time") %>%
           log_("Patient has survived the all treatment cycles") %>%
           set_attribute(key="position", value=8) %>%
           set_attribute(key="qalys",  mod = "+", value=function() func.qaly(get_attribute(exp.sim, "position"), 
                                                                             get_attribute(exp.sim, "Tx1.Time"), 
                                                                             get_attribute(exp.sim, "Tx2.Time"), 10)) %>%
           release(resource="Fu2", amount=1) %>%
           set_attribute(key="cost", mod="+", value=function()func.exp.costs(get_attribute(exp.sim, "Tx1.Cycles"),
                                                                             get_attribute(exp.sim, "Tx1.Time"),
                                                                             get_attribute(exp.sim, "Tx1.Complications"),
                                                                             get_attribute(exp.sim, "Tx2.Cycles"),
                                                                             get_attribute(exp.sim, "Tx2.Time"),
                                                                             get_attribute(exp.sim, "Tx2.Complications")))
         
  )


# Visualize to check whether the defined model structure is ok
plot(bsc.model)
plot(exp.model)

## Section 5: Simulation ----

# Simulation settings
set.seed(5678);       # random number seed for reproducibility
n.patients <- 10000;    # number of patients to simulate 
mon.patients <- 5;    # level of monitoring (see add_generator)

# Define simulation for the best standard care (bsc)
bsc.sim <- simmer() %>%
  add_resource(name="Tx1", capacity=Inf, mon=F) %>%
  add_resource(name="Tx2", capacity=Inf, mon=F) %>%
  add_resource(name="Fu1", capacity=Inf, mon=F) %>%
  add_resource(name="Fu2", capacity=Inf, mon=F) %>%
  add_generator(name_prefix="Patient", trajectory=bsc.model, distribution=at(rep(x=0, times=n.patients)), mon=mon.patients)

# Define simulation for the experimental care (exp)
exp.sim <- simmer() %>%
  add_resource(name="Tx1", capacity=Inf, mon=F) %>%
  add_resource(name="Tx2", capacity=Inf, mon=F) %>%
  add_resource(name="Fu1", capacity=Inf, mon=F) %>%
  add_resource(name="Fu2", capacity=Inf, mon=F) %>%
  add_generator(name_prefix="Patient", trajectory=exp.model, distribution=at(rep(x=0, times=n.patients)), mon=mon.patients)

# Run the BSC simulation
bsc.sim %>% 
  run()

# Run the EXP simulation
exp.sim %>% 
  run()

# Get the outcomes for the monitored attributes for bsc
bsc.out <- get_mon_attributes(bsc.sim);             # retrieve the monitor object
bsc_attributes <- getMultipleAttributes(c("Alive", "Tx1.Event", "Tx1.Complications", "Tx2.Event", "cost", "qalys"), bsc.out)
#Save data to file
write.csv(bsc_attributes, "Data/bsc-results-10000.csv")

# Get the outcomes for the monitored attributes for exp
exp.out <- get_mon_attributes(exp.sim);             # retrieve the monitor object
exp_attributes <- getMultipleAttributes(c("Alive", "Tx1.Event", "Tx1.Complications", "Tx2.Event", "cost", "qalys"), exp.out)
#Save data to file
write.csv(exp_attributes, "Data/exp-results-10000.csv")
