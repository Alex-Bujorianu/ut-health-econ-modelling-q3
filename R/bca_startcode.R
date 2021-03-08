<<<<<<< HEAD
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
source("getSingleAttribute.R", echo=T);
source("getMultipleAttributes.R", echo=T);


## Section 2: Data analysis ----

# Load the dataset
load("trial_dataset.RData");

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

func.age <- function() {
  mean_age <- 60; # mean age of patients is 60 years
  age <- rnorm(mean=mean_age)
  return(age)
}

func.tx1.response <- function() {
  prob_tx1_response <- 1;
  response <- ifelse(prob_tx1_response >= runif(1), 1, 0); #1 is response, 0 is non-response
  return(response)
}

func.tx2.response <- function() {
  prob_tx2_response <- 1;
  response <- ifelse(prob_tx2_response >= runif(1), 1, 0); #1 is response, 0 is non-response
  return(response)
}
## Section 3: Supportive functions ----

# Function for determining the event to happen
Tx1.event <- function() {
  
  #Randomly select whether the patient dies with a 10% probability or not
  event <- ifelse(runif(1) < 0.10, 2, 1);
  
  return(event);                                                                                                  # A return value equal to 0 skips the branch and continues to the next activity.
  
} # Function for defining the event during a cycle of Tx1

# Functions for determining the time-to-events
Tx1.time <- function(Tx1.Event) {
  
  return(30);
  
} # Function for defining the time spent on a cycle of Tx1 




## Section 4: Discrete event simulation model ----

# Define the model structure for the current practice, i.e. best standard care (BSC)
bsc.model <- trajectory() %>%
  
  # Initialization
  set_attribute(key="Alive", value=1) %>%                                                                          # define an attribute to check whether the patient is alive
  
  # First-line treatment
  set_attribute(key="Tx1.Event", value=function() Tx1.event()) %>%                                                 # select the event to happen in this treatment cycle          
  branch(option=function() get_attribute(bsc.sim, "Tx1.Event"), continue=c(T, F),
         
         # Event 1: Full cycle
         trajectory() %>%
           set_attribute(key="Tx1.Time", value=function() Tx1.time(get_attribute(bsc.sim, "Tx1.Event"))) %>%       # determine how long the cycle will last
           seize(resource="Tx1", amount=1) %>%                                                                     # occupy a place in first-line treatment
           timeout_from_attribute(key="Tx1.Time") %>%                                                              # stay in first-line treatment for the determined time
           release(resource="Tx1", amount=1) %>%                                                                   # leave first-line treatment
           rollback(amount=6, times=Inf),                                                                          # go back for another cycle (Hint: look at plot trajectory)
         
         # Event 2: Death
         trajectory() %>%
           set_attribute(key="Tx1.Time", value=function() Tx1.time(get_attribute(bsc.sim, "Tx1.Event"))) %>%       # determine how long the cycle will last
           seize(resource="Tx1", amount=1) %>%                                                                     # occupy a place in first-line treatment
           timeout_from_attribute(key="Tx1.Time") %>%                                                              # stay in first-line treatment for the determined time
           release(resource="Tx1", amount=1) %>%                                                                   # leave first-line treatment
           set_attribute(key="Alive", value=0)                                                                     # update that the patient has died
           
         
  ) # branch first-line treatment

# Visualize to check whether the defined model structure is ok
plot(bsc.model);
  




## Section 5: Simulation ----

# Simulation settings
set.seed(5678);       # random number seed for reproducibility
n.patients <- 100;    # number of patients to simulate 
mon.patients <- 2;    # level of monitoring (see add_generator)

# Define simulation for the best standard care (bsc)
bsc.sim <- simmer() %>%
  add_resource(name="Tx1", capacity=Inf, mon=F) %>%
  add_generator(name_prefix="Patient", trajectory=bsc.model, distribution=at(rep(x=0, times=n.patients)), mon=mon.patients)

# Run the BSC simulation
bsc.sim %>% 
  run()

# Get the outcomes for the monitored attributes
bsc.out <- get_mon_attributes(bsc.sim);             # retrieve the monitor object
getSingleAttribute("Alive", bsc.out);               # get patient-level outcomes for the attribute of interest
View(getMultipleAttributes(c("Alive"), bsc.out));   # get outcomes for multiple outcomes at the same time


=======
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
#setwd("C:/path/to/your/AHEM files/");

# Load funtions for extracting monitored attributes
source("getSingleAttribute.R", echo=T);
source("getMultipleAttributes.R", echo=T);


## Section 2: Data analysis ----

# Load the dataset
load("trial_dataset.RData");

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

func.age <- function() {
  mean_age <- 60; # mean age of patients is 60 years
  age <- rnorm(mean=mean_age)
  return(age)
}

## Section 3: Supportive functions ----

# Function for determining the event to happen
Tx1.event <- function() {
  
  #Randomly select whether the patient dies with a 10% probability or not
  event <- ifelse(runif(1) < 0.10, 2, 1);
  
  return(event);                                                                                                  # A return value equal to 0 skips the branch and continues to the next activity.
  
} # Function for defining the event during a cycle of Tx1

# Functions for determining the time-to-events
Tx1.time <- function(Tx1.Event) {
  
  return(30);
  
} # Function for defining the time spent on a cycle of Tx1 




## Section 4: Discrete event simulation model ----

# Define the model structure for the current practice, i.e. best standard care (BSC)
bsc.model <- trajectory() %>%
  
  # Initialization
  set_attribute(key="Alive", value=1) %>%                                                                          # define an attribute to check whether the patient is alive
  
  # First-line treatment
  set_attribute(key="Tx1.Event", value=function() Tx1.event()) %>%                                                 # select the event to happen in this treatment cycle          
  branch(option=function() get_attribute(bsc.sim, "Tx1.Event"), continue=c(T, F),
         
         # Event 1: Full cycle
         trajectory() %>%
           set_attribute(key="Tx1.Time", value=function() Tx1.time(get_attribute(bsc.sim, "Tx1.Event"))) %>%       # determine how long the cycle will last
           seize(resource="Tx1", amount=1) %>%                                                                     # occupy a place in first-line treatment
           timeout_from_attribute(key="Tx1.Time") %>%                                                              # stay in first-line treatment for the determined time
           release(resource="Tx1", amount=1) %>%                                                                   # leave first-line treatment
           rollback(amount=6, times=Inf),                                                                          # go back for another cycle (Hint: look at plot trajectory)
         
         # Event 2: Death
         trajectory() %>%
           set_attribute(key="Tx1.Time", value=function() Tx1.time(get_attribute(bsc.sim, "Tx1.Event"))) %>%       # determine how long the cycle will last
           seize(resource="Tx1", amount=1) %>%                                                                     # occupy a place in first-line treatment
           timeout_from_attribute(key="Tx1.Time") %>%                                                              # stay in first-line treatment for the determined time
           release(resource="Tx1", amount=1) %>%                                                                   # leave first-line treatment
           set_attribute(key="Alive", value=0)                                                                     # update that the patient has died
           
         
  ) # branch first-line treatment

# Visualize to check whether the defined model structure is ok
plot(bsc.model);
  




## Section 5: Simulation ----

# Simulation settings
set.seed(5678);       # random number seed for reproducibility
n.patients <- 100;    # number of patients to simulate 
mon.patients <- 2;    # level of monitoring (see add_generator)

# Define simulation for the best standard care (bsc)
bsc.sim <- simmer() %>%
  add_resource(name="Tx1", capacity=Inf, mon=F) %>%
  add_generator(name_prefix="Patient", trajectory=bsc.model, distribution=at(rep(x=0, times=n.patients)), mon=mon.patients)

# Run the BSC simulation
bsc.sim %>% 
  run()

# Get the outcomes for the monitored attributes
bsc.out <- get_mon_attributes(bsc.sim);             # retrieve the monitor object
getSingleAttribute("Alive", bsc.out);               # get patient-level outcomes for the attribute of interest
View(getMultipleAttributes(c("Alive"), bsc.out));   # get outcomes for multiple outcomes at the same time


>>>>>>> ae77e7d214f6276ddda30db3049dd1e159c16375
