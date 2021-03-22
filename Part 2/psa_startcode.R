## Section 1: Initialization ----
# Clear the workspace
rm(list=ls()); gc();

# Load the required packages
library(parallel);
library(doSNOW);

# Set the working directory
setwd("C:/path/to/your/AHEM files/");

# Load funtions for extracting monitored attributes
source("getSingleAttribute.R", echo=T);
source("getMultipleAttributes.R", echo=T);



## Section 2: Simulation function ----

runPSA <- function(n.patients, n.runs, free.cores=1, seed=1234) {
  
  # Input information:
  # - n.patients    the number of patients to be simulated in each simulation
  # - n.runs        the number of probabilistic sensitivity analysis runs to be performed
  # - free.cores    the number of CPU cores that cannot be used by the function
  # - seed          random seed value used for reproducibility
  
  # Set random number seed for reproducibility
  set.seed(seed);
  
  # Set up the CPU-cluster
  cl <- makeCluster(detectCores()-free.cores);
  registerDoSNOW(cl);
  clusterExport(cl, c("getSingleAttribute", "getMultipleAttributes"));
  
  # Multi-threaded/parallel simulations
  results <- parSapply(cl, 1:n.runs, function(run) {
  
    ## Initialization
    
    # Load the required packages 
    # - notice "simmer.plot" is not required any longer, so remove the "plot(...trajectory...)" 
    #   statements from the code, if still present
    library(simmer);
    library(fitdistrplus);
    
    
    
    ## Data analysis
    
    # Load the dataset
    load("trial_dataset.RData");
    
    "...parameter definitions..."
    
    "...data analysis..."
    
    
    
    ## Supportive functions
    
    "...function definitions..."
    
    
    
    ## BSC Model
    
    "...bsc.model..."
    
    
    
    ## EXP Model
    
    "...exp.model..."
    exp.model <- trajectory()%>%
      
      # Initialization: do not forget to initialise these cycle attributes or the patients will all die.
      set_attribute(key="Tx1.Cycles", value = 0) %>%
      set_attribute(key="Tx2.Cycles", value = 0) %>%
      set_attribute(key="Tx1.Complications", value=0) %>% #0 for no complications, 1 for minor, 2 for major
      set_attribute(key="Tx2.Complications", value=0) %>%
      set_attribute(key="position", value=0) %>%
      set_attribute(key="qalys", value=0) %>%
      set_attribute(key="Alive", value=1) %>%                                                                          # define an attribute to check whether the patient is alive
      
      # First-line treatment
      set_attribute(key="Tx1.Event", value=function() Tx1.Event(cycles = get_attribute(exp.sim, "Tx1.Cycles"))) %>%                                                 # select the event to happen in this treatment cycle          
      branch(option=function() get_attribute(exp.sim, "Tx1.Event"), continue=c(T, F, F, T, T), #don't forget 5th is false
             
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
               set_attribute(key="Alive", value=0),                                                                     # update that the patient has died
             
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
               set_attribute(key="Alive", value=0)
      ) %>%
      #Second line treatment
      set_attribute(key="Tx2.Event", value=function() Tx2.Event(cycles = get_attribute(exp.sim, "Tx2.Cycles"))) %>%                                                 # select the event to happen in this treatment cycle          
      branch(option=function() get_attribute(exp.sim, "Tx2.Event"), continue=c(T, F, F, T, T), #the 5th has to be false
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
               set_attribute(key="Alive", value=0),                                                                   # update that the patient has died
             
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
               release(resource="Fu2", amount=1)
      )
    
    # Visualize to check whether the defined model structure is ok
    plot(exp.model)
    
    ## Simulations
    
    # Simulation settings
    # - notice that we must not set the random seed value here again as we did in the normal analysis
    # - notice that we must not set the number of patients to be simulated here as we did in the 
    #   normal analysis, because this setting is already provided by the "n.patients" argument in the
    #   "runPSA" function
    mon.patients <- 2;    # level of monitoring (see add_generator)
    
    "...bsc.sim..."
    
    "...exp.sim..."
    
    # Define simulation (exp)
    exp.sim <- simmer() %>%
      add_resource(name="Tx1", capacity=Inf, mon=F) %>%
      add_resource(name="Tx2", capacity=Inf, mon=F) %>%
      add_resource(name="Fu1", capacity=Inf, mon=F) %>%
      add_resource(name="Fu2", capacity=Inf, mon=F) %>%
      add_generator(name_prefix="Patient", trajectory=exp.model, distribution=at(rep(x=0, times=n.patients)), mon=mon.patients)
    
    # Run the exp simulation
    exp.sim %>% 
      run()
    
    
    
    ## Outcomes
    
    "...bsc.out..."
    
    "...exp.out..."
    
    # Calculate average outcomes
    costs.bsc <- "...write your code...";
    costs.exp <- "...write your code...";
    effect.bsc <- "...write your code...";
    effect.exp <- "...write your code...";
    
    # Remove large object to save memory
    rm(bsc.model, exp.model, bsc.sim, exp.sim, bsc.out, bsc.exp);
    
    # Return outcomes of interest, e.g. costs and effects
    return(c(costs.bsc=costs.bsc,
             costs.exp=costs.exp,
             effect.bsc=effect.bsc,
             effect.exp=effect.exp));
    
  })
  
  ## Return results ====
  
  return(results)
  
} # funtion runPSA



## Section 3: Run simulations

psa.out <- runPSA(n.patients=100, n.runs=10);
