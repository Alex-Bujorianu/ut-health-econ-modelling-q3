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
    
    
    
    ## Simulations
    
    # Simulation settings
    # - notice that we must not set the random seed value here again as we did in the normal analysis
    # - notice that we must not set the number of patients to be simulated here as we did in the 
    #   normal analysis, because this setting is already provided by the "n.patients" argument in the
    #   "runPSA" function
    mon.patients <- 2;    # level of monitoring (see add_generator)
    
    "...bsc.sim..."
    
    "...exp.sim..."
    
    
    
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
