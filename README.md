# ut-health-econ-modelling-q3
A repo to share code for the assignment on the UT health economic modelling course. By Alex and Rusheel.

To run the model on your computer: 

1. Copy the URL and open it as a new project in RStudio (select the From Version Control option).
2. That’s it! The file paths are all relative so there is no need to fiddle around with setwd. The data is exported to the Data directory in csv or .RData format.
3. For the PSA data analysis, you can simply open the HTML file in your browser.
4. The BCA startcode file in Part 1 contains the simulation for 10,000 patients and is meant to determine patient-level uncertainty. The PSA startcode in Part 2 contains the simulation for determining parameter uncertainty and performing a probablistic sensitivity analysis. 
5. The data analysis for each part is generally performed in a separate file, namely step-1.4.R and Dx-distribution-fitting.R. The parameters are loaded into the global environment for the simulation code using the source function. NOTE: For the PSA code, all the required parameters have to be exported using the clusterExport function otherwise the multithreaded code won’t work. Be aware of this in case you try to edit the code.
