# TO DO

1. In both BCA and PSA start code, followup1.time returns fixed values. Did we fit a distribution to this?
2. We should ideally determine the age by gender as the male and female ages have different distributions.
3. Tx1.time and Tx2.time in BCA return fixed values. Isn’t this patient rather than parameter uncertainty?
4. Likewise, in BCA we do not use beta distributions for the probability of death. 
5. The way we are modelling the uncertainty in costs is that the function returns different values every time it is called, meaning that for every patient, and for every cycle, the value will be a little different. This is not necessarily wrong, but it’s worth mentioning. We might ideally want the function to only return different values  *per simulation run*. 
6. Are Tx1.Time and Tx2.Time supposed to return 30 or draw from a distribution?

Once we have resolved the issues above, we can run the PSA simulation with a large number of runs, and repeat the BCA similation if needed.
