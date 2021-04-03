# TO DO

1. <s>In both BCA and PSA start code, followup1.time returns fixed values. Did we fit a distribution to this?</s> We fitted a distribution for the *probability* of death during followup (beta_death_followup followup1.event) but not for the time. We assume this period is fixed.
2. <s>We should ideally determine the age by gender as the male and female ages have different distributions.</s> Done!
3. <s>Tx1.time and Tx2.time in BCA return fixed values. Isn’t this patient rather than parameter uncertainty?</s> Done! We implemented this uncertainty in the BCA simulation as well, because it’s a characteristic of the patient rather than the input parameters.
4. <s>Likewise, in BCA we do not use beta distributions for the probability of death.</s> Done.
5. The way we are modelling the uncertainty in costs is that the function returns different values every time it is called, meaning that for every patient, and for every cycle, the value will be a little different. This is not necessarily wrong, but it’s worth mentioning. We might ideally want the function to only return different values  *per simulation run*. 
6. <s>Are Tx1.Time and Tx2.Time supposed to return 30 or draw from a distribution?</s> From what we understood from the PDF document, the duration of a treatment is fixed at 30 days. We only fitted distributions for time-to-death and time-to-major-complication because those make the patient stop treatment.
7. The way we handle major complications is incorrect! The patient should stop treatment with Tx1 or Tx2 if they get major complications (moving into followup 1 or palliative care). They should not continue with it.

Once we have resolved the issues above, we can run the PSA simulation with a large number of runs, and repeat the BCA simulation if needed.
