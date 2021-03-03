# IsotopeMixingModel_GSB
Code for running the Giant Sea Bass mixing model in JAGS

For my second thesis chapter I am exploring the trophic ecology of Giant Sea Bass (GSB) using stable isotope analyses. As a part of this project I wanted to determine what proportion of GSB diet is derived from kelp versus phytoplankton primary production sources, and how/if those proportions change with size. I ran a stable isotope mixing model with a continuous fixed effect of length using the MixSIAR package in R developed in part by my advisor, Brice Semmens, and my former labmate, Brian Stock. MixSIAR relies on JAGS to run the mixing models under the hood. To make sure I had an in depth understanding of the model and to demonstrate my understanding of JAGS, I created my own JAGS script to run the model using the MixSIAR functions/output as a starting point.

This repository includes the code and data to run that mixing model.
