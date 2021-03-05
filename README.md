# IsotopeMixingModel_GSB
Code for running the Giant Sea Bass mixing model in JAGS

For my second thesis chapter I am exploring the trophic ecology of Giant Sea Bass (GSB) using stable isotope analyses. As a part of this project I wanted to determine what proportion of GSB diet is derived from kelp versus phytoplankton primary production sources, and how/if those proportions change with size. I ran a stable isotope mixing model with a continuous fixed effect of length and random effect of GSB individual using the MixSIAR package in R developed in part by my advisor, Brice Semmens, and my former labmate, Brian Stock. MixSIAR relies on JAGS to run the mixing models under the hood. By exploring the code underlying the MixSIAR functions and model output I was able to create a script to run the model directly. This allows me to tinker with the model, and ensure that I have in depth understanding of what the model is doing. 

This repository includes the code and data to run that mixing model directly in JAGS and also in MixSIAR.
