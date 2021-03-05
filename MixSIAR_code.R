#Kayla Blincow
#Update: 3/2/2021

#The goal of this script is to run a Bayesian mixing model to look at relative
#contributions of different primary production sources to GSB diets using MixSIAR.

#Basing this off the alligator example from Stock PeerJ publication

#clear my workspace 
rm(list = ls())

#load packages
library(tidyverse)
library(MixSIAR)
library(R2jags)

# load mixture (consumer) data

mix <- load_mix_data(filename="FinalGSBBulk.csv",
                     iso_names=c("d13C", "d15N"),
                     factors="tag_ID",
                     fac_random=TRUE,
                     fac_nested=NULL,
                     cont_effects="TotalLength")

# load source data
source <- load_source_data(filename= "PPSources.csv",
                           source_factors=NULL,
                           conc_dep=FALSE,
                           data_type="means",
                           mix=mix)

# load TEF data
discr <- load_discr_data(filename="gsb_TEF.csv", mix=mix)


# Define model structure and write JAGS model file
model_filename <- paste0("MixSIAR_model_length.txt")
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

# Run the JAGS model
jags.mod <- run_model(run="short", mix, source, discr, model_filename, 
                      alpha.prior=1, resid_err, process_err)

# Process diagnostics, summary stats, and posterior plots
output_JAGS(jags.mod, mix, source,
            output_options = list(
              summary_save = TRUE,
              summary_name = "ModelSummary/summary_statistics",
              sup_post = TRUE,
              plot_post_save_pdf = TRUE,
              plot_post_name = "ModelSummary/Posteriors/posterior_density",
              sup_pairs = TRUE,
              plot_pairs_save_pdf = TRUE, 
              plot_pairs_name = "ModelSummary/pairs_plot",
              sup_xy = TRUE,
              plot_xy_save_pdf = TRUE,
              plot_xy_name = "ModelSummary/traceplot",
              gelman = TRUE,
              heidel = FALSE,
              geweke = FALSE,
              diag_save = TRUE,
              diag_name = "ModelSummary/diagnostics",
              indiv_effect = FALSE,
              plot_post_save_png = FALSE,
              plot_pairs_save_png = FALSE,
              plot_xy_save_png = FALSE)
)
graphics.off()


