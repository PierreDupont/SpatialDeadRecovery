# SpatialDeadRecovery
This repository contains R and NIMBLE scripts associated with the article "Integrating dead recoveries into pen-population spatial capture-recapture" by [Pierre Dupont](https://www.researchgate.net/profile/Pierre-Dupont-6), [Cyril Milleret](https://www.researchgate.net/profile/Cyril-Milleret), [Madhieh Tourani](https://www.researchgate.net/profile/Mahdieh-Tourani/research), [Henrik Brøseth](https://www.researchgate.net/profile/Henrik-Broseth) and [Richard Bischof](https://www.researchgate.net/profile/Richard-Bischof-2), published in Ecosphere:

[Dupont, P., Milleret, C., Tourani, M., Brøseth, H., & Bischof, R. (2021). Integrating dead recoveries in open‐population spatial capture–recapture models. Ecosphere, 12(7), e03571.](https://esajournals.onlinelibrary.wiley.com/doi/10.1002/ecs2.3571).




The **"R SCRIPTS"** folder contains example R scripts for the simulation and wolverine analysis:

- "_Simulations_SpatialDeadRecovery.R_" contains R and Nimble code necessary to simulate an OPSCR dataset and run all three models compared in the article.
- "_Wolverines_SpatialDeadRecovery.R_" contains R and NImble code to load the different wolverine datasets and run the different OPSCR models using NIMBLE. It also contains code to load the model outputs and check the parameter estimates.





The **"WOLVERINE DATA"** folder contains the wolverine dataset arranged as NIMBLE input objects:

- "_modelCode_" : the NIMBLE model specification.

- "_nimConstants_": a list containing NIMBLE constants used in the model (object dimensions such as number of individuals, number of detectors, number of years ...).

- "_nimParams_": a character vector with the names of the parameters to return MCMC samples for.

- "_nimInits_": a list with initial values for the different parameters in the model.

- "_nimData_": a list with NIMBLE input data.
    - "_z_": matrix of known individual states
    - "_y.alive_": a 3D array containing the number of detections for each individual and year (only for the number of detectors where the individual was detected (-1 codes for no detection); see Turek et al. (2020) for more details).
    - "_y.dead_": only for the OPSC2R and OPSCR+DR; a matrix with individuals in rows and years in columns containing the index of the detector where the individual was recovered (for the OPSC2R) or a binary indicator denoting if the individual was recovered (for the OPSCR+DR).
    - "_lowerHabCoords_" and "upperHabCoords": matrices denoting the upper and lower coordinates of the habitat cells where activity centers can be placed.
    - "_detCounties_": a vector denoting the county where each detector is situated.
    - "_detCovs_": a matrix denoting the GPS track length, the distance to the nearest road and the average percentage of snow cover associated with each detector.
    - "_detResponse_": a matrix denotingif the indivdiual was detected in the previous year or not.
    - "_denCounts_": a vector denoting the average number of wolverine dens associated to each habitat cell.
    - "_trials_": the number of sub-detectors associated to each detector; see Milleret et al. (2018) for more details.
    - "_detector_.xy": x and y coordinates of the detectors where individuals can be detected alive.
    - "_detector.dead.xy_": x and y coordinates of the detectors where individuals can be recovered.
    - "_detectorIndex_": indices of the detectors in the neigborhood of each habitat cell. Used in the LESS approach; see Turek et al. (2020) and Milleret et al. (2019) for more details.
    - "_nDetectorsLESS_": the number of detectors in the neigborhood of each habitat cell.
    - "_habitatIDDet_": Matrix representation of the habitat; each cell with a number > 0 corresponds to a unique habitat cell.
    - "_detectorIndex.dead_" indices of the dead recovery locations in the neigborhood of each habitat cell.
    - "_nDetectorsLESS.dead_":the number of dead recovery locations in the neigborhood of each habitat cell.
    - "_yDets_": 3D array denoting  the index of the detectors where an individual was detected alive each year (-1 codes for no detection).
    - "_nbDetections_": a matrix containing the number of detectors where an individual was detected alive each year.



The **"WOLVERINE OUTPUT"** folder contains the MCMC samples and process results.
Each file contains a lits (named "_results_") with objects:
- "_sims.list_" : a list containing MCMC samples for all parameters in the model.
- "_mean_": a list containing the  mean of the posterior estimates.
- "_sd_": a list containing the standard deviation of the posterior estimates.
- "_q2.5_": a list containing the 2.5% quantile of the posterior estimates.
- "_q25_": a list containing the 25% quantile of the posterior estimates.
- "_q50_": a list containing the 50% quantile of the posterior estimates.
- "_q75_": a list containing the 75% quantile of the posterior estimates.
- "_q97.5_": a list containing the 97.5% quantile of the posterior estimates.
- "_overlap0_" a list indicating if 0 falls in the 95% credible interval for all parameters in the model.
- "_f_": a list containing the proportions of the posterior with the same sign as the mean.
- "_Rhat_": a list containing the potential scale reduction factor (Gelman-Rubin diagnostic).
- "_n.eff_": a list containing the estimated effective sample sizes for all parameters in the model.
- "_colnames.sims_": a list containing the names of all parameters in the model.



The "**FUNCTIONS**" folder contains R scripts for NIMBLE custom functions used in the different OPSCR models.
