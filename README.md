# SpatialDeadRecovery
This repository contains R and NIMBLE scripts associated with the article "Integrating dead recoveries into pen-population spatial capture-recapture" (https://esajournals.onlinelibrary.wiley.com/doi/10.1002/ecs2.3571).


The "R SCRIPTS" folder contains example R scripts for the simulation and wolverine analysis:
-the "Simulations_SpatialDeadRecovery.R" file contains R and Nimble code necessary to simulate an OPSCR dataset and run all three models compared in the article.
-the "Wolverines_SpatialDeadRecovery.R" file contains R and NImble code to load the different wolverine datasets and run the different OPSCR models using NIMBLE. It also contains code to load the model outputs and check the parameter estimates.



The "DATA" folder contains the wolverine dataset arranged as NIMBLE input objects:

- "modelCode" : the NIMBLE model specification.

- "nimConstants": a list containing NIMBLE constants used in the model (object dimensions such as number of individuals, number of detectors, number of years ...).

- "nimParams": a character vector with the names of the parameters to return MCMC samples for.

- "nimInits": a list with initial values for the different parameters in the model.

- "nimData": a list with NIMBLE input data.
    - "z": matrix of known individual states
    - "y.alive": a 3D array containing the number of detections for each individual and year (only for the number of detectors where the individual was detected (-1 codes for no detection); see Turek et al. (2020) for more details).
    - "y.dead": only for the OPSC2R and OPSCR+DR; a matrix with individuals in rows and years in columns containing the index of the detector where the individual was recovered (for the OPSC2R) or a binary indicator denoting if the individual was recovered (for the OPSCR+DR).
    - "lowerHabCoords" and "upperHabCoords": matrices denoting the upper and lower coordinates of the habitat cells where activity centers can be placed.
    - "detCounties": a vector denoting the county where each detector is situated.
    - "detCovs": a matrix denoting the GPS track length, the distance to the nearest road and the average percentage of snow cover associated with each detector.
    - "detResponse": a matrix denotingif the indivdiual was detected in the previous year or not.
    - "denCounts": a vector denoting the average number of wolverine dens associated to each habitat cell.
    - "trials": the number of sub-detectors associated to each detector; see Milleret et al. (2018) for more details.
    - "detector.xy": x and y coordinates of the detectors where individuals can be detected alive.
    - "detector.dead.xy": x and y coordinates of the detectors where individuals can be recovered.
    - "detectorIndex": indices of the detectors in the neigborhood of each habitat cell. Used in the LESS approach; see Turek et al. (2020) and Milleret et al. (2019) for more details.
    - "nDetectorsLESS": the number of detectors in the neigborhood of each habitat cell.
    - "habitatIDDet": Matrix representation of the habitat; each cell with a number > 0 corresponds to a unique habitat cell.
    - "detectorIndex.dead" indices of the dead recovery locations in the neigborhood of each habitat cell.
    -"nDetectorsLESS.dead":the number of dead recovery locations in the neigborhood of each habitat cell.
    - "yDets": 3D array denoting  the index of the detectors where an individual was detected alive each year (-1 codes for no detection).            
    _ "nbDetections": a matrix containing the number of detectors where an individual was detected alive each year.



The "WOLVERINE OUTPUT" folder contains the MCMC samples and process results.
Each file contains a lits (named "results") with objects:
- "sims.list" : a list containing MCMC samples for all parameters in the model.
- "mean": a list containing the  mean of the posterior estimates.
- "sd": a list containing the standard deviation of the posterior estimates.
- "q2.5": a list containing the 2.5% quantile of the posterior estimates.
- "q25": a list containing the 25% quantile of the posterior estimates.
- "q50": a list containing the 50% quantile of the posterior estimates.
- "q75": a list containing the 75% quantile of the posterior estimates.
- "q97.5": a list containing the 97.5% quantile of the posterior estimates.
- "overlap0" a list indicating if 0 falls in the 95% credible interval for all parameters in the model.
- "f": a list containing the proportions of the posterior with the same sign as the mean.
- "Rhat": a list containing the potential scale reduction factor (Gelman-Rubin diagnostic).
- "n.eff": a list containing the estimated effective sample sizes for all parameters in the model.
- "colnames.sims": a list containing the names of all parameters in the model.



The "FUNCTIONS" folder contains R scripts for NIMBLE custom functions used in the different OPSCR models.
