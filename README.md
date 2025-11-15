# PathLength

# Neural-trajectory-dynamics – PathLength  
Scripts for the analyses and plots of the manuscript: *[Neural trajectory dynamics during human motor learning]*  

## Overview  
This repository contains the MATLAB scripts used for computing neural trajectory dynamics and calculating path-length metrics in relation to behavioural / experimental variables, as reported in the accompanying paper.  

## Repository structure  

### `analyses/`  
Contains scripts to load the neural and behavioural data, compute neural trajectories (e.g., via PCA, manifold embedding, state-space reconstruction), compute path length metrics across time/conditions/epochs, and perform statistical comparisons.  
Key scripts might include:  
- `MainCode` — loads data, calculates gradient Maps, Smooths data & calculates projections.  
- `NetDisplacementSlidingWindow.m` — computes the Path Length.
- `PathLengthSlidingWindow.m` — computes the Net Displacement.  
- `PcslidingWindow.m` —  computes Trajectories for each principal component.
- `contrastBrainMaps.m` —  computes task-transition BOLD contrast maps, builds the representational similarity matrix, and performs MDS.
- `headMotion.m` —  extracts framewise displacement and DVARS, and fits/plots mixed-effects models relating FD, DVARS, and path length to phase and day.
- `individualDifferences.m` — quantifies individual differences in neural trajectory compression and tests how strongly they predict each subject’s brain–behaviour coupling across task conditions.

### `plots/`  
Contains scripts to generate all figures reported in the paper. Typical scripts:  
- `pAveraged t-values within each network barplot` —  
- `FCM & AffinityM & GradiantsM` — 
- `PathLengthPlots.m` —  
- `NetDisplacementPlots.m` —
- `PCs Plots.m` —
- `PCs Spatial topography.m` —
- `Representational similarity matrix & UMAP.m` —
- `Trajectory&BehavioralCorrelations.m` —
- `behavioralDataPlots.m` —
- `contrastBrainMap.m` —
- `headMotion.m` —
- `temporal trajectory.m` —

## Getting started  
### Prerequisites  

### `Prerequisites/`
Contains all software requirements and dependencies needed to run the analyses in this repository.  
Key components might include:  
- `MATLAB R2024b` — primary environment for all analyses.  
- `Signal Processing Toolbox (24.2)` — required for preprocessing and time-series operations.  
- `Statistics and Machine Learning Toolbox (24.2)` — used for statistical tests (e.g., LME models, regressions).  
- `Bioinformatics Toolbox (24.2)` — required for data-handling functions in neuroimaging workflows.  
- `GIfTI library (gifti-main)` — reads/writes GIfTI surface files.  
- `SUIT toolbox` — performs cerebellar isolation, normalization, and reconstruction.  
- `SurfStat` — used for surface-based statistical analysis.  
- `JASP (0.19.3)` — used for additional statistical tests.  


