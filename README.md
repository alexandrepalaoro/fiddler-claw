# Allometric constraint and the modulation of weapon evolution by mating system in fiddler crabs
Authors: Cristian L. Klunk, Jônatas J. de Florentino, Daniel S. Caetano, Michael S. Rosenberg & Alexandre V. Palaoro <br>
Contact about code and analyses: alexandre.palaoro@gmail.com
---

### This readme has been divided in three parts. First, we will talk about file structure, then the code, the dataset.

##### File structure:

We have three folders: "code", "data", and "figures".
The "code" folder contains the codes and the Rmarkdown files. The codes end with ".R", while the Rmarkdown is ".Rmd" and the file knitted from it is in ".html".
The "data" folder contains all the data require to run our analyses. Inside, you will find the datasets in ".csv",  ".rds" files required for bootstrap analyses (see below), and ".nex" files that contains the phylogenetic trees.
The "figures" folder contains a couple of figures. We did not upload them all because they are heavy. For all the figures, check the knitted html file. 


##### Code:

There are two codes: "RateMatrixFunctions.R" and "fiddler_crab_analyses.R". <br>

- <i>RateMatrixFunctions.R</i> - This code was taken from a Slater & Friscia (2019, see reference at the end) to run the Ovaskainen test. It is called right at the beginning of the code because it provides ancillary functions not yet implemented in the <i>ratematrix</i> package.

- <i>evol_analyses.R</i> - Here you will find the code we used. There are multiple analyses with differente formats, but it is all explained within the code and the Rmarkdown that comes along.


##### Dataset:

We are uploading two data sets. One contains the data with all species - called "medidas_full.csv". The second data set is named "def_linear_01-08-25.csv". The main difference is that "medidas full.csv" has all species, while the other data set contains only the species we are able to classify in either burrow or surface mating system.

The other two files "mating_short_all_spp.csv" and "def_mating_short.csv" only contains a column with the species and another column with the mating systems. I know this could all be done in the same data frame, but I only learned to do that after I started the analyses. So, it's legacy data but still used.

We thus used all files for analyses. The dataset containing all species were used for ancestral reconstructions (as explained in the rmarkdown and our code), while the other dataset was used for the remainder of the analyses. Each row of both files is an individual. Thus, we used the averages (or standard deviation) for the analyses. 
"NA" cells represent individuals that we could not reliably estimate the measure. For instance, individuals that did not have a tubercle are listed as NA in that column. 

For the phylogenetic tree, please check the paper.


METADATA OF def_linear_01-08-25.csv

In the columns we have the variables, in rows we have the individuals. Empty cells (NA) denote individuals that did not show the morphology being measured.  

COLUMN A: sp - the species. <br>
COLUMN B: ind - the name of the file we measured with the number of the individual. <br>
COLUMN C: carapace - the width of the carapace. Measured as the distance between the widest portion of the carapace. Unit: cm. <br>
COLUMN D: claw_size - the length of the claw. Measured from the joint between the propodus and the carpus to the distal point of the propodus. Unit: cm <br>
COLUMN E: tip_adv - Mechanical advantage for the tip of the dactyl (out.lever2). A division between in.lever and out.lever2. Unitless. <br>
COLUMN F: tub_adv - Mechanical advantage for the tubercle of the dactyl (out.lever). A division between in.lever and out.lever. Unitless. <br>
COLUMN G: mating - the mating system the species was classified as. Either 'burrow' or 'surface'. <br>
COLUMN H: manus - the distance between the joint of the propodus with the carpus to a linear projection of the base of the dactyl joint. This is a linear distance, so the projection forms a perpendicular line with the measure. Unit: cm.  <br>
COLUMN I: in.lever - distance between where the apodeme attached on the dactyl to the fulcrum. Unit: cm. <br>
COLUMN K: out.lever1 - distance between the fulcrum and the tubercle. Unit: cm. <br>
COLUMN L: out.lever1 - distance between the fulcrum and the tip of the dactyl. Unit: cm. <br>


METADATA OF mating_short_all_spp.csv <br>

We used this file to calculate the likely states of species we categorized as a "mixed" mating system. <br>
To do so, lines and the first column (A) are species. <br>
The second column (B) contains whether a species is categorized as burrow mating system (1) or not (0); species that are labeled as 0.5 are mixed species. <br>
The third column (C) contains whether a species is categorized as surface mating system (1) or not (0); species with 0.5 are categorized as mixed mating system. <br>

## References

Slater, G. J., & Friscia, A. R. (2019). Hierarchy in adaptive radiation: a case study using the Carnivora (Mammalia). Evolution, 73(3), 524-539. 

## Packages 

The code was run in R software v4.5.1. <br>
Packages used: <br>
- scatterplot3d(v.0.3-44) <br>
- pander(v.0.6.6) <br> 
- sp(v.2.2-0) <br>
- RColorBrewer(v.1.1-3) <br>
- latticeExtra(v.0.6-30) <br>
- lattice(v.0.22-7) <br>
- corncob(v.0.4.2) <br>
- performance(v.0.15.0) <br>
- scales(v.1.4.0) <br>
- viridis(v.0.6.5) <br>
- viridisLite(v.0.4.2) <br>
- ratematrix(v.1.2.4) <br>
- phylolm(v.2.6.5) <br>
- ellipse(v.0.5.0) <br> 
- vioplot(v.0.5.1) <br>
- zoo(v.1.8-14) <br>
- sm(v.2.2-6.0) <br>
- mvMORPH(v.1.2.1) <br>
- subplex(v.1.9) <br>
- corpcor(v.1.6.10) <br>
- bbmle(v.1.0.25.1) <br>
- lubridate(v.1.9.4) <br>
- forcats(v.1.0.0) <br>
- stringr(v.1.5.1) <br>
- dplyr(v.1.1.4) <br>
- purrr(v.1.1.0) <br>
- readr(v.2.1.5) <br>
- tidyr(v.1.3.1) <br>
- tibble(v.3.3.0) <br>
- ggplot2(v.3.5.2) <br>
- tidyverse(v.2.0.0) <br>
- geiger(v.2.0.11) <br>
- phytools(v.2.4-4) <br>
- maps(v.3.4.3) <br>
- ape(v.5.8-1) <br>
## Sharing/access Information

The file structure and files can be seen and downloaded from: <br>
