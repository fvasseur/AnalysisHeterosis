# AnalysisHeterosis
This is the R script associated with the study of phenotypic non-linearity at the origin of the emergence of heterosis in Arabidopsis thaliana

This R script describes the analysis of trait values mesured on inbred and hybrid individuals.
Traits can be downloaded from Dryad repository.
 
Traits were measured with the same protocol in two experiments that took place at Max Planck Institute for Development Biology (Tuebingen, Germany) between 2013 and 2015. 
"idExp" represents experiment identifier: GT01, which is the first experiment and which contains 451 inbreds (n=2), and GT05, the second experiment which contains 450 hybrids (n=4) and 16 inbreds (n=4). 
Identifier numbers are 1001 genomes identifiers (http://1001genomes.org/).
Detailed protocol for trait measurement can be found in "Image-based methods for phenotyping growth dynamics and fitness in Arabidopsis thaliana" (https://doi.org/10.1101/208512).

Units of traits measured:
- AgeAtReproduction : duration (days) between the emergence of the first two leaves and the end of reproduction, when fruits are drying
- FruitNumber: total number of fruits (siliques) per individual, measured at the end of reproduction
- GrowthRate: Average growth rate (mg d-1), measured as the ratio of final rosette dry mass over plant lifespan
- VegetativeDryMass: rosette dry mass (mg) measured at the growth inflection point (when growth rate is maximum), which was estimated by fitting a sigmoid growth curve on plant dry mass over time.
 
