# Characterization of Population Structure from Genetic Variation in Lactose Intolerance-Related Genes

## Overview

This project explores the genetic variation associated with lactose intolerance across different human populations. Specifically, it focuses on the genetic analysis of four key genes (LCT, FUT2, FUT3, and MCM6) known to be linked to lactose intolerance. By analyzing genetic data from representative populations within European (EAU) and African (AFR) super-populations, the study aims to characterize population structure and assess whether super-population assignment can be predicted based on genetic variation in these genes.

## Objectives

- Analyze the genetic variation in lactose intolerance-related genes across different populations.
- Perform clustering analysis to explore population structure.
- Use classification models (KNN and Random Forest) to predict super-population assignment from genetic data.

## Key Technologies

- **R**: For data analysis, clustering, and classification.
- **vcfR**: For reading and processing VCF (Variant Call Format) files.
- **PCA (Principal Component Analysis)**: To reduce dimensionality and explore the genetic structure of the populations.
- **Hierarchical Clustering**: For identifying population clusters based on genetic data.
- **KNN (k-Nearest Neighbors) and Random Forest**: Classification models used to predict super-population assignment.

## Data Description

The dataset was obtained from the 1000 Genomes Project and includes genotype information across 2,500 individuals from 26 global populations. The focus of this analysis is on four specific populations: Finnish in Finland (FIN) and British in England and Scotland (GRB) from the European super-population, and Yoruba in Ibadan, Nigeria (YRI) and Esan in Nigeria (ESN) from the African super-population.

## Analysis Workflow

1. **Data Preprocessing and Filtering**:
   - VCF files were read into R and filtered to include only SNPs (Single Nucleotide Polymorphisms) passing quality thresholds.
   - Genotypes were numerically encoded, and Hardy-Weinberg Equilibrium (HWE) tests were conducted to filter out variants with significant deviations.

2. **Clustering Analysis**:
   - **PCA**: Performed to visualize the basic structure of the data and assess the clustering of populations.
   - **Hierarchical Clustering**: Applied to explore population clusters using Manhattan distance and various linkage methods.

3. **Classification Models**:
   - **KNN**: Implemented to classify super-populations based on genetic data. The model's accuracy was evaluated across different values of k.
   - **Random Forest (RF)**: Used to classify super-populations with parameter tuning to optimize model performance. The model's accuracy and robustness were assessed using a confusion matrix and out-of-bag (OOB) error rate.

## Results

- **PCA and Clustering**: The PCA revealed weak clustering of populations based on genetic variation. Hierarchical clustering, particularly using complete linkage, identified some distinct clusters but failed to align with super-populations effectively.
- **Classification**: The KNN model showed a 77% accuracy in predicting super-populations, while the optimized Random Forest model achieved a 98% accuracy, highlighting its superior performance in handling genetic data.

## Discussion

- The clustering methods struggled to clearly differentiate populations based on genetic data alone, suggesting the need for more sophisticated approaches or additional genetic markers.
- The Random Forest model demonstrated robustness in classifying populations, underscoring the potential of machine learning methods in genetic studies.

## Limitations and Future Directions

- The analysis was limited by the focus on a few genes and did not include linkage disequilibrium (LD) analysis.
- Future work should expand the dataset, include more genetic markers, and consider environmental factors that influence lactose intolerance.
- Exploring other machine learning techniques and conducting cross-validation will help improve model generalizability.

## References

1. Itan, Y., et al. (2010). A worldwide correlation of lactase persistence phenotype and genotypes. *BMC Evolutionary Biology*.
2. Anguita-Ruiz, A., et al. (2020). Genetics of lactose intolerance: An updated review. *Nutrients*.
3. Tishkoff, S. A., et al. (2007). Convergent adaptation of human lactase persistence in Africa and Europe. *Nature Genetics*.
4. Vuorisalo, T., et al. (2012). High lactose tolerance in North Europeans: a result of migration, not in situ milk consumption. *Perspectives in Biology and Medicine*.
5. Silanikove, N., et al. (2015). The Interrelationships between Lactose Intolerance and the Modern Dairy Industry. *Nutrients*.
6. Ransome-Kuti, O., et al. (1971). Absorption of lactose by various Nigerian ethnic groups. *Pediatric Research*.
7. Fan, S., et al. (2019). African evolutionary history inferred from whole genome sequence data of 44 indigenous African populations. *Genome Biology*.

