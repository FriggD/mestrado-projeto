### Column-wise missingness
Column-wise missingness refers to the proportion of missing values within each column of a dataset, which is crucial for understanding the extent and pattern of missing data. This metric can be calculated by dividing the number of missing entries in a column by the total number of entries in that column. Analyzing column-wise missingness helps in selecting appropriate imputation techniques tailored to each feature, enhancing data quality and model performance.
Calculation Method

    Formula: Column-wise missingness can be calculated using the formula: [ \text{Missingness} = \frac{\text{Number of Missing Values}}{\text{Total Number of Values}} \times 100 ]
    Example: If a column has 10 missing values out of 100 total entries, the column-wise missingness is 10%10%.

Importance in Data Imputation

    Guided Imputation: Techniques like Column-wise Guided Data Imputation (cGDI) utilize this metric to select the best imputation method for each feature, improving accuracy(Petrozziello & Jordanov, 2016).
    Pattern Recognition: Understanding missingness patterns aids in diagnosing data quality issues and selecting appropriate treatment methods(Bannon, 2015).

While column-wise missingness provides valuable insights, it is essential to consider the overall context of missing data, including the mechanisms behind it, such as Missing at Random (MAR) or Missing Not at Random (MNAR), which can significantly influence the choice of imputation methods(Finucane et al., 2024).

### k-nearest SNP masking
The concept of k-nearest SNP masking involves using the k-nearest neighbors (KNN) algorithm to identify and mask single-nucleotide polymorphisms (SNPs) that may interfere with genetic analysis, particularly in gene expression studies. This approach is crucial for improving the accuracy of genetic data interpretation by removing SNPs that could cause hybridization issues in microarray analyses. The calculation of k-nearest SNP masking involves determining the k-nearest neighbors for each SNP and evaluating their impact on classification or expression analysis. This process can be enhanced by optimizing the number of neighbors (k) and the weighting of SNPs based on their relevance to the classification task.
K-Nearest Neighbors Algorithm

    The KNN algorithm is used to classify data points based on the majority class of their k-nearest neighbors.
    In the context of SNP masking, KNN can help identify SNPs that are noise by iteratively adjusting their weights until convergence, where most non-informative SNPs have zero weights(Gu & Ding, 2021).

SNP Masking Approach

    SNP masking involves removing SNP-affected probes to mitigate hybridization problems in microarray analyses.
    This approach is particularly useful in detecting gene expression differences among genetically unique individuals or strains(Walter et al., 2007) (Walter et al., 2007).

Optimization of K in KNN

    The choice of k is crucial and can be optimized based on factors like class imbalance and feature scoring metrics.
    Methods such as constrained variable-wise optimized k (VWOK) and fixed-k derived with principal components analysis (kPCA) can improve feature detection and classification performance(Dawkins & McKinney, 2024).

While k-nearest SNP masking is a valuable tool for improving genetic analysis, it is important to consider the potential limitations of the KNN algorithm, such as sensitivity to the choice of k and the distance metric used. Additionally, the effectiveness of SNP masking may vary depending on the specific genetic data and analysis context.

### Minor Allele Frequency (MAF) masking 

Minor Allele Frequency (MAF) masking refers to the practice of excluding genetic variants from analysis based on their frequency in a population, specifically those with low MAF. This approach is crucial in genomic studies as it can significantly influence the outcomes of population structure inference and genetic association studies. The calculation of MAF involves determining the frequency of the less common allele at a given genetic locus within a population.
Calculation of Minor Allele Frequency

    Definition: MAF is calculated as the number of copies of the minor allele divided by the total number of alleles at that locus.
    Formula: MAF = (Number of minor alleles) / (Total number of alleles) = (2 * Number of heterozygotes + Number of minor homozygotes) / (2 * Total number of individuals).
    Thresholds: Researchers often set MAF thresholds to filter out rare variants, which can confound analyses. For instance, excluding alleles with MAF < 0.05 is common to enhance data quality(Linck & Battey, 2019)(Linck & Battey, 2017).

Implications of MAF Masking

    Population Structure: MAF thresholds can distort population structure inference, leading to less distinct clusters when stringent cutoffs are applied(Linck & Battey, 2019)(Linck & Battey, 2017).
    Genotype-Environment Interactions: MAF influences the estimation of genetic components in twin studies, potentially inflating estimates of additive genetic variance when shared environmental factors are present(Verhulst & Neale, 2016).

While MAF masking is a widely accepted practice to improve data integrity, it may inadvertently exclude valuable information about rare variants that could be critical for understanding complex traits and diseases. This highlights the need for careful consideration of MAF thresholds in genomic research.


### Haplotype cluster masking
Haplotype cluster masking is a technique used in genetic studies to refine the analysis of haplotypes by grouping them based on similarity and evolutionary relationships. This method enhances the accuracy of genotype predictions at causal loci by focusing on clusters of haplotypes that share common ancestry. The calculation of haplotype cluster masking involves defining haplotype clusters, assessing their similarity, and applying algorithms to manage missing data effectively.
Haplotype Clustering

    Definition: Haplotype clusters are formed by grouping haplotypes based on a defined similarity metric, often centered around a putative ancestral haplotype(Waldron et al., 2006).
    Algorithm Implementation: Algorithms, such as those developed by Teo and Small, allow for the inclusion of missing data and downweigh SNPs with high missingness, facilitating accurate clustering(Teo & Small, 2009).

Calculation Methodology

    Similarity Metrics: Closeness between haplotypes is determined using metrics that measure shared segments around functional mutations(Molitor et al., 2003).
    Markov Models: Novel haplotype cluster Markov models can phase genomic samples rapidly, discarding unlikely haplotype pairs to refine clusters(Ball et al., 2017).

While haplotype cluster masking significantly improves the precision of genetic analyses, it may also introduce biases if the underlying assumptions about haplotype similarity do not hold true across diverse populations. This highlights the importance of validating clustering methods in various genetic contexts.

### Random masking
Random masking is a technique used to obscure sensitive data by applying random values, or masks, to the data before processing. This method enhances security by reducing the correlation between the sensitive data and any observable side-channel information. The calculation of random masking involves generating random numbers and applying them systematically to the data.
Key Aspects of Random Masking
Generation of Random Masks

    Random masks can be generated using various methods, including deterministic random number generators (RNG) that produce random numbers of specified lengths(Perion et al., 2015).
    The masks are often created in sets, with each set containing multiple random numbers that can be applied to the data(Chuck et al., 2020).

Application of Masks

    Masks are applied in a structured manner, often in multiple stages, to ensure that the original data remains obscured throughout the processing(Liardet & Teglia, 2008).
    The application can involve different orders and combinations of masks to enhance security further(Liardet & Teglia, 2008).

Quantifying Masking Strength

    The effectiveness of random masking can be quantified using metrics like quantitative masking strength (QMS), which assesses the resistance of the implementation against side-channel attacks(Eldib et al., 2015).

In contrast, while random masking is effective, it can be labor-intensive and prone to errors during implementation, potentially leading to vulnerabilities if not executed correctly(Eldib et al., 2015).

### Row-wise Missingness
Row-wise missingness refers to the absence of data in specific rows of a dataset, which can significantly impact data analysis and interpretation. Calculating row-wise missingness involves determining the proportion of missing values within each row, allowing researchers to assess the extent of missing data and its potential implications for their analyses. The following sections outline the key aspects of calculating row-wise missingness.
Definition and Importance

    Row-wise missingness indicates how many values are missing in each row of a dataset.
    Understanding this metric is crucial for evaluating data quality and determining appropriate handling methods for missing data(Kostenko, 2024).

Calculation Method

    To calculate row-wise missingness:
        Count the total number of entries in each row.
        Count the number of missing entries in that row.
        Compute the missingness percentage using the formula: [ \text{Missingness} = \left( \frac{\text{Number of Missing Entries}}{\text{Total Entries}} \right) \times 100 ]
    This calculation helps identify rows with significant missing data, guiding decisions on imputation or exclusion(Bannon, 2015).

Implications for Analysis

    High row-wise missingness can lead to biased results and reduced statistical power.
    It is essential to consider the pattern of missingness (e.g., Missing at Random vs. Missing Not at Random) when deciding on treatment methods(Pham et al., 2023)(Kostenko, 2024).

While row-wise missingness provides valuable insights into data quality, it is also important to consider the broader context of missing data patterns and their implications for statistical analysis. Addressing missing data effectively requires a comprehensive understanding of both the extent and the mechanisms behind the missingness.
