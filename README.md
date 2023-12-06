# Variant Prioritization Benchmark (VarPB)

## Introduction
The history of the development of machine learning systems has long been guided by the datasets and benchmarks available to researchers. This spans classic datasets such as ImageNet, which encouraged tremendous improvement in object recognition tasks, and the Stanford Question Answering Dataset (SQuAD), which played a similar role in natural language processing, as well as scientific data repositories like the Protein Data Bank without which the AlphaFold breakthrough would not have been possible. Recent work has established Multiplexed Assays of Variant Effects (MAVEs) in a similar role for tools dedicated to predicting the effects of protein variants. MAVEs provide an experimental measure of the functional impact of a variant (typically an amino acid substitution) and have become popular due in part to their ability to be freely downloaded and redistributed through efforts such as the ProteinGym. However, functional effects are not necessarily pathogenic effects. 

<details><summary>Read more</summary>
One of the most common uses of variant effect prediction tools (AlphaMissense, PrimateAI-3D, EVE, CADD, etc) is to identify pathogenic variants in the context of cancer or rare disease. However, predicting functional effects of protein-altering mutations is not equal to predicting the causative variant in a pool of ~10,000 coding variants per genome. We have observed a poor correlation between the performance of tools tasked with predicting functional variant effects in MAVE datasets versus prioritizing disease-causing variants in individuals with rare diseses. Therefore, we believe the field needs a complement to MAVEs that is specific to medical genetics: an open and freely-redistributable benchmark measuring how well a tool can identify the causal variant in a patient with a rare disease. The challenge with this is that individual-level human genetic variation data is typically kept under restricted access in order to protect the privacy of those individuals. Most exceptions to this, such as individuals in the 1000 Genomes Project, have been widely used in public allele frequency databases and as a result have far fewer ultra-rare variants than a typical human (and no "unique" variants with allele frequency of 0). Since many tools utilize allele frequency in making variant pathogenicity predictions, use of such samples would not replicate real-world use scenarios. 

Here, we release the Variant Prioritization Benchmark (VarPB) that we hope will be used analogously to MAVEs for the medical genetics task of identifying rare disease-causing pathogenic variants in individuals. We have collected tens of thousands of rare disease-causing variants from ClinVar with known modes of inheritance and combine them with over 100 individuals from the Personal Genome Project who have consented to open access and free distribution of their individual genotype data, but who are not present in any major public allele frequency databases (to the best knowledge of the authors). 
</details>

## Methods
<details><summary>Show Methods</summary>
Pathogenic variants were collected from ClinVar using its archive of variant summary files. The January variant summary files were collected for the years 2017, 2018, 2019, 2019, 2021, 2022, and 2023 as well as the October 2023 file (which was the most recent file at the time). Additional years are planned to be added as time passes. For each year, we selected the pathogenic or likely pathogenic missense variants that had been added in that year which had a one-star or greater level of support and which also were associated with an OMIM phenotype. We applied data from OMIM to separate these variants into those with dominant and recessive modes of action. This enabled the creation of time-resolved datasets with hundreds to thousands of variants being used for evaluation in each year. The purpose in creating time-resolved datasets is to allow users to ignore results for tools on particular years if they were trained in a supervised learning paradigm with data that included variants from that year of ClinVar data. 

We utilized the 107 individuals from (1) for whom genomic data could be obtained. These are all samples from the Personal Genome Project and as such are consented for genomic data redistribution. They were each sequenced to a mean coverage of approximately 100x. The samples are overwhelmingly of European ancestry, a weakness of the initial version of this benchmark which we hope to remedy in the future. All 107 genomes were sequenced by Complete Genomics and the processed variant calls were downloaded from the [Personal Genome Project website](https://my.pgp-hms.org/public_genetic_data). 

Variant prediction scores for each test variant and each missense variant in each of the 107 individuals were downloaded either from dbNSFP4.4a or by downloading the pre-computed predictions of the individual tools. All tool scores were normalized to a range of 0 to 1 with values closer to 1 representing greater likelihood of pathogenicity. In cases where a tool emitted multiple scores for the same mutation (for example, if the mutation affected multiple transcripts), the highest (most deleterious) prediction value was used. A total of 52 tools have been collected and benchmarked so far. The scores for 46 tools were collected from dbNSFP4.4a, while scores for Maverick, MAPPIN, AlphaMissense, EVE, PrimateAI-3D, and ESM1b were collected separately. 

For each tool, its performance for each year of ClinVar data was derived as follows. Each pathogenic variant was individually 'spiked into' the set of variants for each individual. Next, the tool would score each missense variant in that individual (including the pathogenic one). Then, variants would be ranked in descending order by the tool's score. The rank of the causal variant within the set would be taken as the raw measure of performance of the tool for that variant in that individual. This rank value was then normalized by dividing it by the total number of variants in that individual for which this tool provided scores. In this way, the normalized score gives a sort of 'rank percentile' for the variant which normalizes out any advantage a tool might get if it only has scores available for a subset of all missense variants. This 'spike-in' process is then repeated for each of the 31,811 pathogenic variants being placed into each of the 107 individuals for a total of over 3.4 million tests per tool. These normalized rank values are then plotted to show the cumulative number of simulated cases that would have been solved within the top, say, 0.01% to 0.3% of variants within an individual. This method of evaluation is meant to show how good of a job the tool is doing at 'picking the needle out of the haystack'. Finally, we take the area under this cumulative cases solved curve in order to provide a final score for each tool. A perfect tool would get an area under the curve score of 1, meaning that it gives top rank to the causal variant every single time. While a poor tool would get a score of 0, meaning that it takes it longer to find any of the causal variants than the window of evaluation here (typically, the tool did not prioritize the causal variant into the top 1% of variants in the individual). 

For tools that give separate scores for dominant vs recessive variants (currently just Maverick and MAPPIN), the genotype of the variant in the individual is used to determine which value should be used. Homozygous variants are assigned the recessive predicted score, while heterozygous variants are assigned the dominant predicted score. Additionally, all possible pairs of heterozygous variants on the same gene also get their recessive scores averaged in order to consider the possibility of compound heterozygous variant effects. For the causal variants, those with dominant effects are spiked-in as heterozygous, while those with recessive effects are spiked-in as homozygous. When normalizing performance, these 'zygosity-aware' tools have their causal variant rank divided by the number of variants plus the number of compound heterozygous pairs that they considered in each individual. Also for each of these zygosity-aware tools, a non-zygosity-aware version is modelled where the dominant and recessive pathogenic scores are summed together and the tool is evaluated exactly the same way as all the other tools. Those are listed as Maverick (No Zygosity) and MAPPIN (No Zygosity) in the results set. 

### Example Results
Results for each individual year generate a "solve curve" like this: ![2022 non-normalized ranks](Figures/2022_ranks20.png)
The x-axis is the top-n ranks of variants inspected and the y-axis is the cumulative percentage of simulated cases that would be solved by variants in those top-n ranks.

We then normalize the ranks as described above, which generates a normalized solve curve for that year: ![2022 normalized ranks](Figures/2022_normalizedRanks.png)
This figure is much the same, but the x-axis is now the percentile of variants instead of absolute rank number in order to not give an unfair advantage to tools that did not have available scores for all variants. 

Finally, we calculate the area under those normalized curves to generate our final score for each tool for each individual year (only the top 40 tools are shown for visual clarity): ~[2022 normalized AUC](Figures/2022_AUC_barplots_normalizedRanks.png)

### VarPB correlation with MAVEs
Our initial motivation for creating VarPB was anecdotal observation that strong performance on MAVE datasets did not equate to strong performance at medical genetics prioritization tasks. So, one of our top priorities with VarPB was to assess its correlation with MAVE performance measures for a variety of tools. We used the performance of 10 pathogenicity prediction tools (AlphaMissense, EVE, gMVP, VARITY_R_LOO, REVEL, SIFT, Polyphen2_HVAR, Polyphen2_HDIV, CADD, and PrimateAI) on the 26 human gene datasets in ProteinGym provided in Figure 3B of Cheng, et al (2) as our measure of MAVE performance and plotted this against performance measured on VarPB: ![varPB vs MAVEs](Figures/MAVE_comparison.png)

We observed a spearman correlation coefficient of 0.202 between these two variables, suggesting that VarPB is indeed measuring a distinct aspect of variant effect prediction performance than this set of MAVEs. 

</details>

## Results
### Overall Results
The current overall result of VarPB is shown below: ![Overall Result](Figures/Overall_Result.png)

This figure shows the average area under the normalized solve curve for each tool measured over the seven years of data (2017 - 2023) with the error bars representing the standard deviation of the tool's performance over the years. Supervised learning approaches are shown in grey, while unsupervised methods are shown in black. 

### Individual Year VarPB Data
<details><summary>Show individual year results</summary>
In this section, the results for each individual year of ClinVar variants are reported. We show this with normalized solve curve, normalized area under the solve curve bar plot, non-normalized solve curve for top 20 ranks, non-normalized area under the solve curve for top 20 ranks bar plot, non-normalized solve curve for top 100 ranks, and non-normalized area under the solve curve for top 100 ranks bar plot for each year. 
#### 2023
![2023 normalized solve curve](Figures/2023_normalizedRanks.png)

#### 2022
![2022 normalized solve curve](Figures/2022_normalizedRanks.png)

#### 2021
![2021 normalized solve curve](Figures/2021_normalizedRanks.png)

#### 2020
![2020 normalized solve curve](Figures/2020_normalizedRanks.png)

#### 2019
![2019 normalized solve curve](Figures/2019_normalizedRanks.png)

#### 2018
![2018 normalized solve curve](Figures/2018_normalizedRanks.png)

#### 2017
![2017 normalized solve curve](Figures/2017_normalizedRanks.png)
</details>

## Discussion
In the current version of the benchmark, Maverick (3) outperforms all other tools in all seven years evaluated. It is a supervised learning system and was trained on a set that included the variants used in evaluation years 2017, 2018, and 2019, so its performance on those years should be acknowledged to be inflated (Maverick's overall performance as shown in the first figure would drop from 0.808 to 0.787 if years 2017, 2018, and 2019 were removed, so it would still be the top ranking tool). 

In general, supervised learning methods largely outperformed the unsupervised methods like AlphaMissense, PrimateAI-3D, and EVE. While those tools exhibit strong performance in MAVE benchmarks, there is clearly still information being conveyed in the supervised learning regime that they have not been able to replicate in the unsupervised methodology. This represents the biggest area for improvement in the field. 

This benchmark includes only missense variants, but future versions should offer tests that include frameshifting and non-frameshifting indels as well as stop-gain and stop-loss variants. Eventually, benchmarks like this should expand beyond just protein-altering variants to consider non-coding variants, structural variants, and tandem repeats. It is the authors' hope that a benchmark like this can grow with the field to continually provide a challenge for tool developers to solve. 

## Data Access
All ClinVar variants used and the missense variants from each test individual, along with the ranks by which each tool solved each combination of test variant and individual, are available for download from [Zenodo](https://doi.org/10.5281/zenodo.10278481). Using this data, individuals can evaluate a new tool and compare it against all the others. Code to perform these evaluations is available in the [code folder](https://github.com/ZuchnerLab/VariantPrioritizationBenchmark/tree/main/Code). Unfortunately, the individual prediction scores of each tool on each variant are not included in this benchmark as they are subject to each individual tool's own copyright restrictions. 

The leaderboard of results data is available in the [results folder](https://github.com/ZuchnerLab/VariantPrioritizationBenchmark/tree/main/Results).

## References
1. Mao Q, Ciotlos S, Zhang RY, Ball MP, Chin R, Carnevali P, Barua N, Nguyen S, Agarwal MR, Clegg T, Connelly A, Vandewege W, Zaranek AW, Estep PW, Church GM, Drmanac R, Peters BA. The whole genome sequences and experimentally phased haplotypes of over 100 personal genomes. Gigascience. 2016 Oct 11;5(1):42. doi: 10.1186/s13742-016-0148-z. PMID: 27724973; PMCID: PMC5057367.
2. Cheng J, Novati G, Pan J, Bycroft C, Žemgulytė A, Applebaum T, Pritzel A, Wong LH, Zielinski M, Sargeant T, Schneider RG, Senior AW, Jumper J, Hassabis D, Kohli P, Avsec Ž. Accurate proteome-wide missense variant effect prediction with AlphaMissense. Science. 2023 Sep 22;381(6664):eadg7492. doi: 10.1126/science.adg7492. Epub 2023 Sep 22. PMID: 37733863.
3. Danzi MC, Dohrn MF, Fazal S, Beijer D, Rebelo AP, Cintra V, Züchner S. Deep structured learning for variant prioritization in Mendelian diseases. Nat Commun. 2023 Jul 13;14(1):4167. doi: 10.1038/s41467-023-39306-7. PMID: 37443090; PMCID: PMC10345112.

