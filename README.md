Copyright (c) 2017 The Broad Institute, Inc.
All rights reserved.

# Bladder Cancer Classifier


The 2017 TCGA paper, ["Comprehensive Molecular Characterization of Muscle-Invasive Bladder Cancer. Cell. 2017 Oct"](https://www.cell.com/cell/fulltext/S0092-8674(17)31056-5) employed unbiased NMF consensus clustering of RNA-seq data (n = 408) to identify five expression-based molecular subtypes: luminal-papillary, luminal-infiltrated, luminal, basal-quamous and neuronal.  The paper reports that these five expression subytpes may stratify response to different treatments.

The 2019 paper, ["The Cancer Genome Atlas Expression Subtypes Stratify Response to Checkpoint Inhibition in Advanced Urothelial Cancer and Identify a Subset of Patients with High Survival Probability, European Urology, volume 75, 961-964 (2019)"](https://www.sciencedirect.com/science/article/abs/pii/S0302283819301605?via%3Dihub) reports on the development of single-patient subytpe classifer based on the TCGA 2017 expression-based molecular subytpes which was used to identify 11 patients in the IMvigor trials that fell into the neuronal subytpe.  Membership in the neuronal subtype was shown to be associated with a high response rate to immunotherapy.

This repository contains the classifier code (written in R), the TCGA 2017 data required by the classifer, example expression data and the classifier's output 

## License



[BSD 3-Clause License](LICENSE)