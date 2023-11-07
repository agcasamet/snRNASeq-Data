# snRNASeq-Data Analysis
Data of snRNA-Seq analysis related with the paper "Tumor-induced alterations in single-nucleus transcriptome of atrophying muscles indicate enhanced protein degradation and reduced oxidative metabolism"

# Summary
Single-nucleus RNA-sequencing of cachectic muscles from tumor-bearing mice.
We report that remote tumor growth in mice impacted the transcriptome of muscle tissue, which exhibited the enrichment of type IIb myonuclear gene signatures, the activation of muscle atrophy-related mechanisms, including EDA2R, and the suppression of pathways associated with muscle contractility and oxidative metabolism. 

# Experimental Design
Syngeneic C57BL/6 mice injected with Lewis Lung Carcinoma (LLC) cells were sacrificed 16 days after tumor inoculation while they experienced moderate cachexia and loss of muscle mass and function. We investigated single-nucleus transcriptomes of the tibialis anterior (TA) muscle from tumor-bearing mice and their non-tumor-bearing controls

# Protocols
- Mice were housed at 22°C and under 50% humidity with 12 hours of light and 12 hours of dark cycles (07:00 – 19:00) and given ad libitum access to a standard rodent chow diet and water in Koc University Animal Research Facility in accordance with institutional policies and animal care ethics guidelines.

- 8-12 weeks old male mice with C57BL/6 background and Lewis lung carcinoma (LLC) cells were used for tumor inoculation. LLC cells were cultured in DMEM (Sigma, no. 5796) with 10% fetal bovine serum (FBS) and penicillin/streptomycin (Invitrogen). 5 x 10<sup>6</sup> LLC cells were injected subcutaneously over the flank while control mice received PBS only. Muscle tissues were harvested at 16 days after LLC inoculation.

- TA muscle samples dissected from 6 mice in each group were combined and processed together. Samples were chopped with dissection scissors and placed into a lysis buffer containing 10 mM Tris-HCl, 10 mM NaCl, 3 mM MgCl2, and 0.1% NP40 in nuclease-free water. Samples were homogenized using a douncer and filtered with 70 µm and 40 µm cell strainer. After centrifugation for 5 min at 500 g at 4 °C, supernatant was discarded and nuclei were resuspended and stained with DAPI and subjected to fluorescence activated cell sorting.

- 10X Genomics applications were used following the manufacturer’s guidelines (Chromium Next GEM Single Cell 3ʹ Reagent Kits) to prepare the libraries, which were sequenced using the Illumina HiSeq X system.

# Data
- Sequencing data was first analyzed and filtered using Cell Ranger (v7.0) Single-Cell Software Suite provided by 10X Genomics. Data was counted and mapped with cellranger count function with --include-introns option for pre-mRNAs (mm10). 

- Further analysis was performed using Seurat (v4.3.0) R (v4.2.2) package on R Studio (v2022.12.0), which filters nuclei, normalizes expression data, and carries out principal component analysis for clustering and Uniform Manifold Approximation and Projection (UMAP) visualization.

- Seurat R package also performs marker assignment of the clusters and gene ontology analysis of their uniquely expressed genes. 

- ClusterProfiler (4.7.1.003) R package used for KEGG pathway enrichment analysis and GSEA software (v4.3.2) used for hallmark gene sets enrichment analysis. 

- You can acces the data from Google Drive link: https://drive.google.com/file/d/1MrI5MTk2NGLLwBz0hyi1XjMg9fr3XSNN/view?usp=sharing

# Code
You can see the code in R script file.

# Cite
Agca S, Domaniku-Waraich A, Bilgic SN, Sucuoglu M, Dag M, Dogan SA, Kir S. Tumor-induced alterations in single-nucleus transcriptome of atrophying muscles indicate enhanced protein degradation and reduced oxidative metabolism. Preprint. DOI: https://doi.org/10.1101/2023.10.26.564119
