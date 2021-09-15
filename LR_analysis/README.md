# Data preparation

Detailed information and code for data preparation can be found in
`cellphonedb_prep.Rmd` or `cellphonedb_prep.html`. 

## Brief summary

In total, three LR interaction analyses were conducted. 

    1. Cluster 1 vs background (all other clusters) for seropositive patients
    2. Cluster 1 vs background (all other clusters) for seronegative patients
    3. Cluster 1 vs cluster 1 between seropositive and seronegatove patients

## CellphoneDB

All analyses were conducted using python 3.8.2 and CellphoneDB version 2.1.7.

Example command:

```
cellphonedb method statistical_analysis \
   --counts-data hgnc_symbol \
   --project-name=RA1to3_cluster_1_vs_all \
   --output-path all_clusters \
   --threads 20 \
meta_RA1to3_cluster_1_vs_all.txt counts_RA1to3_cluster_1_vs_all.txt
```