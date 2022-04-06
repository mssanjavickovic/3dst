# 3dst

This is a public repository for all code connected to 3D spatial transcriptomics analysis in the rheumatoid arthitis synovium. 

Please cite: Vickovic et al. Three-dimensional spatial transcriptomics uncovers cell type localizations in the human rheumatoid arthritis synovium. Communications Biology volume 5, Article number: 129 (2022). 


# Tech workflow
![github-small](https://github.com/klarman-cell-observatory/3dst/blob/master/workflow.png)

# File Structure Overview
All processed files are available at: https://singlecell.broadinstitute.org/single_cell/study/SCP460

![github-small](https://github.com/klarman-cell-observatory/3dst/blob/master/files.png)

We recommed using the `Bulk Download` function and to consult the file descriptions as mentioned bellow. 

#### `*Raw.exp*` files: matrix with raw UMI-filtered counts tsv files with:
`X_Y` barcode (X_Y) coordinate as column names
`GENE` as row names 

#### `*Norm.exp*` files: matrix with normalized counts tsv files with:
`X_Y` barcode (X_Y) coordinate as column names
`GENE` as row names 

#### `*all_inf*` files: files connect (x,y) coordinates to infiltrate annotation regions in `3dst` with:
`SectionID_X_Y` barcode (X_Y) coordinate  
`infiltrate` infiltrate ID 

#### `*HE*` files: jpg H&E image files

#### `*Mask*` files: jpg spots image files

#### `*xy*` files: files connect (x,y) coordinates to 2D centroids of actual ST spots in `3dst` and are used to create a mask during segmentation with:
`X_Y` barcode (X_Y) coordinate  
`pix_x` centroid pixel (x) coordinate in the HE and Mask images
`pix_y` centroid pixel (y) coordinate in the HE and Mask images

#### `*xyz*` files: files connect (x,y) coordinates to z-axis alignments and enable the 3D viewing:
`SectionID_X_Y` Section number followed barcode (X_Y) coordinate  
`V1` rotated and aligned centroid (x) coordinate in the HE and Mask images
`V2` rotated and centroid (y) coordinate in the HE and Mask images

#### `*Cluster*` files: files connect (x,y) coordinates to normalized expression of genes found in each respective cluster with:
`Image ID` eg. RA1_HE_1 which matches the `*HE*` file names
`x` barcode (X) coordinate 
`y` barcode (Y) coordinate  
`GENE` list of genes with matched expression values 

#### `*scRNAseq*` files: files connect (x,y) coordinates to normalized imputed expression of cell types found in each respective X_Y coordinate with:
`Image ID` eg. RA1_HE_1 which matches the `*HE*` file names
`x` barcode (X) coordinate 
`y` barcode (Y) coordinate  
`b_cell` eg. B cell specific score (same format for all 13 tested cell types)

# Pairing HE image and sequencing data
Please refer to our github repo SpoTter. 

# Segmentation
This is [code](./segmentation) for segmenting HE nuclei and cytoplasm. HE image segmentation was performed by combining Ilastik and CellProfiler. The labeled segmentation mask was used to assign the individual spots to the corresponding Cell ID. The output CSV file includes Cell IDs, X and Y position of the cells (centroid) and the corresponding spots.

# Cell typing 
This is [code](./cell_typing) for imputing cell types onto (x,y) spatial positions based on scRNA-seq data taken from [Stephenson et al](https://www.nature.com/articles/s41467-017-02659-x). [GO enrichment](./GO_enrichment) was performed on the selected genes. 

# Clustering and Differential expression (DE) analysis
This is [code](./de_analysis) for clustering and DE analysis between clusters.

# 3d app for viewing the data
This is the [app](https://spatialtranscriptomics3d.shinyapps.io/3DSeTH_RA/) for viewing 3D aligned and normalized ST-based expression. The raw [code](./app) is also available. 

# ST histoCAT integration
histoCAT interactive sessions with ST clusters and cell type annotations are available [here] (https://singlecell.broadinstitute.org/single_cell/study/SCP460)
The corresponding files can be load as "session" into histoCAT and visualized. For more information check histoCAT [documentation](https://github.com/BodenmillerGroup/histoCAT/releases/download/histoCAT_1.76/histoCATmanual_1.76.pdf).

![](ST_histoCAT.gif)
