# 3dst

This is a public repository for all code connected to 3D spatial transcriptomics analysis in the rheumatoid arthitis synovium. 

Please cite: Vickovic et al. Multidimensional transcriptomics uncovers spatial dynamics in the rheumatoid arthritis synovium. bioRxiv. 


# Tech workflow
![github-small](https://github.com/broadinstitute/hdst/blob/master/workflow.png)

# File Structure Overview
All processed files are available at: https://portals.broadinstitute.org/single_cell/study/SCP446

![github-small](https://github.com/broadinstitute/hdst/blob/master/files.png)

We recommed using the `Bulk Download` function and to consult the file descriptions as mentioned bellow. 

#### `*exp_matrix*` files: matrix with raw UMI-filtered counts tsv files with:
`X_Y` barcode (X_Y) coordinate as row names
`GENE` as column names 

#### `*norm.exp*` files: matrix with normalized counts tsv files with:
`X_Y` barcode (X_Y) coordinate as row names
`GENE` as column names 

#### `*all_inf*` files: files connect (x,y) coordinates to infiltrate annotation regions in `3dst` with:
`X_Y` barcode (X_Y) coordinate  
`Inf` infiltrate ID 

#### `*xy*` files: files connect (x,y) coordinates to centroids of actual ST spots in `3dst` and are used to create a mask during segmentation with:
`X_Y` barcode (X_Y) coordinate  
`pix_x` centroid pixel (x) coordinate in the HE and Mask images
`pix_y` centroid pixel (x) coordinate in the HE and Mask images

#### `*HE*` files: jpg H&E image files

#### `*Mask*` files: jpg spots image files

#### `*Cluster*` files: files connect (x,y) coordinates to normalized expression of genes found in each respective cluster with:
`Image ID` eg. RA1_HE_1 which matches the `*HE*` file names
`x` barcode (X) coordinate 
`y` barcode (Y) coordinate  
`GENE` list of genes with matched expression values 

#### `*scRNAseq*` files: files connect (x,y) coordinates to normalized imputed expression of cell types found in each respective X_Y coordinate with:
`Image ID` eg. RA1_HE_1 which matches the `*HE*` file names
`x` barcode (X) coordinate 
`y` barcode (Y) coordinate  
`b_cell` eg. B cell specific score (same format for all 12 tested cell types)

# Alignment
This is [code](./alignment) for aligning HE images to (x,y) barcode coordiantes as given by ST Pipeline ([v.1.5.1](https://github.com/SpatialTranscriptomicsResearch/st_pipeline/releases/tag/1.5.1)). 

# Segmentation
This is [code](./segmentation) for segmenting HE nuclei and cytoplasm. HE image segmentation was performed by combining Ilastik and CellProfiler. The labeled segmentation mask was used to assign the individual spots to the corresponding Cell ID. The output CSV file includes Cell IDs, X and Y position of the cells (centroid) and the corresponding spots.

# Cell typing 
This is [code](./cell_typing) for imputing cell types onto (x,y) spatial positions based on scRNA-seq data taken from ([Stephenson et al](https://www.nature.com/articles/s41467-017-02659-x)). 

# Clustering and Differential expression (DE) analysis
This is [code](./clustering) for clustering and DE analysis between clusters.

# 3d app for viewing the data
This is the ([app](https://spatialtranscriptomics3d.shinyapps.io/3DSeTH_RA/)) for viewing 3D aligned and normalized ST-based expression.
