# detectree

This is the repository for Sebastian Hickman's AI4ER MRes project, titled 'Detecting changes in tall tree height with machine learning, LiDAR, and RGB imagery'.

Its key components are: scripts to read in and tile geospatial data, an implementation of Mask R-CNN from Detectron2 (Wu et al., 2019) to perform tree crown delineation from RGB imagery, scripts to delineate tree crowns from LiDAR data using UAVforestR (T. Swinfield, https://github.com/swinersha/UAVforestR), and scripts to analyse the growth and mortality of identified trees from repeat observations.

The workflow of the project is described by the following image.


## Repository contents
'''
├── LICENSE
├── README.md          <- The top-level README for developers using this project.
|
├── requirements       <- Directory containing the requirement files.
│
├── tiling             <- R scripts to read and tile geospatial data for subsequent analysis.
│   |
│   ├── data_loading   <- Scripts to download tiffs of RGB and LiDAR data, and shapefiles of manual crowns
│   │
│   ├── preprocessing  <- Scripts to convert raw RGB tiffs into tiled pngs, and to convert shapefiles to geojsons
|   |
│   └── tests          <- Scripts for unit tests of your functions
│
├── bayesian optimisation
|   |                  
│   └── bo             <- Notebook to carry out Bayesian optimisation on the parameters of the ITCfast algorithm
|
├── models             <- Notebooks to train and test models
|   |                  
│   ├── mask rcnn      <- Notebook to train, test and evaluate a Mask R-CNN model on manually labelled tree crowns, and then predict on new RGB images.
|   |
|   └── itcfast        <- Notebook to test ITCfast on an area of interest, and evaluate its performance.
│
├── export
|   |                  
│   └── predictions    <- Notebook to export the predictions of the Mask R-CNN model back to a geospatial format (shapefile).
|
└── change analysis    <- Notebooks to analyse changes in individual trees between repeat observations. 
'''
