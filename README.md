# detectree

This is the repository for Sebastian Hickman's AI4ER MRes project, titled 'Detecting changes in tall tree height with machine learning, LiDAR, and RGB imagery'.

Its key components are: scripts to read in and tile geospatial data, an implementation of Mask R-CNN from Detectron2 (Wu et al., 2019) to perform tree crown delineation from RGB imagery, scripts to delineate tree crowns from LiDAR data using UAVforestR (T. Swinfield, https://github.com/swinersha/UAVforestR), and scripts to analyse the growth and mortality of identified trees from repeat observations. The code includes notebooks written in Python and R scripts.

The data used to evaluate and test Mask R-CNN and ITCfast are freely available at https://zenodo.org/record/5090039.

## Workflow

The workflow of the project is described by the following image.

<img width="1000" alt="workflow" src= https://github.com/shmh40/detectree/blob/main/imgs/workflow_120721.PNG > 

## Repository structure
```
├── LICENSE
|
├── README.md          <- The top-level README for developers using this project.
|
├── requirements       <- Directory containing the requirement files.
│
├── tiling             <- R scripts to read and tile geospatial data for subsequent analysis.
│   │
│   └── preprocessing  <- Scripts to convert raw RGB tiffs into tiled pngs, and to convert shapefiles to geojsons
│
├── bayesian_opt
|   |                  
│   └── bo             <- Notebook to carry out Bayesian optimisation on the parameters of the ITCfast algorithm
|                         and script to run and evaluate ITCfast with the suggested parameters. 
|
├── models             <- Notebooks to train and test models
|   |                  
│   ├── mask_rcnn      <- Notebook to train, test and evaluate a Mask R-CNN model.
|   |
|   └── itcfast        <- Notebook to test ITCfast on an area of interest, and evaluate its performance.
│
|
└── change_analysis    <- Notebooks to analyse changes in individual trees between repeat observations. 
```

## Deploying Mask R-CNN

We provide pre-trained model weights, which can be used to directly predict on your area of forest. A notebook is provided to easily deploy our pre-trained model in models/mask_rcnn. Our model was trained on 999 trees in Paracou, French Guiana, which is an area of lowland tropical forest. If you would like to train your own model, using your own manually delineated crowns, this is possible with the training notebook provided in models/mask_rcnn. This may improve performance if you are predicting on a type of forest significantly different to lowland tropical forest.

Here is an example image of the predictions made by Mask R-CNN.

<img width="500" alt="workflow" src= https://github.com/shmh40/detectree/blob/main/imgs/plot_13_285520_583300.jpg > 

## Evaluating Mask R-CNN

Code to evaluate both models is also provided in models/mask_rcnn. This requires some manually delineated tree crowns in your area of interest. The model is evaluated using standard COCO metrics, including Average Precision and Average Recall.

## Customising Mask R-CNN training

If you wish to train your own model, you may wish to alter the hyperparameters used by Mask R-CNN while training. All hyperparameters are easily altered with our notebook. The key hyperparameters you may wish to vary include the depth of the ResNet backbone, the learning rate, and the batch size.

## ITCfast

ITCfast can be deployed and evaluated using the notebooks provided in models/itcfast. The parameters of the model can be easily altered, or optimised to your own data using the scripts in bayesian_opt/bo.

## Change Analysis    

Change analysis between 2014 and 2020 is carried out with the script in the change_analysis directory, which also produces plots given in the project report, similar to this example.

<img width="1000" alt="workflow" src= https://github.com/shmh40/detectree/blob/main/imgs/tree_mortality.PNG > 

## Tiling

The tiling directory provides a script to convert geospatial data, such as GeoTiffs and GeoJSONs into formats supported by Detectron2, png images and JSONs. The outputs of this script can then be passed into Mask R-CNN to train and test the model. 
