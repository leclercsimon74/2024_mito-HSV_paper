![Image S3](https://github.com/leclercsimon74/2024_mito-HSV_paper/blob/main/Supplementary%20Figure%20S3%20(3View)/Supplementary%203View_LD.png)

**S3 Fig. Mitochondria elongate as the infection proceeds.**
**(a)** Serial block face scanning electron microscopy (SBF-SEM) images of noninfected and infected MEF cells at 8 and 12 hpi. Scale bars: 2 μm. **(b)** Quantitative analysis of mitochondrial length (nmito = 202, 64, and 244 for NI, 8, and 12 hpi, respectively). The box plots show the mean (dashed line) and the interquartile range. Statistical significance was determined using the Student’s t-test. The significance values are denoted as **** (p<0.0001), *(p<0.05), or ns (not significant).

https://doi.org/10.1371/journal.ppat.1011829.s003

# [Figure.py](https://github.com/leclercsimon74/2024_mito-HSV_paper/blob/main/Supplementary%20Figure%20S3%20(3View)/Figure.py)
## Libraries:
- `pandas`: used for data manipulation and analysis (reading CSV)
- `seaborn`: used for creating statistical visualizations (boxplots, swarmplots)
- `statannotations`: used for statistical comparisons within boxplots
- `matplotlib`: used for creating visualizations (images, plots)
- `tifffile`: used for reading and writing TIFF image files
- `numpy`: used for numerical computations (arrays)
- `scipy.ndimage`: used for image processing (erosion)
- `mpl_toolkits`: used for adding scalebars to plots

## Image Processing and Visualization (First Part):
### Load Images:
The code loads raw images, segmentation masks for nuclei, mitochondria, and cells from multiple conditions (NI_slice_30, 8hpi_slice_50, 12hpi_slice_68) using `tifffile.imread`.
### Filter Images:
It filters the raw images by applying the cell segmentation mask. This removes anything outside the cell boundaries.
The function rescale_intensity from `skimage.exposure` is used to adjust the image intensity for better visualization.
### Plot Images:
It creates a figure with 3 subplots.

Each subplot displays a filtered raw image for a specific condition (Non-infected, 8hpi, 12hpi).

Scale bars are added to each image using `AnchoredSizeBar`.


## Mitochondrial Length Analysis (Second Part):
### Load Data:
It reads data from a CSV file named `length.csv` and `Volume.csv` using pandas.read_csv. This file contains information about mitochondrial lengths and volumes for each condition.
### Boxplot and Swarmplot:
The code creates a figure with 2 subplots.
It uses `seaborn.boxplot` to visualize the distribution of mitochondrial lengths and volumes for each condition. Customization is applied to remove unnecessary elements (caps, whiskers, medians) and highlight mean values.
A swarmplot is created on top of the boxplot using `seaborn.swarmplot` to show individual data points.
### Statistical Comparison (Optional):
The code utilizes `statannotations.Annotator` to perform statistical comparisons between pairs of conditions (specified in the pairs list) using the chosen test (t-test_ind or Brunner-Munzel). The results are displayed on the boxplot.

## Saving Results:
The script saves the combined image visualizations as `images.svg`.
It saves the boxplot and swarmplot of mitochondrial lengths as `graph.svg`.
