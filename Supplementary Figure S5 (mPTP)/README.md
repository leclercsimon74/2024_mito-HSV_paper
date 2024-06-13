![Figure S5]()

**S5 Fig. Mitochondrial permeability transition pores prefer closed state in infected cells.**
The state of mPTP was assessed by cellular loading of fluorescent Calcein-AM and Co2+ quencher. **(a)** Representative cells showing Calcein (cyan) distribution in the cytosol of noninfected and infected MEF cells at 4, 8, and 12 hpi. The plasma membrane localization is shown (cell, grey line). **(b)** The fluorescent intensity of cytoplasmic Calcein and nuclear viral replication compartment marker, ICP4, in noninfected and infected cells at 4, 8, and 12 hpi (n = 150, 134, 65, and 152, respectively). Statistical significance was determined using the Student´s t-test. The significance values are denoted as **** (p<0.0001), *** (p<0.001), or ns (not significant).

https://doi.org/10.1371/journal.ppat.1011829.s005

# [Data_extraction.py]()
The Python code performs image analysis to extract data on mPTP (Mitochondrial Permeability Transition Pore) experiment.

## Image Loading and Preprocessing
The code uses` tifffile.imread` to read images from two channels (likely representing ICP4 and Mitotracker) for conditions like '4hpi', '8hpi', and '12hpi'.
It performs smoothing on both images using a Gaussian filter (`scipy.ndimage.gaussian_filter`) to reduce noise.
## Nuclei Detection
The code uses the Mitotracker channel to identify potential nuclei locations. It performs local thresholding (`skimage.filters.threshold_local`) followed by binary closing operations to clean the data.
Nuclei are identified and labeled using `skimage.measure.label`, and properties like area, centroid, and solidity are measured using `regionprops_table`.

Small and irregular nuclei are removed based on size and solidity thresholds.
## Cytoplasm Segmentation and Association with Nuclei
The code removes the nuclei area from the Mitotracker channel to create a mask for the cytoplasm.
It uses watershed segmentation (`skimage.segmentation.watershed`) to potentially split loosely connected cells into individual ones.
To associate nuclei with their corresponding cytoplasm, the code calculates the distance from each nucleus centroid and uses it along with a marker mask to define cell boundaries.

## Data Extraction and DataFrame Creation
The code iterates through the identified nuclei and their associated cytoplasm labels.
It calculates intensity values (e.g., 95th percentile) from the ICP4 channel within the nucleus region to assess infection levels.
For each cell, it measures the average intensity from the Mitotracker channel within the cytoplasm to quantify mitochondria presence.

The code calculates the difference in average intensity between the two images (potentially representing mPTP activity) within the same cytoplasm region.
Finally, a DataFrame (`pandas.DataFrame`) is created to store these extracted data points for each identified cell, including:
- Nucleus ID
- Nucleus area and coordinates
- Cytoplasm area and coordinates
- Infection intensity
- Cytoplasm Mitotracker intensity
- Average mPTP intensity before and after treatment
- Difference in mPTP intensity between treatments

## Visualization and Saving:
The code generates visualizations using `matplotlib.pyplot` to show the original images, detected nuclei, cell segmentation, and relationships between extracted data points.
It saves the DataFrame for each condition as a pickle file (.pckl) for further analysis.


*Note:*
The code handles scenarios where a nucleus might not have a clear cytoplasm association.
It iteratively refines cell segmentation using watershed if loose cells are detected.
The processing for NI (Non-Infected) data follows a similar approach but uses a different image path.



# [Figure.py]()

This code analyzes microscopy images from an mPTP (Mitochondrial Permeability Transition Pore) experiment to quantify infection levels and mPTP activity in cells.

## Libraries and Global Variables
Imports necessary libraries like `tifffile`, `pandas`, `matplotlib`, `scipy`, etc., for image processing, data manipulation, and visualization.
Defines global variables like `resize_f` for image resizing, `pxl_size` for scaling the displayed image, `infection_threshold` to identify infected cells, `condition_name_list` for experiment conditions, and color palettes.

## Functions
`make_swarmplot`: Creates a swarmplot with statistical tests for comparing data distributions across different conditions.

`add_scalebar`: Adds a scale bar to the plot based on pixel size and scale bar size in the image.

`pick_cell`: Selects a representative cell based on various intensity and area values from the DataFrame.

`make_cell_img`: Creates a visualization of a single cell with its nucleus, cytoplasm, and intensity signals.


## Looping through Conditions
The code iterates through four conditions: 'NI' (Non-Infected), '4hpi', '8hpi', and '12hpi'.

## Loading Images and Data
For each condition, it loads the corresponding multi-channel microscopy image using `tifffile.imread`.
It reads the pre-processed data (likely containing features like nucleus area, cytoplasm area, infection intensity, etc.) from a pickle file using `pandas.read_pickle`.

## Selecting a Representative Cell
The `pick_cell` function is used to choose a representative cell based on its nucleus and cytoplasm area, infection intensity, and mPTP intensity. This helps visualize an "average" cell for each condition.

## Creating Cell Visualization
The `make_cell_img` function takes the chosen cell data and images as input. It defines a bounding box around the cell based on its cytoplasm coordinates. It extracts sub-images for the specific channels (ICP4 and Mitotracker) and resizes them based on the `resize_f` factor.

Labels for cytoplasm and nucleus are created based on the cell coordinates. The function differentiates between infected and non-infected cells for removing background signal (ICP4) from the mPTP intensity. It combines the nucleus, cytoplasm, mPTP intensity, and optional mitochondrial channel using a custom color merger function. Finally, the function displays the image with a scale bar.

## Data Extraction and Statistical Analysis
The extract_data function extracts specific columns (e.g., mPTP intensity difference or infection intensity) from the DataFrame for each condition, excluding outliers using z-scores.It creates a pandas DataFrame with the extracted data for different conditions. The script defines pairs of conditions for statistical comparisons.

It uses `seaborn.swarmplot` to visualize the data distribution across conditions and performs statistical tests (t-test_ind) using the `statannotations` library for pairwise comparisons. The results are saved as CSV files.

## Saving and Displaying Results
The script creates a figure with subplots for visualizing individual cells and the swarmplots with statistical comparisons.
It saves the figure as an SVG file and displays it using `plt.show()`.
