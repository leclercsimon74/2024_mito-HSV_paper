## [Image analysis](image_analysis.py)

The provided script processes .tif image files to detect nuclei, infected cells, and cytoplasm using different channels in the images. It also visualizes the segmented images and extracts various features into a DataFrame, which is then saved as a pickle file.

Here's an explanation of the key steps and functions used in the script:

### Detailed Breakdown

1. **Imports and Setup:**
    - Import necessary libraries for image processing, analysis, and file handling such as `scipy`, `skimage` or `pandas`

2. **File Selection:**
    - List all .tif files in the current working directory. *Not provided, too large! (140Mo)*

3. **Processing Each Image File:**
    - Loop through each .tif file and read it using `tifffile.imread`.

4. **Nucleus Detection:**
    - The first channel (`img[0]`) is processed to detect nuclei. Steps include histogram equalization, Gaussian smoothing, Otsu thresholding, hole filling, and binary opening. Finally, connected components are labeled.

5. **Infection Detection:**
    - The second channel (`img[1]`) is processed similarly to detect infection markers.

6. **Measure Nucleus vs Infection Marker:**
    - Measure properties of the nuclei regions using `regionprops_table`. If the file name contains 'NI' *Noninfected*, skip infection filtering.

7. **Cytoplasm Detection:**
    - The fourth channel (`img[3]`) is used to detect the cytoplasm using similar steps. Watershed segmentation is applied to associate nuclei with cytoplasm regions.

8. **Visualization:**
    - Create and save a visualization of the segmented nucleus, infection, and labeled cells. *Join in this folder and after this code*

9. **Measure Signal:**
    - Measure various properties (e.g., area, mean intensity) for each nucleus and its corresponding cytoplasm and infection regions. Store these measurements in a DataFrame.

10. **Save Data:**
    - Save the DataFrame as a pickle file for later analysis.


### Example image output

![NI image](https://github.com/leclercsimon74/2024_mito-HSV_paper/blob/main/Figure%206/Biochemical%20microscopy%20assay/NI-M%2BRfig.png)
![8hpi image](https://github.com/leclercsimon74/2024_mito-HSV_paper/blob/main/Figure%206/Biochemical%20microscopy%20assay/8hpi-M+Rfig.png)


## [Data analysis](data_analysis.py)

This script processes and visualizes data from pickle files containing measurements from microscopy images. It calculates several metrics for each condition, generates swarm plots with statistical annotations, and saves the resulting figures. Here's an explanation and breakdown of the key sections of the script:

### Detailed Breakdown

1. **Imports and Setup:**
    - Import necessary libraries for data processing, visualization, and statistical analysis such as `seaborn`, `scipy` or `matplotlib`

2. **File Selection:**
    - List all pickle files in the current working directory.

3. **Parameters and Pairs for Statistical Tests:**
    - Define the dot size for swarm plots, the statistical test to be used, and the pairs of conditions to be compared.

4. **Helper Functions:**

    - **`calculate_data(dic)` Function:**
        - Calculate and correct values for nucleus intensity and other metrics, and remove outliers based on z-score.

    - **`grab_and_show(table, condition, title, ax)` Function:**
        - Extract data for each condition, create swarm plots, and annotate with statistical test results.

    - **`exploratory_figure(dic, cond)` Function:**
        - Generate and save a figure with four subplots showing various metrics for a given condition.

5. **Group Data and Calculate Metrics:**
    - Group the data by condition and file type.

6. **Calculate Corrected Metrics:**
    - Calculate corrected metrics for each condition and file type.

7. **Generate Exploratory Figures:**
    - Generate and save exploratory figures for each condition.

8. **Generate Final Analysis Figure:**
    - Generate and save a final analysis figure comparing different metrics across conditions.

