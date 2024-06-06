Due to the size of the images and the dataset, only example images used as illustrative images in Figure are shown.



### Explanation of [image_analysis.py](https://github.com/leclercsimon74/2024_mito-HSV_paper/blob/main/Figure%204/PLA/image_analysis.py)

This script performs analysis on Proximity Ligation Assay (PLA) images to identify and measure PLA foci, including their distances to mitochondria. The analysis includes segmenting nuclei, applying a watershed algorithm to identify PLA foci, and calculating distances from PLA foci to mitochondria. The results are visualized and saved in a DataFrame. Warning, some confusion on my side with PLA and PCA in the code.

#### Detailed Breakdown:

1. **Imports and Initial Configuration**:
   - The script imports necessary libraries for data manipulation (`numpy`, `pandas`), image processing (`scipy`, `skimage`), and visualization (`matplotlib`).
   - Imports a custom utility module (`utils`).
   - Sets the resolution for the image data (`sampling`).

2. **Function Definition: `PLA_analysis`**:
   This function analyzes images within a specified folder (`condition`), performs segmentation and measurement, and saves results in a DataFrame.

   **Parameters**:
   - `condition`: The folder name containing the images to analyze.
   - `result_df`: A DataFrame to store and merge results from the analysis.

   **Returns**:
   - Updated `result_df` with the analysis results.

3. **Image Processing Pipeline**:
   - **Image Loading**: Loads the images from the specified folder using a utility function (`u.get_img`).
   - **Nucleus Segmentation**: Segments nuclei using Otsu's method.
   - **PLA Segmentation**: Smooths the PLA image using a Gaussian filter and segments PLA foci using a triangular threshold method. The segmented nuclei are removed from the PLA segmentation.
   - **Watershed Segmentation**: Applies the watershed algorithm to separate individual PLA foci and labels them.
   - **Region Properties Measurement**: Measures properties of the segmented PLA foci (e.g., area, centroid, intensity).
   - **Mitochondria Segmentation**: Smooths the mitochondria image and segments mitochondria using Otsu's method. Removes nuclei from the segmentation.
   - **Distance Calculation**: Calculates the distance from each PLA foci to the closest mitochondria.

4. **Visualization**:
   - Creates subplots to visualize different stages of the analysis, such as nucleus segmentation, PLA segmentation, watershed results, mitochondria segmentation, and distance calculation.


5. **Main Execution**:
   - Initializes an empty DataFrame to store results.
   - Calls the `PLA_analysis` function for each condition (e.g., 'NI', '4 hpi', '8 hpi', '12 hpi').


### Key Outputs:
- **Visualizations**: Intermediate visualizations of segmentation and distance measurements for each image.
- **CSV and Pickle Files**: `PLA_foci_distance.csv` and `PLA_foci_distance.pckl` containing detailed measurements and analysis results for all conditions.


