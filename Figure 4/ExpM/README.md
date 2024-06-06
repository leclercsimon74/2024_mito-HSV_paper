Due to the size of the images and the dataset, only example images used as illustrative images in Figure 4 are shown.


### Explanation of the [Image_processing_foci.py](https://github.com/leclercsimon74/2024_mito-HSV_paper/blob/main/Figure%204/ExpM/Image_processing_foci.py)

This script performs image analysis to study the spatial relationship between VAPB (Vesicle-associated membrane protein-associated protein B) and mitochondria in cells. The analysis includes segmenting nuclei, removing nuclear regions from the images, detecting VAPB foci, calculating distances from these foci to mitochondria, and visualizing the results. The final results are stored in a DataFrame and saved for further analysis.

#### Detailed Breakdown:

1. **Imports and Initial Configuration**:
   - The script imports necessary libraries for image processing (`numpy`, `scipy`, `skimage`), data handling (`pandas`), and visualization (`matplotlib`).
   - Imports a custom utility module (`utils`).
   - Sets pixel size and magnification factor to calculate spacing for image processing.

2. **Variables Initialization**:
   - `pixel_size`: Specifies the pixel size in micrometers for Z, X, and Y dimensions.
   - `mag_f`: Magnification factor.
   - `spacing`: Calculates the spacing based on pixel size and magnification factor.
   - `NI_path_list` and `INF_path_list`: Get the list of image file paths for the non-infected (NI) and infected (8 hpi) conditions using a utility function.

3. **DataFrame Initialization**:
   - `final_df`: An empty DataFrame to store the final results.

4. **Processing Non-Infected (NI) Images**:
   - Iterates over each image file in `NI_path_list`.
   - **Image Loading**: Loads the image data in ZCXY format using a utility function.
   - **Image Splitting**: Splits the loaded image into separate channels for nuclei, mitochondria, and VAPB.
   - **Nucleus Segmentation**: Smooths the nucleus image using a Gaussian filter, determines a threshold using Otsu's method, and performs binary operations to clean and fill the segmented nuclei.
   - **Nucleus Removal**: Removes the nucleus regions from the mitochondria and VAPB images.
   - **Smoothing**: Applies Gaussian smoothing to the mitochondria and VAPB images.
   - **VAPB Segmentation**: Thresholds the VAPB image to create a binary mask and labels the VAPB foci.
   - **Foci Properties Measurement**: Measures properties of the labeled VAPB foci (area, centroid, intensity) using `skimage.measure.regionprops_table`.
   - **Distance Calculation**: Calculates the distance from each VAPB foci to the closest mitochondria using Euclidean distance transform.
   - **Result Storage**: Adds the results to `final_df` with additional information about the condition and file name.
   - **Visualization**: Creates subplots to visualize the labeled VAPB foci, dilated VAPB foci, and the distance histogram.

5. **Processing Infected (8 hpi) Images**:
   - Similar steps are repeated for the infected images in `INF_path_list` with slight variations in binary dilation and erosion iterations for nucleus segmentation.

6. **Saving Results**:
   - Saves the final DataFrame `final_df` to CSV and pickle files for later analysis.

### Key Outputs:
- **Visualizations**: Intermediate visualizations of labeled VAPB foci, dilated VAPB foci, and distance histograms for each image.
- **CSV and Pickle Files**: `Foci_distance_analysis.csv` and `Foci_distance_analysis.pckl` containing detailed measurements and analysis results for non-infected and infected conditions.

