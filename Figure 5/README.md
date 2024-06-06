> [!IMPORTANT]
> Please cite the paper if you use the data or the code!

![Image here](https://github.com/leclercsimon74/2024_mito-HSV_paper/blob/main/Figure%205/Figure5%20-%20Manual.png)

**Fig 5. The cristae become thicker and shorter along the progression of the infection.**

**(a)** Representative focused ion-beam scanning electron microscopy (FIB-SEM) images of noninfected and HSV-1-infected MEF cells at 8 hpi. The mitochondrial outer membrane (yellow) and cristae (green) are shown. Scale bars, 0.5 μm. **(b)** The 3D structure of cristae (green) reconstructed from FIB-SEM stacks. See also S4 Movie. The 3D quantitative analysis of **(c)** the maximal thickness and **(d)** the length of cristae calculated using a watershed algorithm to individualize the cristae lamella in noninfected and infected cells (nmito = 8, ncrist = 112 and 121 for NI and 8 hpi, respectively). **(e)** The surface area of segmented cristae in each mitochondria. The box plots show the mean (dashed line) and the interquartile range. Statistical significance was determined using the Student´s t-test. The significance values are denoted as **** (p<0.0001) or ns (not significant).

https://doi.org/10.1371/journal.ppat.1011829.g005


## [Analysis+Figure.py](https://github.com/leclercsimon74/2024_mito-HSV_paper/blob/main/Figure%205/Analysis%2BFigure.py)

The Python script performs a comparative analysis of mitochondrial cristae between non-infected (NI) and infected (8 hours post-infection, hpi) cells. The script processes 3D segmentation images of mitochondria and cristae to compute various metrics and visualize the results. Here is a detailed explanation of the script, including its functionalities and how it operates.

### Script Breakdown

1. **Imports and Initial Setup**:
   - The script imports various libraries such as `tifffile`, `matplotlib`, `seaborn`, `numpy`, `scipy`, and others for image processing, analysis, and visualization.
   - `calculate_local_thickness`, `measure2D`, `angle_between_v1_v2`, `found_intersection`, `get_tang_mask` functions are defined to aid in image analysis.
   
2. **Function Definitions**:
   - `calculate_local_thickness`: Calculates the local thickness of a 3D binary array.
   - `measure2D`: Measures properties of labeled 2D binary images.
   - `angle_between_v1_v2`: Computes the angle between two vectors.
   - `found_intersection`: Finds the intersection point where a vector exits a binary mask.
   - `get_tang_mask`: Retrieves the tangent mask of an intersection point on the edge of a binary image.
   - `extract_data`: Extracts and calculates metrics such as thickness, orientation, surface area, and length of cristae from the given segmented image files.

3. **Data Extraction and Analysis**:
   - Paths to the segmented images of NI and 8 hpi samples are specified.
   - The script iterates through these image paths, applying the `extract_data` function to compute various metrics for each image.
   - Collected metrics include `ratio_volume` (ratio of cristae volume to mitochondria volume), `thickness`, `orientation`, `surface_area`, and `length`.

4. **Visualization**:
   - The script creates a figure comprising several subplots:
     - **Part A**: Shows raw data images with segmentation overlays for both NI and 8 hpi samples.
     - **Part B**: Displays 3D reconstructions of the mitochondria and cristae.
     - **Part C**: Visualizes the computed metrics (thickness, surface area, and length) using swarm and box plots.
   - The `make_swarmplot` function is used to create the swarm and box plots with annotations for statistical significance.

5. **Final Steps**:
   - The script saves the generated figure in both PNG and SVG formats.
   - CSV files for each metric are also saved for further analysis.

