> [!IMPORTANT]
> Please cite the paper if you use the data or the code!

![Image here](https://github.com/leclercsimon74/2024_mito-HSV_paper/blob/main/Figure%204/Figure4_LD.png)

**Fig 4. The number of ER-mitochondria contact sites increase and they are clustered in infection.**
**(a)** Visualization of the mitochondrial structure by the ten-fold robust expansion microscopy (TREx) in noninfected and infected Vero cells at 8 hpi (n = 6). The distributions of an ER protein tyrosine vesicle-associated membrane protein B (VAPB, yellow) located in the ER-mitochondria contact sites and mitochondria labeled with MitoTracker (magenta) are shown. The localization of the nuclear viral replication compartment and the cytoplasmic viral ICP4 protein are presented by HSV-1 EYFP-ICP4 (cyan) and the nucleus by the DAPI stain (blue). Scale bars, 1 μm. **(b)** Violin plots show the distance between the mitochondria and VAPB. **(c)** Proximity ligation analysis (PLA) of contact sites in noninfected and infected Vero cells at 4, 8, and 12 hpi. The PLA signal between VAPB and regulator of microtubule dynamics (RMDN3) is visualized by fluorescent spots (yellow). Mitochondria are labeled with MitoTracker (red), the nucleus with DAPI, and viral replication compartments are visualized by EYFP-ICP4. Scale bars, 10 μm. **(d)** Box plots showing the number of PLA foci and **(e)** volume of the foci per cell (ncells = 27, 27, 28, and 24 for NI, 4, 8, and 12 hpi, respectively). The box plots show the mean (dashed line) and the interquartile range. Statistical significance was determined using the Student´s t-test. The significance values are denoted as **** (p<0.0001), ** (p<0.01), * (p<0.05) or ns (not significant).

Link to the publication: https://doi.org/10.1371/journal.ppat.1011829.g004

## [Graph_making.py](https://github.com/leclercsimon74/2024_mito-HSV_paper/blob/main/Figure%204/Graph_making.py)
### General Description:
This script performs data analysis and visualization of the distance between VAPB foci and mitochondria and the characteristics of PLA (Proximity Ligation Assay) foci across different infection conditions. The script removes outliers, creates violin plots, swarm plots, and boxplots, and applies statistical tests to compare conditions.

### Detailed Description:

1. **Imports and Initial Configuration**:
   - The script imports necessary libraries for data manipulation (`numpy`, `pandas`), statistical analysis (`scipy`), and visualization (`matplotlib`, `seaborn`).
   - It also imports a custom utility module (`utils`) and the `Annotator` class for statistical annotations on plots.
   - Sets the default mathematical text rendering in matplotlib to regular font.

2. **Utility Functions**:
   - `outlier_removal(data, zscore)`: Removes outliers from the data using the z-score method.
   - `make_swarmplot(data, ax, dotsize, pairs)`: Creates a swarm plot with overlaid boxplot and statistical annotations.

3. **Settings and Parameters**:
   - Sets the size of the figure and various parameters for data filtering and visualization.
   - Configures the color palette for the plots using a custom utility function `u.get_color`.

4. **Data Loading and Filtering**:
   - Loads the distance analysis data for VAPB foci (`ExpM/Foci_distance_analysis.pckl`) and filters the data based on area size.
   - Defines pairs of conditions for statistical comparison and the statistical test to be used (`Brunner-Munzel`).
   - Creates a dictionary to store the distances of VAPB foci for non-infected and 8 hours post-infection conditions.
   - Converts the dictionary to a pandas DataFrame and saves it as a CSV file.

5. **Visualization - VAPB Foci to Mitochondria Distance**:
   - Creates a violin plot and a boxplot to visualize the distances of VAPB foci to mitochondria.
   - Adds a horizontal line indicating the contact site threshold and applies statistical annotations to compare conditions.
   - Sets the y-axis scale to a custom log transformation for better visualization.

6. **Data Loading and Filtering - PLA Foci**:
   - Loads the distance analysis data for PLA foci (`PLA/PLA_foci_distance.pckl`) and filters the data based on foci volume.
   - Defines pairs of conditions for statistical comparison and sets the z-score threshold for outlier removal.

7. **Visualization - Number of PLA Foci per Cell**:
   - Groups the filtered data by cell name and counts the number of PLA foci per cell for each condition.
   - Removes outliers and converts the data to a pandas DataFrame.
   - Creates a swarm plot with overlaid boxplot to visualize the number of PLA foci per cell and applies statistical annotations.

8. **Visualization - Mean Distance of PLA Foci to Mitochondria per Cell**:
   - Groups the filtered data by cell name and calculates the mean distance of PLA foci to mitochondria for each condition.
   - Removes outliers and converts the data to a pandas DataFrame.
   - Creates a swarm plot with overlaid boxplot to visualize the mean distance of PLA foci to mitochondria per cell and applies statistical annotations.

9. **Visualization - Mean Volume of PLA Foci per Cell**:
   - Groups the filtered data by cell name and calculates the mean volume of PLA foci for each condition.
   - Removes outliers and converts the data to a pandas DataFrame.
   - Creates a swarm plot with overlaid boxplot to visualize the mean volume of PLA foci per cell and applies statistical annotations.

10. **Saving Figures**:
    - Adjusts the layout of the plots for better visualization.
    - Saves the figures in both PNG and SVG formats.


### Key Outputs:
- **VAPB Foci to Mitochondria Distance Plot**: Visualizes the distances of VAPB foci to mitochondria for non-infected and 8 hours post-infection conditions, with a contact site threshold line.
- **Number of PLA Foci per Cell Plot**: Visualizes the number of PLA foci per cell for different infection conditions.
- **Mean Distance of PLA Foci to Mitochondria per Cell Plot**: Visualizes the mean distance of PLA foci to mitochondria per cell for different infection conditions.
- **Mean Volume of PLA Foci per Cell Plot**: Visualizes the mean volume of PLA foci per cell for different infection conditions.
- **CSV Files**: Saves the processed data used for plotting into CSV files.

Graphs generated (on the side of Figure 4):
![Graphs](https://github.com/leclercsimon74/2024_mito-HSV_paper/blob/main/Figure%204/Figure-graph.png)
