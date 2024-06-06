## [Data_analysis.py](https://github.com/leclercsimon74/2024_mito-HSV_paper/blob/main/Figure%203/Data/Data_analysis.py)

### General Description:
The script analyzes and visualizes the distances between endoplasmic reticulum (ER) and mitochondria in different conditions (non-infected, 8 hours post-infection, and 12 hours post-infection). It also examines the number of contact sites per mitochondria and the area of these contact sites. Statistical analysis is performed to compare these metrics across different conditions.

### Detailed Description:

1. **Imports and Initial Configuration**:
   - The script imports libraries for data manipulation (`numpy`, `pandas`), statistical analysis (`scipy`), visualization (`matplotlib`, `seaborn`), and custom utilities (`utils`).

2. **Utility Functions**:
   - `make_swarmplot(data, ax, dotsize)`: Creates a swarm plot with overlaid boxplot and statistical annotations.
   - `get_data(data, zscore=3)`: Removes outliers from the data using the z-score method and converts it to a pandas DataFrame.

3. **Settings and Parameters**:
   - Sets the size of the dots for the swarm plot and the statistical test to use (`t-test_ind` or `Brunner-Munzel`).
   - Configures the color palette for the plots.

4. **File Paths**:
   - Specifies the file paths for the distance maps and other data for the different conditions (non-infected, 8h post-infection, 12h post-infection).

5. **Data Loading**:
   - Loads the distance data for each condition using `numpy.load`.
   - Optionally, random subsets of the data can be selected for analysis.

6. **Outlier Removal**:
   - Removes outliers from the data using a z-score threshold.

7. **Data Preparation**:
   - Prepares the data for visualization by creating a dictionary with condition names and their respective data.
   - Converts the dictionary to a pandas DataFrame and arranges the columns in a specified order.

8. **Visualization**:
   - **Violin and Box Plots**:
     - Creates a figure with three subplots.
   - **ER-mitochondria distances**
     - Violin plot to show the distribution
     - Adds a horizontal line at 30 nm to indicate the contact site threshold.
     - Applies statistical annotations to compare the conditions.
   - **Number of Contact Sites**:
     - Loads the number of contact sites and mitochondrial surface area data.
     - Normalizes the number of contact sites by the mitochondrial surface area.
     - Creates a swarm plot to visualize the number of contact sites per mitochondria.
   - **Area of Contact Sites**:
     - Loads the contact site area data.
     - Normalizes the data to convert units from nm² to μm².
     - Creates a swarm plot to visualize the area of contact sites.
     - Applies statistical annotations to compare the conditions.

9. **Statistical Analysis**:
   - Compares the conditions using the specified statistical test and annotates the plots with the results.

10. **Saving Figures**:
    - Saves the figures in both PNG and SVG formats.

### Key Outputs:
- **ER-Mitochondria Distance Plot**: Visualizes the distances between ER and mitochondria for the three conditions, with a threshold line for contact sites.
- **Number of Contact Sites Plot**: Visualizes the number of contact sites per mitochondria surface area.
- **Area of Contact Sites Plot**: Visualizes the area of contact sites.
- **CSV Files**: Saves the processed data used for plotting into CSV files.


Graphs being generated with the code:
![graphs](https://github.com/leclercsimon74/2024_mito-HSV_paper/blob/main/Figure%203/Data/ER-mitochondria%20distance.png)
