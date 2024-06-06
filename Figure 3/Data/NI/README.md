### General Description:
The script processes segmented 3D microscopy images of mitochondria, nucleus, and endoplasmic reticulum (ER). It generates 3D meshes for visualization, calculates distances between mitochondria and ER, identifies contact sites, and produces visualizations and figures to highlight these features.

### Detailed Description:

1. **Imports and Initial Configuration**:
   - The script imports necessary libraries for 3D image processing, visualization, and numerical computations.
   - It sets parameters like voxel spacing *(30, 8, 8 nm, ZXY)* and color codes for different cellular structures.

2. **File Paths and Parameters**:
   - Paths to the segmented 3D images of mitochondria, nucleus, and ER are specified.
   - Parameters like maximum ER distance *(700nm, empirically determined)* and contact site distance *(30 nm)* are defined.

3. **Helper Functions**:
   - `extend(img)`: Pads the image with a border of zeros to avoid edge effects during processing.
   - `add_mesh(img, color)`: Converts a binary image to a 3D mesh using the marching cubes algorithm and applies a specified color.

4. **Reading and Processing Images**:
   - The script reads the segmented images using `tifffile`.
   - It binarizes the images (converts labeled images to boolean masks).
   - The `extend` function is applied to the binary images to prevent edge effects.

5. **Generating 3D Meshes**:
   - It generates 3D meshes for the nucleus and ER and adds them to a list of meshes.
   - For mitochondria, it calculates the distance transform to the nearest ER voxel and identifies contact sites.

6. **Distance Calculation and Contact Site Identification**:
   - The script calculates distances from mitochondria to the nearest ER using `distance_transform_edt`.
   - It colors the vertices of the mitochondrial mesh based on their distance to the ER. Contact sites are highlighted in red.

7. **Visualization**:
   - The meshes are visualized using `open3d.visualization.draw_geometries`.
   - 2D slices of the original images are displayed with contours indicating the mitochondria, nucleus, and ER.

8. **Contour and Scalebar Visualization**:
   - It generates figures showing 2D slices of the original image with contours for mitochondria, nucleus, and ER.
   - A scalebar is added to the figure for reference.

9. **Colorbar Creation**:
   - A colorbar is created to show the distance scale used in the visualization.

10. **Contact Site and Surface Area Calculation**:
    - The script labels individual mitochondria and calculates their surface area.
    - It identifies and counts contact sites between mitochondria and ER.
    - The surface area of these contact sites is calculated.

11. **Saving Results**:
    - Lists of contact site areas, numbers, and mitochondrial surface areas are saved as numpy arrays.


### Key Outputs:
- 3D visualizations of the nucleus, ER, and mitochondria with distance-based coloring.
- 2D slice figures showing original images with structure contours.
- Colorbar for the distance scale.
- Arrays of contact site areas, counts, and mitochondrial surface areas.

Output example:
![Figure 3a](https://github.com/leclercsimon74/2024_mito-HSV_paper/blob/main/Figure%203/Data/NI/NI_3DView.png)
