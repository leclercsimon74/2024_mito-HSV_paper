*Please cite this paper if you use the following code!*

The final figure involves merging multiple SVG files and manually reorganizing the network.


**GRO-Seq.py**:

This program performs various data processing and analysis tasks related to gene expression data. Here is an overview of what the program does:

1. Imports necessary libraries: `gzip`, `json`, `pandas`, `networkx`, `matplotlib`, `plotly`, and `numpy`.
2. Sets global parameters such as the value of `fc` and `folder`.
3. Defines several functions:
   - `grab_merge_data`: Extracts specific data from a DataFrame based on gene names or gene synonyms.
   - `hist_from_dict_value`: Generates a histogram from values in a dictionary.
   - `select_junction`: Filters a DataFrame based on a list of protein names.
   - `grab_data`: Extracts specific data from a protein atlas based on gene names or gene synonyms.
   - `remove_link`: Removes links from a networkx graph based on a DataFrame and a specified value.
4. Opens a protein atlas file in JSON format and reads the data into a variable called `data`.
5. Reads gene expression data from an Excel file into a DataFrame called `groseq_data`.
6. Filters the `groseq_data` DataFrame based on certain conditions and assigns the filtered data to variables `up_4h`, `down_4h`, `up_8h`, and `down_8h`.
7. Defines a list of data fields to extract from the protein atlas.
8. Calls the `grab_merge_data` function to extract information for the up- and down-regulated proteins at 4h and 8h time points.
9. Combines the extracted information for up- and down-regulated proteins into dictionaries (`up_4h_info`, `down_4h_info`, `up_8h_info`, and `down_8h_info`).
10. Combines the information for all regulated proteins into dictionaries (`eight_h_info`, `for_h_info`, and `all_info`).
11. Creates dictionaries (`unique_4h_info` and `unique_8h_info`) containing information for uniquely regulated proteins at 4h and 8h time points.
12. Stores the combined information in a DataFrame called `data_df` and saves it as a pickle file.
13. Generates histogram figures and saves them as images.
14. Create a Venn diagram to compare gene expression at different logFC values and saves it as an image.
15. Performs additional processing steps involving protein name mapping.

Require files (can be downloaded):
- 9606.protein.info.v11.5.txt.gz
- 9606.protein.links.full.v11.5.txt.gz
- proteinatlas.json.gz

Gro-Seq dataset:
- getDiffExpression_VERO_ChlSab1.1_noadj_rpkmAdded_GeneInfoAdded_EXCEL2.xlsx
