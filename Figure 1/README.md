> [!IMPORTANT]
> Please cite the paper if you use the data or the code!


The final figure involves merging multiple SVG files and manually reorganizing the network.

![Figure 1, Low definition](https://github.com/leclercsimon74/2024_mito-HSV_paper/blob/main/Figure%201/Figure1_v2_LD.png)

**Fig 1. HSV-1 infection alters the host transcriptome.**

Global run-on sequencing (GRO-seq) analysis of nascent RNA levels of mitochondrial and mitochondria-associated proteins in infected Vero cells at **(a)** 4 and **(b)** 8 hpi. The genes encoding for mitochondrial proteins (square-shaped nodes, bolded titles) and nonmitochondrial cellular interactor proteins (round-shaped nodes) are shown. The upregulated (orange) or downregulated (green) transcripts in response to infection are visible together with unregulated interacting transcripts (grey). The node size correlates with a logarithmic fold change (logFC) of regulation and the thickness of black lines between nodes is proportional to the interaction in the STRING database (https://string-db.org/). The main functions of the proteins are denoted with yellow circles with their size proportional to the number of interactors. **(c)** GO term classification for mitochondrial processes in the infected cells at 4 and 8 hpi. The color bar indicates upregulation (orange-red) or downregulation (green-blue), and low change of gene transcription (vertical stripes). The non-significant enrichment is also shown (*).

https://doi.org/10.1371/journal.ppat.1011829.g001

> [!IMPORTANT]
> Python code to generate (part) of the figure

[GRO-Seq_GH.py](GRO-Seq_GH.py):

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

[make_table.py](make_table.py):

- Data Loading: It loads JSON and Excel data (proteinatlas.json.gz and getDiffExpression_VERO_ChlSab1.1_noadj_rpkmAdded_GeneInfoAdded_EXCEL2.xlsx) containing protein information and gene expression data, respectively.
- Data Processing:
   - Filters gene expression data based on certain criteria like fold change and false discovery rate at different time points (4 hours and 8 hours post-infection).
   - Extracts relevant information about proteins from the loaded JSON file based on gene names and synonyms.
- Analysis:
   - Counts the number of proteins associated with a specific subcellular location (e.g., Mitochondria) that are upregulated or downregulated at different time points.
   - Generates Venn diagrams to compare gene expression between different time points.
   - Generates bar plots to visualize the count of proteins associated with different biological processes, categorized by their regulation status and time points.
   - Analyzes protein functions and biological processes using Gene Ontology (GO) annotations and visualizes them using a Sankey diagram.
   - Creates a heatmap to display the log-fold changes of proteins associated with specific GO terms at different time points.

Overall, this code aims to provide insights into how gene expression changes over time during a specific biological process, with a focus on proteins associated with mitochondria.
