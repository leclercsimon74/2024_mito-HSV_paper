![Figure S2](https://github.com/leclercsimon74/2024_mito-HSV_paper/blob/main/Supplementary%20Figure%20S2%20(GRO-Seq%20on%20Hela)/Supplementary%20Figure%20S2.png)

**S2 Fig. Analysis of nascent RNA levels of mitochondrial and mitochondria-associated proteins in infected Hela cells.**
**(a)** A Venn comparison of GRO-Seq dataset showing gene regulation of mitochondrial protein-ssociated gene expression in infected Hela cells at 4 and 8 hpi. The number of regulated proteins at 4 (green) and 8 hpi (pink), or at both time points (blue) is shown. The size of the circle is proportional to the number of regulated proteins. **(b)** The upregulation and downregulation of mitochondrial genes at 4 and 8 hpi clustered according to their predicted functional biological process during infection. The horizontal line at value 0 represents gene expression in noninfected cells, and positive and negative values represent the upregulation and downregulation of gene expression, respectively.

https://doi.org/10.1371/journal.ppat.1011829.s002

# make_table.py
Require the [protein atlas](https://www.proteinatlas.org/download/proteinatlas.json.gz)

## Libraries:

- `pandas`: used for data manipulation and analysis (creating DataFrames, sorting, etc.)
- `numpy`: used for numerical computations (arrays, mathematical functions)
- `matplotlib`: used for creating visualizations (heatmaps, bar charts, etc.)
- `json`: used for working with JSON data (loading and parsing)
- `gzip`: used for working with gzip compressed files (reading proteinatlas.json.gz)
- `os`: used for interacting with the operating system (file paths)


## Data Loading and Preprocessing:

### Load Protein Data:
The code reads protein data from a JSON file named `proteinatlas.json.gz` using gzip.open and json.loads.
### Load Gene Expression Data:
It reads gene expression data from an Excel file named `"GRO-seq HSV1"+os.sep+"getDiffExpression_JKL_Hela112018_HSV1_0_12h_rpkm_added_condensed_EXCEL_sorting.xlsx"` using pandas.read_excel.
It adds a "Gene name" column by splitting the "Annotation/Divergence" column.
## Gene Selection and Filtering:
It retrieves a list of gene names from the protein data using a list comprehension.

The code filters the gene expression data based on fold change (FC) and p-value thresholds. It selects genes with:
- FC > fc (user-defined threshold, set to 0.585, 1, or 2 in the code)
- p-value < 0.05
- Maximum RPKM > 0.5

## Extracting Information:
A list named to_extract specifies the information to be extracted for each gene from the protein data: "Gene", "Molecular function", "Biological process", "Subcellular location", and "Uniprot".

Function `grab_merge_data`: This function takes gene expression data (filtered DataFrame) and a list of gene names as input.
- It iterates through each gene in the expression data.
- It searches for the corresponding gene in the protein data (considering synonyms).
- If a match is found, it extracts the specified information and stores it in a dictionary.
- The function returns a dictionary containing the extracted information for each identified gene.

## Data Analysis and Visualization:
The code creates separate DataFrames for up-regulated and down-regulated genes at each time point (4h, 8h, 12h).

It combines the up-regulated and down-regulated genes from all time points into a single DataFrame.

The code focuses on genes located in the mitochondria (`Mitochondria`). It counts the number of up- and down-regulated genes in the mitochondria for each time point.
### Venn Diagrams:
The code uses matplotlib_venn to create Venn diagrams to visualize the overlap between genes regulated at different time points (full dataset and mitochondria focus).

### Dual Bar Plots (Mito Only):
The code analyzes the biological processes of genes located in the mitochondria.
It categorizes genes by their biological process and creates dual bar charts to show the count of up-regulated and down-regulated genes for each process at different time points (4h and 8h).
### Similar Analysis on Full Dataset:
The code performs a similar analysis on the full dataset (not limited to mitochondria) to identify biological processes affected by the treatment.
## Saving Results:

The code saves the processed data (full gene information) as a CSV file.
It saves the results of the subcellular localization analysis (counts and tables) as CSV files.
It saves the generated visualizations (Venn diagrams and bar charts) as PNG and SVG images.
