![Supplementary Figure S1](https://github.com/leclercsimon74/2024_mito-HSV_paper/blob/main/Supplementary%20Figure%201%20(GRO-SEQ)/FigureS1.png)

**S1 Fig. Cellular and mitochondrial proteins regulated in response to viral infection.**
**(a)** A Venn diagram of the GRO-Seq dataset showing gene regulation of mitochondrial protein-associated gene transcription in infected Vero cells. The number of regulated proteins at 4 (green) and 8 hpi (pink), or at both time points (blue) is shown. The size of the circle is proportional to the number of regulated proteins. **(b)** The upregulation and downregulation of mitochondrial genes at 4 and 8 hpi clustered according to their predicted functional protein-protein interactions. GO terms of their functional classes of major biological processes detected in infection are shown. The horizontal line at value 0 represents gene transcription in noninfected cells, and positive and negative values represent the upregulation and downregulation of the gene transcription.

Published on: https://doi.org/10.1371/journal.ppat.1011829.s001



[FigureS1.py](FigureS1.py) <- Summarize by ChatGPT:
- Imports necessary libraries/modules such as utils, gzip, json, pandas, and numpy.
- Need to download on [String](https://string-db.org/cgi/download.pl) the protein info (9606.protein.info.v11.5.txt.gz)
- Need to download on [protein atlas](https://www.proteinatlas.org/download/proteinatlas.json.gz) the protein atlas as a json file (proteinatlas.json.gz)
- The script reads in data from an Excel file (getDiffExpression_VERO_ChlSab1.1_noadj_rpkmAdded_GeneInfoAdded_EXCEL2.xlsx) containing gene expression data.
- It filters the gene expression data based on certain criteria, such as fold change (fc) and false discovery rate (FDR), to identify upregulated and downregulated genes at 4 and 8 hours post-infection (hpi) with HSV1.
- The script defines a function grab_merge_data to extract and merge relevant information about genes/proteins from the JSON data and the gene expression DataFrame.
- It extracts information about upregulated and downregulated proteins at 4 and 8 hpi, focusing on attributes such as gene name, molecular function, biological process, subcellular location, etc.
- The script performs further analysis, such as counting the number of proteins localized in mitochondria and generating a table summarizing the upregulated and downregulated proteins at 4 and 8 hpi.
- It visualizes the results using Venn diagrams to compare gene expression between 4 and 8 hpi, as well as bar plots to show the biological processes associated with upregulated and downregulated proteins in mitochondria.

Overall, the script provides a comprehensive analysis of gene expression data, with a focus on proteins related to mitochondria and their regulation during HSV1 infection.
