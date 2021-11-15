### Disclaimer

This tool should be considered a hypothesis generating tool with the following caveats:
* **2021-11**: Color by identity for the amplicon protocol includes missing positions, inflating the number of unique barcodes found in the population. Other methods for clustering barcodes are under investiagtion.  
* **2021-08**: Branch lengths have not been adjusted for estimated GW mutation rates. Do not consider the date on the tree options as an accurate representation of divergence time between samples.


### Coloring options

**Gene_Protocol**: Barcode identity using 4 mitochondrial genes (cytB, CO3, ND3-5) and the amplicon protocol for samples that overlap. For more information on the previous barcode groups, see [Ribado et al. 2021, PLOS NTDs](https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0009609).

**Amplicon_Protocol**: Barcode identity using the new amplicon protocol targeting 80% of the mitochondrial genome. Currently colored by exact identity of all variant sites after filtering.

**Case_GPS**: The provided coorindated for each worm, when available. Defaults to country GPS coordinates when case coorindates were not provided. Cases with overlapping coorindates will be shown as a pie chart for the user defined coloring option. 


### VCF Information

Specific information on the underlying dataset and filtering parameters. The number of variants at the end of filtering are used to determine 'Amplicon_Protocol' identity. 

Summary|Value
-----|-----
