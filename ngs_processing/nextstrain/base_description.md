### Disclaimer

This tool should be considered a hypothesis generating tool with the following caveats:
* **2021-11**: Color by identity for the amplicon protocol includes missing positions, inflating the number of unique barcodes found in the population. Other methods for clustering barcodes are under investigation.  
* **2021-08**: Branch lengths have not been adjusted for estimated GW mutation rates. Do not consider the date on the tree options as an accurate representation of divergence time between samples.


### Coloring options

**Gene_Protocol**: Barcode identity using 4 mitochondrial genes (cytB, CO3, ND3-5). For more information on gene protocol barcode groups, see [Ribado et al. 2021, PLOS NTDs](https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0009609).

**Amplicon_Protocol**: Barcode identity using the new amplicon protocol targeting 80% of the mitochondrial genome. Currently colored by exact identity of all variant sites after filtering.

**Yearly_Sibship**: Sibship inference per year determined with Colony (v2.0.6.7). Genotypes for each sample were coded as present/absent haplotypes at 130 nuclear microsatellite loci. Samples included in this category met allele-calling parameters for >= 70% of 130 microsatellite loci.

**ParentalOffspring_Pairs**: Parentage inference determined with Colony (v2.0.6.7). See `Yearly_Sibship` for more details on processing.  

**Case_GPS**: The provided coordinates for each worm, when available. Samples were given country GPS coordinates when case coordinates were unavailable. Cases with overlapping coordinates are shown as a pie chart for the user defined coloring option. 


### VCF Information

Specific information on the underlying dataset and filtering parameters. The number of variants at the end of filtering are used to determine 'Amplicon_Protocol' identity. 

Summary|Value
-----|-----