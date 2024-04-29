### Disclaimer

The Guinea worm instance of Nexstrain should be considered a hypothesis generating tool with the following caveats:
* **2024-04**: Amplicon barcode colors have been updated to account for more observed barcodes in the specimen library. Barcodes observed in less than 20 specimens have been reduced to limit the number of colors. For further exploration of less common barcodes, use the "Amplicon_Barcode" filtering option to show specific barcode groupings. ß
* **2024-04**: After further investigations to cluster barcodes with missing positions to more common barcodes, barcodes with few missing positions are generally not genetically identical to other barcodes with all variant positions available. To avoid incorrect putative linkages with clustering, researchers recommend considering pairwise relatedness (not highlighted in Nextstrain) while balancing pairwise barcode completeness and epidemiological investigation evidence when interpreting more rarely observed barcodes. 
* **2021-11**: Color by identity for the amplicon protocol includes missing positions, potentially inflating the number of unique barcodes found in the population. Other methods for clustering barcodes are under investigation.  
* **2021-08**: Branch lengths have not been adjusted for estimated GW mutation rates. Estimates may also be biased due to missing collection date information (given 01-01-{collection year} as default when missing). Do not consider the date on the tree options as an accurate representation of divergence time between samples. Use the `Sampled_Year` option for the correct case year per sample.  


### Coloring options

**Amplicon_Protocol**: Barcode identity using the updated tiled amplicon protocol targeting 80% of the mitochondrial genome. Currently colored by exact identity of all variant sites after filtering.

**Gene_Protocol** (Optional): Barcode identity using 4 mitochondrial genes (cytB, CO3, ND3-5). For more information on gene protocol barcode groups, see [Ribado et al. 2021, PLOS NTDs](https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0009609).

**Yearly_Sibship** (Optional): One hundred and thirty nuclear microsatellite sites were tested for informativeness. Analyses with 10 different sets of randomly chosen loci indicated that the presence/absence of some loci influenced the pedigree assignment of some groups. Kinship analysis was repeated after assessing the Polymorphism Information Content (PIC) of each locus using the method of Botstein et al. (1980).
Results in the current analysis were obtained using loci with a PIC >= 0.5, which was considered by Botstein et al. (1980) to be the minimum value for a highly informative marker (i.e.,
a marker with high discriminatory capacity). Sibship inference per year determined with Colony (v2.0.6.7). Genotypes for each sample were coded as present/absent haplotypes at loci with PIC >= 0.5. Loci included in this analysis are a subset of the loci that provided successful allele calls in >= 75% of the 384-sample trial run. Samples included in this category met allele-calling parameters for >= 70% of 130 microsatellite loci. Individuals were assigned to a a cluster if ML probability of inclusion >= 0.90. 

**ParentalOffspring_Pairs** (Optional): Parentage inference determined with Colony (v2.0.6.7). See `Yearly_Sibship` for more details on processing.  

**Case_GPS**: The provided coordinates for each worm, when available. Samples were given country GPS coordinates when case coordinates were unavailable. Cases with overlapping coordinates are shown as a pie chart for the user defined coloring option. 


### VCF Information

Specific information on the underlying dataset and filtering parameters. The number of variants at the end of filtering are used to determine 'Amplicon_Protocol' identity. 

Summary|Value
-----|-----
