# Conditional frequency spectra

## Description
This is more or less a one-stop shop for every result generated in [Patel et al., 2024](https://www.biorxiv.org/content/10.1101/2024.06.15.599126v1). This Snakefile can generate conditional frequency spectra from empirical data as well as from theoretical models, using either SLiM simulations or fastDTWF, for both generic two-population models and the [Jouganous et al. 2017 out-of-Africa model](https://pubmed.ncbi.nlm.nih.gov/28495960/) (ft. some slight modifications).
### Useful bits of data
* `data/gwas` contains GWAS summary statistics for 106 complex traits analyzed in this paper. The summary statistics were curated by Yuval Simons, Hakhamanesh Mostafavi, and Julie Zhu based on the [Neale lab UK Biobank GWAS](https://www.nealelab.is/uk-biobank)
* `data/distributions` contains all frequency spectra analyzed in this paper if you want to tinker with some distributions
 * Note: I would use distributions generated with [fastDTWF](https://github.com/jeffspence/fastDTWF) over those generated with [SLiM](https://messerlab.org/slim/); fastDTWF numerically computes frequency spectra, which is inherently less noisy than approximating a frequency spectrum from 2000 simulations
### Useful bits of code
* `scripts/probability_dist.py` is really nice if you need to do any kind of manipulations with (discretized) frequency spectra or other probability distributions
* `scripts/empirical_cfs.py` will generate conditional frequency spectra for whatever empirical data you want to analyze

## Requirements
### for analyzing empirical data
* `data/CADD_bestfit` - B-statistics, [downloadable from the Sella lab GitHub](https://github.com/sellalab/HumanLinkedSelectionMaps/blob/master/Bmaps/CADD_bestfit.tar.gz) (note: these are in GRCh37 so everything must be in GRCh37 unless you want to convert the B-maps)
* `data/combined_gwas.txt` - table of GWAS summary statistics with columns SNP (chr:bp:ref:alt), effect, and trait_idx
* `data/snp_ID_conversion_table.txt` - table mapping SNP (chr:bp:ref:alt) to rsID
* `data/1KG_sample_info.txt` - table mapping 1K Genomes sample IDs to population labels
* `data/1KGenomes` - 1K Genomes GRCh37 VCFs for each chromosome
* `data/variation_feature.txt.gz` - table mapping rsIDs to ancestral allele states, [downloadable from Ensembl](https://ftp.ensembl.org/pub/release-113/mysql/homo_sapiens_variation_113_38/)
* `data/freq_WB` - table mapping SNP (chr:bp:ref:alt) to alternate allele frequencies in UK Biobank White British (or whatever your GWAS population is)
### for analyzing simulations
* SLiM simulations need to be run separately (see `scripts/slim`) before running the Snakefile (using Snakemake to spawn 200k+ jobs is just bad practice)