# GSEAvalidator
R package to check correlations in GSEA pathways using MS-based proteomics experimental data.

Requires you to first run the MSstats pipeline up to the dataProcess function. The package then takes the protein level output of the dataProcess function and a JSON GSEA file to determine pathways with highly correlated proteins.

In the future this package will include additional functions to determine if MS proteomics experiments would be a good fit for causal inference.
