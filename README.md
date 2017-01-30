#Calculating site-specific evolutionary rates at the amino-acid or codon level yields similar rate estimates

This repository contains all the scripts and data to reproduce the results of:

D. K. Sydykova, C. O. Wilke (preprint). Calculating site-specific evolutionary rates at the amino-acid or codon level yields similar rate estimates. PeerJ Preprints 5:e2739v1. [https://peerj.com/preprints/2739v1/] (https://peerj.com/preprints/2739v1/)

##Contents:

- `mech_codon` contains results for the alignments simulated with the dN/dS model. 

This directory contains:

	- `assigned_rates` contains true site-wise dN/dS for all simulated alignments.
	- `filtered_sites` contains information on all sites with no substitutions for each simulated alignment. 
	- `inferred_rates` contains inferred site-wise dN/dS for all simulated alignments.
	- `processed_rates` contains tables with all site-wise rates (true dN/dS, inferred dN/dS, and inferred Rate4Site) for all simulated alignments.
	- `r4s_rates` contains inferred site-wise Rate4site rates for all simulated alignments. 
	 
- `mut_sel` contains results for the alignments simulated with the mutation-selection (MutSel) model. MutSel alignments were simulated by Spielman et al. (2016). True site-wise dN/dS for their alignments can be found in their repository [https://github.com/sjspielman/dnds_1rate_2rate] (https://github.com/sjspielman/dnds_1rate_2rate)

This directory contains:

	- `filtered_sites` contains information on all sites with no substitutions for each simulated alignment. 
	- `processed_rates` contains tables with all site-wise rates (true dN/dS, inferred dN/dS, and inferred Rate4Site) for all simulated alignments.
	- `r4s_rates` contains inferred site-wise Rate4site rates for all simulated alignments. 
	
`natural_prot`

`plots`

`src`