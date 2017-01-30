#Calculating site-specific evolutionary rates at the amino-acid or codon level yields similar rate estimates

This repository contains all the scripts and data to reproduce the results of:

D. K. Sydykova, C. O. Wilke (preprint). Calculating site-specific evolutionary rates at the amino-acid or codon level yields similar rate estimates. PeerJ Preprints 5:e2739v1. doi: [10.7287/peerj.preprints.2739v1] (10.7287/peerj.preprints.2739v1)

##Contents:

`mech_codon` contains site-wise rates for the alignments simulated with the dN/dS model. This directory also contains information on sites in the simulated alignments with no substitutions. 

This directory contains:
	1. `assigned_rates` contains true site-wise dN/dS for all simulated alignments.
	2. `filtered_sites` contains information on all sites with no substitutions for each simulated alignment. 
	3. `inferred_rates` contains inferred site-wise dN/dS for all simulated alignments.
	4. `processed_rates` contains tables with all site-wise rates (true dN/dS, inferred dN/dS, and inferred Rate4Site) for all simulated alignments.
	5. `r4s_rates` contains inferred site-wise Rate4site rates for all simulated alignments. 
	 
`mut_sel` contains site-wise rates for the alignments simulated with the mutation-selection (MutSel) model. MutSel alignments were simulated by Spielman et al. (2016). True site-wise dN/dS for their alignments can be found in their repository [https://github.com/sjspielman/dnds_1rate_2rate] (https://github.com/sjspielman/dnds_1rate_2rate)

This directory contains:
	1. `filtered_sites` contains information on all sites with no substitutions for each simulated alignment. 
	4. `processed_rates` contains tables with all site-wise rates (true dN/dS, inferred dN/dS, and inferred Rate4Site) for all simulated alignments.
	5. `r4s_rates` contains inferred site-wise Rate4site rates for all simulated alignments. 
	
`natural_prot`

`plots`

`src`