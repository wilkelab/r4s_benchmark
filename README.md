#Calculating site-specific evolutionary rates at the amino-acid or codon level yields similar rate estimates

This repository contains all the scripts and data to reproduce the results of:

D. K. Sydykova, C. O. Wilke (preprint). Calculating site-specific evolutionary rates at the amino-acid or codon level yields similar rate estimates. PeerJ Preprints 5:e2739v1. [https://peerj.com/preprints/2739v1/] (https://peerj.com/preprints/2739v1/)

##Contents:

* `mech_codon` contains results for the alignments simulated with the dN/dS model.

	+ `assigned_rates` contains true site-wise dN/dS.
	
	+ `filtered_sites` contains information on all sites without any amino acid substitutions for each simulated alignment. 
	
	+ `inferred_rates` contains inferred site-wise dN/dS.
	
	+ `processed_rates` contains tables with all site-wise rates: true dN/dS, inferred dN/dS, and inferred Rate4Site.
	
	+ `r4s_rates` contains inferred site-wise Rate4site rates. 
	 
* `mut_sel` contains results for the alignments simulated with the mutation-selection (MutSel) model. MutSel alignments were simulated by Spielman et al. (2016). True site-wise dN/dS for their alignments can be found in their repository [https://github.com/sjspielman/dnds_1rate_2rate] (https://github.com/sjspielman/dnds_1rate_2rate)

	+ `filtered_sites` contains information on all sites without any amino acid substitutions for each simulated alignment. 
	
	+ `processed_rates` contains tables with all site-wise rates: true dN/dS, inferred dN/dS, and inferred Rate4Site.
	
	+ `r4s_rates` contains inferred site-wise Rate4site rates. 
	
* `natural_prot` contains data and results for the natural alignments from Spielman and Wilke (2013) and Meyer and Wilke (2015). The data we used can be found at [https://github.com/sjspielman/mammalian_gpcr_selection](https://github.com/sjspielman/mammalian_gpcr_selection) and [https://github.com/ausmeyer/hiv_structural_determinants](https://github.com/ausmeyer/hiv_structural_determinants), respectively. 

	+ `aln` contains HIV-1 and GPCR sequences used in our analysis
	
		+ `aligned_seqs` contains GPCR amino acid sequences we aligned.
		
		+ `back_translated_aln` codon alignments that were translated back from amino acid alignments.
		
		+ `raw_aln` contains raw FASTA files from the repositories mentioned.
		
		+ `reforematted_aln` contains nucleotide alignments with sequence IDs reformatted. These were used as input for `HyPhy`. 
	
	+ `filtered_sites` contains information on all sites without any amino acid substitutions for each alignment. 

	+ `inferred_dNdS` contains site-wise inferred dN/dS.
	
	+ `processed_rates` contains tables with site-wise inferred dN/dS and inferred Rate4Site. 
	
	+ `r4s_rates` contains inferred site-wise Rate4site rates. 
	
	+ `trees` contains trees inferred from amino acid alignment for each protein. This directory also contains trees with reformatted sequence IDs to be used as input for `HyPhy`.
	
* `plots` contains final figures used in the publication.

* `src` contains all of the scripts used to analyse the data and plot the figures. The usage of each script is described in the analysis section. 

##Analysis:
	
###dN/dS model

###MutSel model

###Natural proteins
	