#Calculating site-specific evolutionary rates at the amino-acid or codon level yields similar rate estimates

This repository contains all the scripts and data to reproduce the results of:

D. K. Sydykova, C. O. Wilke (preprint). Calculating site-specific evolutionary rates at the amino-acid or codon level yields similar rate estimates. PeerJ Preprints 5:e2739v1. [https://peerj.com/preprints/2739v1/] (https://peerj.com/preprints/2739v1/)

##Contents

`mech_codon` contains results for the alignments simulated with the dN/dS model.

* `assigned_rates` contains true site-wise dN/dS.
	
* `filtered_sites` contains information on all sites without any amino acid substitutions for each simulated alignment. 
	
* `inferred_rates` contains inferred site-wise dN/dS.
	
* `processed_rates` contains tables with all site-wise rates: true dN/dS, inferred dN/dS, and inferred Rate4Site.
	
* `r4s_rates` contains inferred site-wise Rate4site rates. 
	 
`mut_sel` contains results for the alignments simulated with the mutation-selection (MutSel) model. MutSel alignments were simulated by Spielman et al. (2016). True site-wise and inferred dN/dS for their alignments can be found in their repository [https://github.com/sjspielman/dnds_1rate_2rate] (https://github.com/sjspielman/dnds_1rate_2rate)

* `filtered_sites` contains information on all sites without any amino acid substitutions for each simulated alignment. 
	
* `processed_rates` contains tables with all site-wise rates: true dN/dS, inferred dN/dS, and inferred Rate4Site.
	
* `r4s_rates` contains inferred site-wise Rate4site rates. 
	
 `natural_prot` contains results for the natural alignments from Spielman and Wilke (2013) and Meyer and Wilke (2015). The data we used can be found at [https://github.com/sjspielman/mammalian_gpcr_selection](https://github.com/sjspielman/mammalian_gpcr_selection) and [https://github.com/ausmeyer/hiv_structural_determinants](https://github.com/ausmeyer/hiv_structural_determinants), respectively. 

* `aln` contains HIV-1 and GPCR sequences used in our analysis

	+ `aligned_seqs` contains GPCR amino acid sequences we aligned.
	
	+ `back_translated_aln` codon alignments that were translated back from amino acid alignments.
	
	+ `raw_aln` contains raw FASTA files from the repositories mentioned.
	
	+ `reforematted_aln` contains nucleotide alignments with sequence IDs reformatted. These were used as input for `HyPhy`. 
	
* `filtered_sites` contains information on all sites without any amino acid substitutions for each alignment. 

* `inferred_dNdS` contains site-wise inferred dN/dS.

* `processed_rates` contains tables with site-wise inferred dN/dS and inferred Rate4Site. 

* `r4s_rates` contains inferred site-wise Rate4site rates. 

* `trees` contains trees inferred from amino acid alignment for each protein. This directory also contains trees with reformatted sequence IDs to be used as input for `HyPhy`.
	
`plots` contains final figures used in the publication.

`src` contains all of the scripts used to analyze the data and plot the figures. The usage of each script is described in the analysis section below. 

##Analysis
	
###dN/dS model

The analysis in this section requires [https://github.com/sjspielman/dnds_1rate_2rate] (https://github.com/sjspielman/dnds_1rate_2rate) in the same directory as the current repository.

1. Copy trees from [https://github.com/sjspielman/dnds_1rate_2rate] (https://github.com/sjspielman/dnds_1rate_2rate) using the command line `cp ../dnds_1rate_2rate/trees/n*_bl*.tre ./trees/`.

2. Simulate alignments using `./src/write_run_sim_aln.sh`. This script will write `run_sim_aln.sh` to simulate dN/dS alignments. 

3. Translate simulated nucleotide alignments to amino acids using `./src/write_run_translate_aln.sh`.

4. Infer site-wise dN/dS with `HyPhy` using the script in `dnds_inference/submit_run_inference.sh` directory in `src`. This script was copied from [https://github.com/sjspielman/dnds_1rate_2rate] (https://github.com/sjspielman/dnds_1rate_2rate) and modified for this analysis.

5. Infer site-wise Rate4Site scores using `./src/write_run_r4s_mech_codon.sh`. This script will write `run_r4s_mech_codon.sh` which uses `r4s_pipeline.sh` to run Rate4Site on simulated alignments. 

6. Concatenate all rates into a table with `./src/concatenate_mech_codon_rates.r`. 

###MutSel model

The analysis in this section requires [https://github.com/sjspielman/dnds_1rate_2rate] (https://github.com/sjspielman/dnds_1rate_2rate) in the same directory as the current repository.

1. Translate simulated nucleotide alignments from Spielman et al. (2016) to amino acids using `./src/write_run_translate_aln.sh`.

2. Infer site-wise Rate4Site scores using `./src/write_run_r4s_mut_sel.sh`. This script will write `run_r4s_mut_sel.sh` which uses `r4s_pipeline.sh` to run Rate4Site on simulated alignments. 

3. Concatenate all rates into a table with `./src/concatenate_mut_sel_rates.r`. 

###Natural proteins

1. Align amino acid sequences using `./src/write_run_align_natural_prot.sh`.

2. Back translate amino acid alignments into codon alignments with `./src/run_translate_aln_aa_to_codon.sh`. This script requires original nucleotide sequences.

3. Infer trees from the amino acid sequences with RAxML. The script `./src/write_run_raxml.sh` will write `run_raxml.sh` which will run the inference. 

4. Infer site-wise dN/dS with `HyPhy` using the script in `dnds_inference/submit_run_inference_nat_prot.sh` directory in `src`. This script was copied from [https://github.com/sjspielman/dnds_1rate_2rate] (https://github.com/sjspielman/dnds_1rate_2rate) and modified for this analysis.

5. Infer site-wise Rate4Site scores using `./src/write_run_r4s_natural_prot.sh`. This script will write `run_r4s_natural_prot.sh` which will run Rate4Site on natural alignments. 

6. Concatenate all rates into a table with `./src/concatenate_natural_prot_rates.r`. 
