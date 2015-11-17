/*
    Run as HYPHYMP CPU=<#> autoFEL1.bf

    This script runs a 1-rate FEL to compute dN/dS. Uses the MG94xHKY85 model with F1x4 state frequencies. It also assumes provided branch lengths are the optimized branch lengths.
*/

HYPHYDIR = "batchfiles/"; // replace "placeholder" with FULL PATH to working directory 
WDIR=".";                   // also here
datafile="ENSGT00640000091523.Euteleostomi.001.nt.fas";              // input fasta file, WITHOUT tree
treefile="ENSGT00640000091523.Euteleostomi.001.nwk";                  // input tree file
nucfitfile="nuc.fit";               // whatever you want
outfile="test.txt";                   // output file name

inputRedirect = {};
inputRedirect["01"]="Universal";        // Genetic Code
inputRedirect["02"]="New Analysis";     // New analysis
inputRedirect["03"]=WDIR+datafile;      // Alignment file
inputRedirect["04"]="Default";          // Use HKY85 and MG94xHKY85.
inputRedirect["05"]=WDIR+treefile;      // Tree file
inputRedirect["06"]=WDIR+nucfitfile;    // Save nucleotide fit to..
inputRedirect["07"]="Estimate dN/dS only";  // Estimate dN/dS w/out branch corrections
inputRedirect["08"]="One rate FEL";     // 1-rate FEL (dS constant across sites)
inputRedirect["09"]="0.01";             // FPR for LRT
inputRedirect["10"]=WDIR+outfile;       // output file


ExecuteAFile (HYPHYDIR + "QuickSelectionDetection_sjs.bf", inputRedirect);
