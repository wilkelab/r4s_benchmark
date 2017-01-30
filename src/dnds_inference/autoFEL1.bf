HYPHYDIR = "placeholder/batchfiles/";
WDIR="placeholder/";
datafile="temp.fasta";
treefile="temp.tre";
nucfitfile="nuc.fit";
outfile="fel.txt";

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
