# Teleosts
These different steps allow the analysis of the human orthology genes present in different species, and their statistics. 
The programme returns a csv file containing the statistical tests (chi2 and hypergeometry) for each species, but also 2 files containing the list of all human genes and the presence or absence of orthology in the different species. 1 file for duplicate orthologies and 1 file for singleton orthologies. 


## Step
1. Download in csv format for different species the orthologs of human genes encoding a protein on Biomart Ensembl (https://www.ensembl.org/biomart/martview/c92f5711966e86d6f9fd7bc989766704) 
2. Put all files uploaded in a file marked "DATA".
3. Download the files _launch_analysis.sh_, _step1.py_ and _step2.R_ and put them in the parent directory of the DATA directory.
4. Start the program _launch_analysis.sh_: `./launch_analysis.sh`.
