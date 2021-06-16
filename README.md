# Blanck Group Code Repository

Code Package A
- Performs pairwise alignment functionality between sequencing read of interest and known V and J regions.
- Performs analysis of V-J junction to determine whether the read of interest contains productive junction (in-frame and without stop codons)
- Instructions for use:
    1. Instantiate object corresponding to immune receptor. Arguments passed must include read of interest (sequence), and path to known V and J regions (dbpath). Known region files must be fasta files and be names by their immune receptor and V or J region. E.g. "TRAV.fasta". vthreshold and jthreshold optional arguments determine cutoff values for pairwise alignment scoring. 
    2. Calling .run() method on the immune receptor object will perform full pairwise alignment and junction analysis. Returns pandas dataframe of results. 

Code Package B
- VDJRecord class performs the follwing functions:
    1. Uses Code Package A objects to perform VDJ analysis on path containing sequence read files. Further filters results based on match length and percent match. 
    2. Analyses physico-chemical properties of amino acid sequences of CDR3 domains. 
- VDJAnalyze class performs kaplan-meier surival analysis

Code Package C
- Script for performing Kaplan-Meier and Cox regression survival analysis on physicochemical properties of CDR3 domains. 

Code Package D
- Tools for complementarity scoring of CDR3 domains against potential antigen sequences. 
- Script takes directory of results generated by Code Package B VDJRecord object and csv file of gene mutations. 
- Can take directory of multiple results folders. 

Code Package E
- Script for Pearson correlation of complementarity scores generated by Code Package D against gene expression values. 
- Script for Cox regression survival analysis of gene expression. 

Code Package F
- Script for sample-dropout simulation of p-values 
