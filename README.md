Genome Analysis Pipeline
-------------------------

About
-----
The pipeline is for performing genome analysis 

After completing the make command kindly configure the dna_template.txt and rna_template.txt to point to the correct path for the tools

Submit a list of sample id's along with DNA/RNA in a tab separated list to create_config.pl with common path to create directory and the path for directory containing the fastq files

Once the config file for each sample is created a shell script within each working directory is created for initiating the analysis.

Steps:
------
(Refer to example provided)
1. make -f Makefile
2. perl target/create_config.pl target/sample_list . src/main/example/
3. cd ERR1025650
4. sh ERR1025650.sh

---
Author
Ankur
