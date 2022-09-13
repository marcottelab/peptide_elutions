# peptide_elutions

Analysis accompanying "Alternative proteoforms and proteoform-dependent assemblies in humans and plants"

Claire D. McWhite, Wisath Sae-Lee, Yawning Yuan, Anna L. Mallam, 
Nicolas A. Gort Freitas, Silvia Ramundo, Masayuki Onishi, and Edward M. Marcotte


The proteoform analysis comprises of 3 parts: processing, scoring, and vizualizing.  

### Input data
**1. peptide information from fractionation experiment**<Enter>
	
If the fractionation experiment was analyzed with MSFragger, run the following script to combine data from each fraction into one single file.<Enter> 

Example code: python3 /scripts/format_MSFragger_files.py --root_folder **folder where your results from MSFragger analysis are located** --fractionation_name **name of your fractionation experiment** --output_file **name of outputfile wide format** --fraction_order **name of outputfile fraction order file**

Expected result:<Enter> 
![alt text](https://user-images.githubusercontent.com/32718019/187560693-c5e8851d-a7cc-4705-bf01-0b6e575f1673.png)
	
**2. In-silico digest peptides**
	
Example code: python2.7 /scripts/trypsin.py --input_file /test/uniprot_human.fasta --output_file uniprot_human_digested.csv --miss 2 --positions True
	
Expected result:<Enter>
	
![alt text](https://user-images.githubusercontent.com/32718019/188028137-ccdc1511-13e7-40ff-883f-f5075daf1ed1.png)



### Processing
Process peptide files for Gaussian fitting 
Script: peptide_identification_single_frac.R <Enter> 

Input file: <Enter> 
![alt text](https://user-images.githubusercontent.com/32718019/187560693-c5e8851d-a7cc-4705-bf01-0b6e575f1673.png)

Example code: Rscript-4.0.3 /scripts/peptide_identification_single_frac.R --elut_wide_file /data/pivot_test.csv --fraction_order /data/fraction_order_test.csv --peps /data/uniprot_human_digested.csv --seqlen /data/seq_length_homo_sapiens.tsv --spec human --output_file /data/short_tidy_unique_MB_sup_SEC.csv  

Expected result:<Enter> 
![alt text](https://user-images.githubusercontent.com/32718019/187561244-b0e6ed26-6e5c-4f65-ab26-d1386462185d.png)


### Scoring
**1. Identify peaks from peptide elution profile of each protein in a fractionation experiment using Gaussian Mixture Model. Multiple peaks suggest the existence of proteoforms or intact proteins eluting with different binding partners.** <Enter> 
	
Script: Gaussian_fitting.R <Enter> 
	
Input file: <Enter> 
![alt text](https://user-images.githubusercontent.com/32718019/186778488-8172fdfc-f8d8-400b-89ac-76dec4752308.png) <Enter> 

Example code: Rscript-4.0.3 /scripts/Gaussian_fitting.R --input_file example/short_tidy_unique_anna_hekSEC2.csv --simple_AdapGauss /scripts/simple_AdaptGauss.R --output_file example/short_tidy_unique_anna_hekSEC2_peaks.csv <Enter> 

Expected result:<Enter> 
	
![alt text](https://user-images.githubusercontent.com/32718019/186778945-6d2824fb-8350-4787-825d-6908834f9f9a.png)

**2. Calculate terminal bias score in order to prioritize proteins to inspect manually.**

Input file 

1.peptide file (Same as input file for Gaussian fitting step)

Input file: <Enter> 
	
![alt text](https://user-images.githubusercontent.com/32718019/186778488-8172fdfc-f8d8-400b-89ac-76dec4752308.png) <Enter> 

2.peak file (From previous Gaussian fitting step)
	
![alt text](https://user-images.githubusercontent.com/32718019/186778945-6d2824fb-8350-4787-825d-6908834f9f9a.png)

Example code: Rscript-4.0.3 /scripts/terminal_bias.R --input_file /test/short_tidy_unique_anna_hekSEC2_78.csv --peaks /peaks_short_tidy_unique_anna_hekSEC2_78.csv --output_file /test/terminal_bias_short_tidy_unique_anna_hekSEC2_78.csv <Enter>

Expected result:<Enter> 
![alt text](https://user-images.githubusercontent.com/32718019/189236527-15524bef-9682-4d47-a436-69203dcd0a1f.png)

This example shown in the expected result above is the terminal bias score (abslog2fc) for each Gaussian peak for PUR2 (see Figure 4A in the manuscript). The higher score suggests the existence of the proteoform, but would require further manual inspection. Terminal bias score is used to narrow down the list of proteins to inspect manually.
	
For example, this histogram shows the distribution of terminal bias scores for the proteins in a size exclusion fractionation of HEK293T lysate. As shown in Figure 4A, both the full length and short proteoforms of PUR2 can be detected based on their peptide elution profile. The terminal bias score for PUR2 is indicated by the red arrow. For the initial inspection, we recommend examining proteins with terminal bias score > 3, which corresponds to ~150 proteins from this experiment, for example.

![My Image](example_files/tb_dist.png)

### Visual categorizing

Visualize proteins through a Shiny app. The example below demonstrates how to view peptide elution profile for PUR2 from HEK293T cell fractionation. Run elution_viewer_simplified.Rmd in Rstudio.
	
Input file 
1. The same input file as the first step in scoring.
	
![alt text](https://user-images.githubusercontent.com/32718019/186778488-8172fdfc-f8d8-400b-89ac-76dec4752308.png) <Enter> 
	
2. Domain information

Information on domains were taken from Interpro (https://www.ebi.ac.uk/interpro/) and MobiDB (https://mobidb.bio.unipd.it/). The format of the input table for peptide_elution viewer is shown below. For convenience, we provide domain information files fro human, Arabidopsis, and Chlamydomonas in data_files (domaindis_setup.txt).
	
![alt text](https://user-images.githubusercontent.com/32718019/188288836-550608be-d1a2-405c-982b-19b7508fe2cc.png)

Expected result:<Enter> 

![My Image](example_files/PUR2_elution_viewer.png)

In order to identify the breakpoint for the proteoform and inspect the peptides observed for a particular proteoform further, you can highlight the peptides of interest in the Shiny app as shown below:

![My Image](example_files/PUR2_elut_viewer_zoom.png)


The information on peptides can be exported from the Shiny app. By looking at the last observed peptide for a proteoform, you can identify the breakpoint for the proteoform. For this PUR2 proteoform, the break point is residue 433. 

![image](https://user-images.githubusercontent.com/32718019/189445819-22f64d2d-ee79-4d26-8383-92862f4d10a7.png)

