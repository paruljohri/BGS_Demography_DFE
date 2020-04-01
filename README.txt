Towards an Evolutionarily Appropriate Null Model: Jointly Inferring Demography and Purifying Selection

#Contact pjohri1@asu.edu for questions.

This folder has scripts and statistics related to Johri et al 2020 Genetics (https://www.genetics.org/content/early/2020/03/09/genetics.119.303002.abstract). Following are the list of contents:

a) Folder: Statistics.zip
1. Files that have statistics calculated for forward simulations with background selection in demographic non-equilibrium:
	demo_disc_5_SingExon_sumstats_50.txt
	demo_disc_5_SingExon_sumstats_pi.txt
2. File that has statistics calculated for forward simulations with background selection in demographic equilibrium (500 bp):
	eqm_disc_5_500bp_sumstats_50.txt
	eqm_disc_5_500bp_sumstats_pi.txt
3. File that has statistics calculated for forward simulations with background selection under demographic equilibrium (1 kb):
	eqm_disc_5_1kb_sumstats_50.txt
	eqm_disc_5_1kb_sumstats_pi.txt
4. File that has statistics calculated for forward simulations with background selection under demographic equilibrium (5 kb):
	eqm_disc_5_5kb_sumstats_50.txt
	eqm_disc_5_5kb_sumstats_pi.txt
5. File that has statistics calculated for forward simulations with background selection under demographic equilibrium (10 kb):
	eqm_disc_5_10kb_sumstats_50.txt
	eqm_disc_5_10kb_sumstats_pi.txt

b)Folder: Dpgp3_data.zip
1. File with names of the 94 exons used in this study, along with their pairwise divergence and rates of recombination:
	sing_exon_500_2000_4000_prop_div.txt
2. Folder with aligned fasta files of all 94 exons with their ancestral sequence. All alignements begin with a 4kb intergenic region and end with the single-exon functional element:
	SingleExon94_Fasta_phastCons0.8
3. Folder with .ms files that contain derived allele (=1) information in .ms format.
	SingleExon94_DerivedAllele_phastCons0.8
3. Folder with .fixed files that contain alignement locations where fixed differences were found between the ancestor and 76 individuals from Zambia.
	SingleExon94_FixedSubstitutions_phastCons0.8

c)Folder: Scripts
1. Slim script for simulations with background selection and demographic non-equilibrium:
	demo_disc_5_SingExon_osg.slim
2. Script to calculate sliding window statistics from simulated data:
	statistics_slidingwindow_pylibseq_SingExon_osg.py
3. Script to fit the recovery of nucleotide diversity w.r.t distance from selected sites:
	get_pirecovery_statistics.R
4. Script to calculate statistics in neutral, linked and functional regions around simulated functional regions:
	statistics_bigwindow_pylibseq_SingExon_osg.py
5. Script to calculate statistics in neutral, linked and exons in dpgp3 data (using pylibseq):
	statistics_bigwindow_pylibseq_dpgp3.py
	get_final_statistics_dpgp3.R
6. Script to calculate LD statistics in neutral, linked and exons in dpgp3 data which incorporates missing data:
	statistics_bigwindow_mystats_dpgp3.py

D) Mathematica notebook to show the calculation of analytical solution to calculate the reduction in nucleotide diversity due to BGS caused by a functional element:
calculateB_Equation3_solution

E)Python script to show the calculation of analytical solution to calculate the reduction in nucleotide diversity due to BGS caused by a functional element:
calculate_B_analytically_Eq3
