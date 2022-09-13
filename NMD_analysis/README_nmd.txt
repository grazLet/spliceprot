last update:13/09/2020



You can Download NMD Classifier at https://sourceforge.net/projects/transcriptome-analysis/files/NMD_Classifier.tar.gz/download

The config.txt used by NMD_analysis.pl script was filled with the path of following files:
"our_dataset.gtf", "ensembl_anottation.gtf", "genomic_chromossome_fasta_seq/*.fa", "output_folder/"

The "our_dataset.gtf" was made by using SpliceProt intermediary tables: variant_file, consenso_file, cluster_file and anotation_file'
And we construct a script to generate GTF file for each specie.

human_GTFfromVariants.v3.pl - for human
Rat_GTFfromVariants.v3.pl - for Rat
Mmu_GTFfromVariants.v3.pl - for mouse

Our in house script "removeNMD.v3.pl" was used to apply NMD-tagets filter in our data set

removeNMD.v3.pl NMD_Classifier_output.txt SpliceProt_fasta.fa > New_fasta_withoutNMD.fa
