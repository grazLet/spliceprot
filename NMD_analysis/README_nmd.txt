You can Download NMD Classifier at https://sourceforge.net/projects/transcriptome-analysis/files/NMD_Classifier.tar.gz/download



The config.txt used by NMD_analysis.pl script was filled with the path of following files:
"our_dataset.gtf", "ensembl_anottation.gtf", "genomic_chromossome_fasta_seq/*.fa", "output_folder/"

Our in house script was used to apply NMD-tagets filter in our data set

removeNMD.v3.pl NMD_Classifier_output.txt SpliceProt_fasta.fa > New_fasta_withoutNMD.fa
