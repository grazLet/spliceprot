last update:15/09/2020

###after the selection of the best hits of psl alignment file we need to...

##remove some remanescent redundances 

perl ../Scripts/removeduplicates.pl alinhamento_hg38_release99_reps.psl

###Next scripts was used to generate files that we use to create SpliceProt 2.0 repository.

##  ENSG vs ENST dictionary

# remember to download the fasta file of cdna for each specie 

perl dicioENSTandENSG.pl ../alinhamentos/Homo_sapiens.GRCh38.cdnaANDncrna.fa ENSTversusENSG

##Criando a tabela clusters 

perl FormatSQLtable.BLAST9ver.v3.1.pl Arquivos/ENSTversusENSG.tsv alinhamentos/alinhamento_hg38_release99_reps.psl.out.semrep.psl alinhamentos/alinhamento_hg38_release99.blast9 clusters_release99.csv

## Dump & Restore of database

psql -U leticia -f Homo_sapiens_hg38.sql Homo_sapiens_release100

## populating the cluster table

\copy clusters (id_clusters, chr, seq_acc, start_chr, end_chr, start_est, end_est, strand, identity, cluster_id) FROM '/home/leticia/Documentos/Doutorado/SpliceProt/SpliceProt_Ensembl_99/clusters_release99.csv' with delimiter as E'\t' CSV HEADER;
s
## creating the file to TSL table

wget ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.chr.gtf.gz

gunzip Homo_sapiens.GRCh38.100.chr.gtf.gz

perl buildTSLtable.pl alinhamentos/Homo_sapiens.GRCh38.cdnaANDncrna.fa Arquivos/Homo_sapiens.GRCh38.99.chr.gtf > TSL_release99.csv

## populating the TSL table

\COPY TSL (id_tsl, nm_transcrito, tsl) FROM '/home/leticia/Documentos/Doutorado/SpliceProt/SpliceProt_Ensembl_100/Homo_sapiens/Arquivos/TSL_release100.csv' with delimiter as E'\t' CSV HEADER;

## creating the file to clust_seq table
perl buildTSLndSeq.pl alinhamentos/Homo_sapiens.GRCh38.cdnaANDncrna.fa Arquivos/Homo_sapiens.GRCh38.99.chr.gtf

## modifications of clust_seq table

\copy clust_seq_2 (cluster_id, sequence, seq_acc, tissue, tsl) FROM '/home/leticia/Documentos/Doutorado/SpliceProt/SpliceProt_Ensembl_99/Arquivos/clust_seq_release99.csv' with delimiter as E'\t' CSV HEADER;

update smart_est.clust_seq_2 set cluster_id = substring(cluster_id from 5 for 15); #humano 
update smart_est.clust_seq_2 set cluster_id = substring(cluster_id from 8 for 15); #camundongo
ALTER TABLE smart_est.clust_seq_2 ALTER COLUMN cluster_id TYPE integer USING cluster_id::integer;

## criando o arquivo para a tabela select_consenso 

#humano perl Select_consenso_human.pl ENST_seq_TSL.csv > select_consenso.csv

#camundongo perl Select_consenso_mouse.pl ../alinhamentos/Arquivos/ENST_seq_TSL.csv > select_consenso.csv
