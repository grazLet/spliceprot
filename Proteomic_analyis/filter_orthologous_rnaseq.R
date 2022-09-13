library("tximport")
library("readr")
library("dplyr")
library('biomaRt')

#humano
csv.dir <- "B:/Doutorado/ortologos_flavia/resultados_finais/rnaseq/humano/quant/"
samps <- read_delim("B:/Doutorado/ortologos_flavia/resultados_finais/rnaseq/humano/samples_human_liver.csv", 
    delim = ";", escape_double = FALSE, trim_ws = TRUE
)
#head(samps)

samps$condition <- factor(samps$condition)
#arquivos
files <- paste0(csv.dir, samps$sample_id)
files <- paste0(files, '/quant.sf')
names(files) <- samps$sample_id
quants <- read_tsv(files)
#head(quants)
quants$Name <- gsub("\\..*","",quants$Name)

quants_mean <- quants %>%
  group_by(Name) %>%
  dplyr::summarize(TPM_Mean = mean(TPM, na.rm=TRUE)
)

target <- c("ENST00000368235","ENST00000380191","ENST00000539909",
    "ENST00000305218","ENST00000229329","ENST00000458591",
    "ENST00000340124","ENST00000308696","ENST00000305208",
    "ENST00000353703","ENST00000372839","ENST00000207870",
    "ENST00000436784","ENST00000265517","ENST00000280551",
    "ENST00000520121","ENST00000651975","ENST00000274364",
    "ENST00000230792","ENST00000304858","ENST00000435765",
    "ENST00000292644","ENST00000457587","ENST00000425206",
    "ENST00000277165","ENST00000398805","ENST00000336079",
    "ENST00000360160"
)

alvos <- dplyr::filter(quants_mean, Name %in% target)

alvos_expressos <- alvos %>% 
  dplyr::filter(alvos$TPM_Mean > 0
)

mart = useMart(host="useast.ensembl.org", 
               biomart="ENSEMBL_MART_ENSEMBL", 
               dataset="hsapiens_gene_ensembl"
)

G_list <- getBM(filters= "ensembl_transcript_id", attributes= c("ensembl_transcript_id","hgnc_symbol"),values=alvos_expressos$Name,mart= mart)

unicos_humano_genes <- data.frame(unique(G_list$hgnc_symbol))

write.csv2(alvos_expressos,"B:/Doutorado/ortologos_flavia/resultados_finais/rnaseq/humano/transcritos_ortologos_expressos_humano.csv")

write.csv2(G_list,"B:/Doutorado/ortologos_flavia/resultados_finais/rnaseq/humano/lista-genes_transcritos_ortologos_expressos_humano.csv")

#camundongo
csv.dir <- "B:/Doutorado/ortologos_flavia/resultados_finais/rnaseq/camundongo/quant/"

samps <- read_delim("B:/Doutorado/ortologos_flavia/resultados_finais/rnaseq/camundongo/samples_mouse_liver.csv", 
                    delim = ";", escape_double = FALSE, trim_ws = TRUE
)

#head(samps)
samps$condition <- factor(samps$condition)
#arquivos
files <- paste0(csv.dir, samps$sample_id)
files <- paste0(files, '/quant.sf')
names(files) <- samps$sample_id
quants <- read_tsv(files)
#head(quants)
quants$Name <- gsub("\\..*","",quants$Name)

quants_mean <- quants %>%
  group_by(Name) %>%
  dplyr::summarize(TPM_Mean = mean(TPM, na.rm=TRUE)
)

target <- c("ENSMUST00000029708","ENSMUST00000223396","ENSMUST00000045376",
    "ENSMUST00000025918","ENSMUST00000032419","ENSMUST00000021335",
    "ENSMUST00000034136","ENSMUST00000027432","ENSMUST00000073049",
    "ENSMUST00000018470","ENSMUST00000131288","ENSMUST00000039610",
    "ENSMUST00000112543","ENSMUST00000159809","ENSMUST00000036382",
    "ENSMUST00000162562","ENSMUST00000029805","ENSMUST00000047923",
    "ENSMUST00000058578","ENSMUST00000166581","ENSMUST00000068603",
    "ENSMUST00000025065","ENSMUST00000020630","ENSMUST00000030769",
    "ENSMUST00000060805","ENSMUST00000030078","ENSMUST00000091190"
)

alvos <- dplyr::filter(quants_mean, Name %in% target)

alvos_expressos <- alvos %>% 
  dplyr::filter(alvos$TPM_Mean > 0
)

mart = useMart(host="useast.ensembl.org", 
               biomart="ENSEMBL_MART_ENSEMBL", 
               dataset="mmusculus_gene_ensembl"
)

G_list <- getBM(filters= "ensembl_transcript_id", attributes= c("ensembl_transcript_id","mgi_symbol"),values=alvos_expressos$Name,mart= mart)

unicos_camundongo_genes <- data.frame(unique(G_list$mgi_symbol))

write.csv2(alvos_expressos,"B:/Doutorado/ortologos_flavia/resultados_finais/rnaseq/camundongo/transcritos_ortologos_expressos_camundongo.csv")

write.csv2(G_list,"B:/Doutorado/ortologos_flavia/resultados_finais/rnaseq/camundongo/lista-genes_transcritos_ortologos_expressos_camundongo.csv")

#rato

csv.dir <- "B:/Doutorado/ortologos_flavia/resultados_finais/rnaseq/rato/quant/"
samps <- read_delim("B:/Doutorado/ortologos_flavia/resultados_finais/rnaseq/rato/samples_rat_liver.csv", 
    delim = ";", escape_double = FALSE, trim_ws = TRUE
)

#head(samps)
samps$condition <- factor(samps$condition)
#arquivos
files <- paste0(csv.dir, samps$sample_id)
files <- paste0(files, '/quant.sf')
names(files) <- samps$sample_id
quants <- read_tsv(files)
#head(quants)
quants$Name <- gsub("\\..*","",quants$Name)

quants_mean <- quants %>%
  group_by(Name) %>%
  dplyr::summarize(TPM_Mean = mean(TPM, na.rm=TRUE)
)

target <- c("ENSRNOT00000025986","ENSRNOT00000024952","ENSRNOT00000016709",
    "ENSRNOT00000028743","ENSRNOT00000018734","ENSRNOT00000040548",
    "ENSRNOT00000077275","ENSRNOT00000024306","ENSRNOT00000025045",
    "ENSRNOT00000016981","ENSRNOT00000019106","ENSRNOT00000074595",
    "ENSRNOT00000014631","ENSRNOT00000064809","ENSRNOT00000018796",
    "ENSRNOT00000064091","ENSRNOT00000035017","ENSRNOT00000066968",
    "ENSRNOT00000023628","ENSRNOT00000016450","ENSRNOT00000060568",
    "ENSRNOT00000059458","ENSRNOT00000092078"
)

alvos <- dplyr::filter(quants, Name %in% target)

alvos_expressos <- alvos %>% 
  dplyr::filter(alvos$TPM_Mean > 0
)

mart = useMart(host="useast.ensembl.org", 
  biomart="ENSEMBL_MART_ENSEMBL", 
  dataset="rnorvegicus_gene_ensembl"
)

G_list <- getBM(filters= "ensembl_transcript_id", attributes= c("ensembl_transcript_id","external_gene_name"),values=alvos_expressos$Name,mart= mart)

unicos_rato_genes <- data.frame(unique(G_list$external_gene_name))

write.csv2(alvos_expressos,"B:/Doutorado/ortologos_flavia/resultados_finais/rnaseq/rato/transcritos_ortologos_expressos_rato.csv")

write.csv2(G_list,"B:/Doutorado/ortologos_flavia/resultados_finais/rnaseq/rato/lista-genes_transcritos_ortologos_expressos_rato.csv")
