library(RColorBrewer)
library(VennDiagram)
library(dplyr)
library(sqldf)
library(Biostrings)
library(readxl)

setwd("B:/Doutorado/PatternLab/Outputs_Certos_Maio/")

#### lendo output PatternLab ####
splice_non_A <-read.table("PXD020656--spliceprot-COMPLETO-SEMNMD-leticia-11052022/peptides-PXD020656--spliceprot-COMPLETO-SEMNMD-leticia-11052022.txt", sep = "\t")
colnames(splice_non_A) <- c("Sequence", "Spectrum count", "Protein count", "The Protein Loci")
splice_non_A <- splice_non_A[!grepl("contaminant",splice_non_A$`The Protein Loci`),]
splice_non_A <- splice_non_A[!grepl("Reverse",splice_non_A$`The Protein Loci`),]

#### novo filtro por PeptideScan >2 ####

Scans <- read_excel("PXD020656--spliceprot-COMPLETO-SEMNMD-leticia-11052022/Scans-PXD020656--spliceprot-COMPLETO-SEMNMD-leticia-11052022.xlsx", 
                    col_types = c("numeric", "text", "text", 
                                  "text", "text", "text", "text", "text", 
                                  "numeric", "numeric", "text", "text", 
                                  "numeric", "numeric", "numeric", 
                                  "text", "text", "numeric", "text", 
                                  "numeric", "text", "text"))
#View(Scans)

teste2 <- Scans %>% filter(Scans$PrimaryScore >= 2)

common <- intersect(splice_non_A$Sequence,teste2$PeptideSequence)
common <- data.frame(common)
colnames(common) <- c('Sequence')

final <- sqldf('
select  *
from    splice_non_A two
        inner join common one
          on one.Sequence = two.Sequence')

a <- final[complete.cases(final), ]

splice_non_A <- data.frame(a$Sequence,a$`Spectrum count`,a$`Protein count`,a$`The Protein Loci`)
colnames(splice_non_A) <- c("Sequence", "Spectrum count", "Protein count", "The Protein Loci")
splice_non_A$Sequence <- gsub("\\[15.9949]","", splice_non_A$Sequence)
splice_non_A$Sequence <- gsub("\\[]","", splice_non_A$Sequence)

#splice_non_A$`The Protein Loci` <- gsub("ENST","", splice_non_A$`The Protein Loci`) #para versao com ENS no ProteinLoci

splice_non_REDUND <- splice_non_A %>% group_by(Sequence, 
                                               `The Protein Loci`,
                                               `Protein count`) %>% 
  summarize(`Spectrum count` = sum(`Spectrum count`))

splice_non_REDUND$Mesma.Isoforma <- rep(0, nrow(splice_non_REDUND))
# 2. descobre se tem mesma isoforma ou nÃ£o
for(x in 1:nrow(splice_non_REDUND)){
  print(paste0("Linha ",x," do data.frame"))
  b <- as.character(splice_non_REDUND$`The Protein Loci`[x])
  b <- splice_non_REDUND$`The Protein Loci`[x]
  b <- unlist(strsplit(b, "[|]"))
  b <- b[-1]
  s <- 0
  for(i in 1:length(b)){
    s<- s+ sum(!grepl(substr(b[i], 2,10 ), b[-i]))
  }
  splice_non_REDUND$Mesma.Isoforma[x] <- s
}
# 3. Filtra o data.frame pelas linhas em que deu mais do que 1 match nas isoformas
df_spliceprot_complete <- splice_non_REDUND[splice_non_REDUND$Mesma.Isoforma==0,]

spliceprot <- subset(df_spliceprot_complete, df_spliceprot_complete$`Protein count`==1)

#### swissprot ####

swiss <- read.table("PXD020656-swiss-leticia-13042022/peptides-PXD020656-swiss-leticia-13042022.txt", sep = "\t")
colnames(swiss) <- c("Sequence", "Spectrum count", "Protein count", "The Protein Loci")
swiss <- swiss[!grepl("contaminant",swiss$`The Protein Loci`),]
swiss <- swiss[!grepl("Reverse",swiss$`The Protein Loci`),]

#### novo filtro por PeptideScan >2 ####

Scans <- read_excel("PXD020656-swiss-leticia-13042022/Scans-PXD020656-swiss-leticia-13042022.xlsx", 
                    col_types = c("numeric", "text", "text", 
                                  "text", "text", "text", "text", "text", 
                                  "numeric", "numeric", "text", "text", 
                                  "numeric", "numeric", "numeric", 
                                  "text", "text", "numeric", "text", 
                                  "numeric", "text", "text"))

#View(Scans)

teste2 <- Scans %>% filter(Scans$PrimaryScore >= 2)

common <- intersect(swiss$Sequence,teste2$PeptideSequence)
common <- data.frame(common)
colnames(common) <- c('Sequence')

final <- sqldf('
select  *
from    swiss two
        join common one
          on one.Sequence = two.Sequence')

a <- final[complete.cases(final), ]

swiss <- data.frame(a$Sequence,a$`Spectrum count`,a$`Protein count`,a$`The Protein Loci`)

colnames(swiss) <- c("Sequence", "Spectrum count", "Protein count", "The Protein Loci")
swiss$Sequence <- gsub("\\[15.9949]","", swiss$Sequence)
swiss$Sequence <- gsub("\\[]","", swiss$Sequence)

swiss_non_REDUND <- swiss %>% group_by(Sequence, 
                                       `The Protein Loci`,
                                       `Protein count`) %>% 
  summarize(`Spectrum count` = sum(`Spectrum count`))


swissprot <- subset(swiss_non_REDUND, swiss_non_REDUND$`Protein count`==1)

#### encontrando as exclusivas do SpliceProt ####

common <- intersect(spliceprot$Sequence, swissprot$Sequence)

common <- data.frame(common)

require(plyr)

df <- join_all(list(spliceprot,swissprot), by = 'Sequence', type = 'inner')

redudantes_spliceANDswiss <- group_by(df, df$Sequence) %>% slice(1)

write.csv2(redudantes_spliceANDswiss, file = "PXD020656-redundantes-SPLICE-COMPLETO-swiss.csv")

#### comparando os 548 com os 2976

df = spliceprot %>% anti_join(swissprot,by="Sequence") #as 543

exclusivas_swissprot = swissprot %>% anti_join(spliceprot,by="Sequence") #exclusivas swiss

write.csv2(exclusivas_swissprot, file = "PXD020656-exclusivas_swissprot-SPLICEPROT-COMPLETO.csv")

exclusivas_spliceprot <- df

write.csv2(exclusivas_spliceprot, file = "PXD020656-exclusivas_spliceprot-SPLICE-COMPLETO-swiss.csv")

#### filtrando as exclusivas do spliceprot completo contra o banco completo do swissprot ####

swiss_banco_completo = Biostrings::readAAStringSet("B:/Doutorado/fastas/uniprot-created_[19860101+TO+20211117]+organism__mus+musculus_-filte--.fasta")

Header = names(swiss_banco_completo)
Sequence = paste(swiss_banco_completo)
swiss_banco_completo <- data.frame(Header, Sequence)

b <- sqldf('
select  *
from    exclusivas_spliceprot two
        left join swiss_banco_completo one
          on one.Sequence like "%" || two.Sequence || "%"
')

c <- b[complete.cases(b), ]

teste <- intersect(exclusivas_spliceprot$Sequence, c$Sequence)

teste <- data.frame(teste)

#agora falta ver oq nao tem em comum entre c e teste
proteotipicos_spliceprot <- exclusivas_spliceprot[!(exclusivas_spliceprot$Sequence %in% teste$teste),]
write.csv(proteotipicos_spliceprot,"PXD020656--spliceprot-COMPLETO-SEMNMD-leticia-11052022/proteotipicos-PXD020656-SPLICE-COMPLETO.csv")

#sair do R e rodar este script no último output 
#perl ~/disk1/leticia.costa/PatternLab/vin_filter_substrings_proteotipicos.pl proteotipicos-PXD020656-SPLICE-COMPLETO.csv > proteotipicos-PXD020656-SPLICE-COMPLETO-final.csv

#abrindo o novo arquivo de proteotípicos exclusivos do SpliceProt Completo

proteotipicos_spliceprot_final <- read_delim("PXD020656--spliceprot-COMPLETO-SEMNMD-leticia-11052022/proteotipicos-PXD020656-SPLICE-COMPLETO-final.csv", 
                                             delim = "\t", escape_double = FALSE, 
                                             col_names = FALSE, trim_ws = TRUE)
#mudando o nome das colunas do arquivo
colnames(proteotipicos_spliceprot_final) <- c('TheProteinLoci','Sequence','SpectrumCount')

#### openprot ####
open <- read.table("PXD020656-openprot-leticia-17042022/peptides-PXD020656-openprot-leticia-17042022.txt", sep = "\t")
colnames(open) <- c("Sequence", "Spectrum count", "Protein count", "The Protein Loci")
open <- open[!grepl("contaminant",open$`The Protein Loci`),]
open <- open[!grepl("Reverse",open$`The Protein Loci`),]

Scans <- read_excel("PXD020656-openprot-leticia-17042022/Scans-PXD020656-openprot-leticia-17042022.xlsx", 
                    col_types = c("numeric", "text", "text", 
                                  "text", "text", "text", "text", "text", 
                                  "numeric", "numeric", "text", "text", 
                                  "numeric", "numeric", "numeric", 
                                  "text", "text", "numeric", "text", 
                                  "numeric", "text", "text"))
#View(Scans)

teste2 <- Scans %>% filter(Scans$PrimaryScore >= 2)

common <- intersect(open$Sequence,teste2$PeptideSequence)
common <- data.frame(common)
colnames(common) <- c('Sequence')

final <- sqldf('
select  *
from    open two
        inner join common one
          on one.Sequence = two.Sequence')

a <- final[complete.cases(final), ]

open <- data.frame(a$Sequence,a$`Spectrum count`,a$`Protein count`,a$`The Protein Loci`)
colnames(open) <- c("Sequence", "Spectrum count", "Protein count", "The Protein Loci")
open$Sequence <- gsub("\\[15.9949]","", open$Sequence)
open$Sequence <- gsub("\\[]","", open$Sequence)

open_non_REDUND <- open %>% group_by(Sequence, 
                                     `The Protein Loci`,
                                     `Protein count`) %>% 
  summarize(`Spectrum count` = sum(`Spectrum count`))


openprot <- subset(open_non_REDUND, open_non_REDUND$`Protein count`==1)

#### encontrando as exclusivas do OpenProt ####

common <- intersect(openprot$Sequence, swissprot$Sequence)

common <- data.frame(common)

require(plyr)

df <- join_all(list(openprot,swissprot), by = 'Sequence', type = 'inner')

#redudantes_openANDswiss <- group_by(df, df$Sequence) %>% slice(1)

#write.csv2(redudantes_spliceANDswiss, file = "redundantes-SPLICE-COMPLETO-swiss.v2.csv")

#### comparando os 548 com os 2976

df = openprot %>% anti_join(swissprot,by="Sequence") #as 543

exclusivas_swissprot_open = swissprot %>% anti_join(openprot,by="Sequence") #exclusivas swiss

write.csv2(exclusivas_swissprot_open, file = "PXD020656-exclusivas_swissprot-OPENPROT.csv")

exclusivas_openprot <- df

write.csv2(exclusivas_openprot, file = "PXD020656-exclusivas_openprot-swiss.csv")

#### filtrando as exclusivas do openprot contra o banco completo do swissprot ####

b <- sqldf('
select  *
from    exclusivas_openprot two
        left join swiss_banco_completo one
          on one.Sequence like "%" || two.Sequence || "%"
')

c <- b[complete.cases(b), ]

teste <- intersect(exclusivas_openprot$Sequence, c$Sequence)

teste <- data.frame(teste)

#proteotípicos openprot
proteotipicos_openprot <- exclusivas_openprot[!(exclusivas_openprot$Sequence %in% teste$teste),]
write.csv(proteotipicos_openprot,"PXD020656-openprot-leticia-17042022/proteotipicos-PXD020656-openprot.csv")

#sair do R e rodar este script no último output 
#perl vin_filter_substrings_proteotipicos.pl PXD020656-proteotipicos-PXD020656-openprot.csv > PXD020656-proteotipicos-PXD020656-openprot-final.csv

#abrindo o novo arquivo de proteotípicos exclusivos do SpliceProt Completo
library(readr)
proteotipicos_openprot_final <- read_delim("PXD020656-openprot-leticia-17042022/proteotipicos-PXD020656-openprot-final.csv", 
                                             delim = "\t", escape_double = FALSE, 
                                             col_names = FALSE, trim_ws = TRUE)
#mudando o nome das colunas do arquivo
colnames(proteotipicos_openprot_final) <- c('TheProteinLoci','Sequence','SpectrumCount')

#### graficos ####

#### venn openprot vs swissprot ####
mycol <- c("yellow", "red")

venn.diagram(
  x = list(openprot$Sequence, swissprot$Sequence),
  category.names = c("OpenProt" , "SwissProt"),
  #category.names = c("" , ""),
  filename = 'B:/Doutorado/PatternLab/Outputs_Certos_Maio/Figuras/PXD020656/mouse-OPENPROT-swiss.png',
  output=TRUE,
  main = "MOUSE OPENPROT vs SWISSPROT",
  main.cex = 0.1,
  
  # Output features
  imagetype="png" ,
  height = 300 , 
  width = 400 , 
  resolution = 600, #dpi
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = mycol,
  
  # Numbers
  cex = .3,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.3,
  #cat.pos = 1,
  #cat.fontface = "bold",
  #cat.default.pos = "text",
  #cat.fontfamily = "sans",
)

#### venn spliceprot vs swissprot ####
mycol <- c("yellow", "red")

venn.diagram(
  x = list(spliceprot$Sequence, swissprot$Sequence),
  category.names = c("Splice" , "Swiss"),
  #category.names = c("" , ""),
  filename = 'B:/Doutorado/PatternLab/Outputs_Certos_Maio/Figuras/PXD020656/mouse-SPLICE-COMPLETO-swiss.png',
  output=TRUE,
  main = "MOUSE-SPLICE-COMPLETO vs SWISSPROT",
  main.cex = 0.2,
  
  # Output features
  imagetype="png" ,
  height = 300 , 
  width = 400 , 
  resolution = 600, #dpi
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = mycol,
  
  # Numbers
  cex = .3,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.3,
  #cat.pos = 1,
  #cat.fontface = "bold",
  #cat.default.pos = "text",
  #cat.fontfamily = "sans",
)

#### venn spliceprot vs openprot - proteotipicos exclusivos (resultado final) ####

#### venn spliceprot vs swissprot ####
mycol <- c("yellow", "red")

venn.diagram(
  x = list(exclusivas_spliceprot$Sequence, exclusivas_openprot$Sequence),
  category.names = c("Splice" , "Open"),
  #category.names = c("" , ""),
  filename = 'B:/Doutorado/PatternLab/Outputs_Certos_Maio/Figuras/PXD020656/mouse-SPLICE-COMPLETO-VS-OPENPROT.png',
  output=TRUE,
  main = "HUMAN-SPLICE-COMPLETO vs OPENPROT",
  main.cex = 0.2,
  
  # Output features
  imagetype="png" ,
  height = 300 , 
  width = 400 , 
  resolution = 600, #dpi
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = mycol,
  
  # Numbers
  cex = .3,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.3,
  #cat.pos = 1,
  #cat.fontface = "bold",
  #cat.default.pos = "text",
  #cat.fontfamily = "sans",
)

#### venn spliceprot vs swissprot filtro final ####
mycol <- c("yellow", "red")

venn.diagram(
  x = list(proteotipicos_spliceprot_final$Sequence, proteotipicos_openprot_final$Sequence),
  category.names = c("Splice" , "Open"),
  #category.names = c("" , ""),
  filename = 'B:/Doutorado/PatternLab/Outputs_Certos_Maio/Figuras/PXD020656/mouse-SPLICE-COMPLETO-VS-OPENPROT-filtrofinal.png',
  output=TRUE,
  main = "MOUSE-SPLICE-COMPLETO vs OPENPROT FILTRO FINAL",
  main.cex = 0.1,
  
  # Output features
  imagetype="png" ,
  height = 300 , 
  width = 400 , 
  resolution = 600, #dpi
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = mycol,
  
  # Numbers
  cex = .3,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.1,
  #cat.pos = 1,
  #cat.fontface = "bold",
  #cat.default.pos = "text",
  #cat.fontfamily = "sans",
)

#### venn proteotipicos spliceprot vs swissprot ####
mycol <- c("yellow", "red")

venn.diagram(
  x = list(proteotipicos_spliceprot_final$Sequence, swissprot$Sequence),
  category.names = c("Splice" , "Swiss"),
  #category.names = c("" , ""),
  filename = 'B:/Doutorado/PatternLab/Outputs_Certos_Maio/Figuras/PXD020656/mouse-PROTEOTIPICOS-SPLICE-COMPLETO-swiss.png',
  output=TRUE,
  main = "MOUSE-PROTEOTIPICOS-SPLICE-COMPLETO vs SWISSPROT",
  main.cex = 0.2,
  
  # Output features
  imagetype="png" ,
  height = 300 , 
  width = 400 , 
  resolution = 600, #dpi
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = mycol,
  
  # Numbers
  cex = .3,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.3,
  #cat.pos = 1,
  #cat.fontface = "bold",
  #cat.default.pos = "text",
  #cat.fontfamily = "sans",
)