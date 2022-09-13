library(RColorBrewer)
library(VennDiagram)
library(dplyr)
library(sqldf)
library(Biostrings)
library(readxl)
library(readr)

#### human ####

human <- read.table("B:/Doutorado/PatternLab/Outputs_Certos_Maio/SpliceProt_Digerido/PXD008720-spliceprot-DIGERIDO-SEMNMD-17052022/old/peptides-PXD008720-spliceprot-DIGERIDO-SEMNMD-17052022-old.txt", sep = "\t")
colnames(human) <- c("Sequence", "Spectrum count", "Protein count", "The Protein Loci")
human <- human[!grepl("contaminant",human$`The Protein Loci`),]
human <- human[!grepl("Reverse",human$`The Protein Loci`),]

#### novo filtro por PeptideScan >2 ####

Scans <- read_excel("Scans-PXD008720-spliceprot-DIGERIDO-SEMNMD-17052022-old.xlsx", 
                    col_types = c("numeric", "text", "text", 
                                  "text", "text", "text", "text", "text", 
                                  "numeric", "numeric", "text", "numeric", 
                                  "numeric", "numeric", "text", "text", 
                                  "numeric", "text", "numeric", "text", 
                                  "text"))

#View(Scans)

teste2 <- Scans %>% filter(Scans$PrimaryScore >= 2.5 & Scans$DeltaCN > 0.05)

common <- intersect(human$Sequence,teste2$PeptideSequence)
common <- data.frame(common)
colnames(common) <- c('Sequence')

final <- sqldf('
select  *
from    human two
        inner join common one
          on one.Sequence = two.Sequence')

a <- final[complete.cases(final), ]

human <- data.frame(a$Sequence,a$`Spectrum count`,a$`Protein count`,a$`The Protein Loci`)
colnames(human) <- c("Sequence", "Spectrum count", "Protein count", "The Protein Loci")
human$Sequence <- gsub("\\[15.9949]","", human$Sequence)
human$Sequence <- gsub("\\[]","", human$Sequence)

#human$`The Protein Loci` <- gsub("ENST","", human$`The Protein Loci`) #para versao com ENS no ProteinLoci

human_non_REDUND <- human %>% group_by(Sequence, 
                                       `The Protein Loci`,
                                       `Protein count`) %>% 
  summarize(`Spectrum count` = sum(`Spectrum count`))

human_non_REDUND$Mesma.Isoforma <- rep(0, nrow(human_non_REDUND))
# 2. descobre se tem mesma isoforma ou nÃ£o
for(x in 1:nrow(human_non_REDUND)){
  print(paste0("Linha ",x," do data.frame"))
  b <- as.character(human_non_REDUND$`The Protein Loci`[x])
  b <- human_non_REDUND$`The Protein Loci`[x]
  b <- unlist(strsplit(b, "[|]"))
  b <- b[-1]
  s <- 0
  for(i in 1:length(b)){
    s<- s+ sum(!grepl(substr(b[i], 2,10 ), b[-i]))
  }
  human_non_REDUND$Mesma.Isoforma[x] <- s
}
# 3. Filtra o data.frame pelas linhas em que deu mais do que 1 match nas isoformas
human_df_spliceprot_complete <- human_non_REDUND[human_non_REDUND$Mesma.Isoforma==0,]

human_spliceprot <- subset(human_df_spliceprot_complete, human_df_spliceprot_complete$`Protein count`==1)

#### mouse ####
mouse <- read.table("../Clean_PXD020656-splice-leticia-18112021_peptides.txt", sep = "\t")
colnames(mouse) <- c("Sequence", "Spectrum count", "Protein count", "The Protein Loci")
mouse$Sequence <- gsub("\\[15.9949]","", mouse$Sequence)
mouse$Sequence <- gsub("\\[]","", mouse$Sequence)

mouse_splice_non_REDUND <- mouse %>% group_by(Sequence, 
                                               `The Protein Loci`,
                                               `Protein count`) %>% 
  summarize(`Spectrum count` = sum(`Spectrum count`))

mouse_splice_non_REDUND$Mesma.Isoforma <- rep(0, nrow(mouse_splice_non_REDUND) )
# 2. descobre se tem mesma isoforma ou não
for(x in 1:nrow(mouse_splice_non_REDUND)){
  print(paste0("Linha ",x," do data.frame"))
  b <- as.character(mouse_splice_non_REDUND$`The Protein Loci`[x])
  b <- mouse_splice_non_REDUND$`The Protein Loci`[x]
  b <- unlist(strsplit(b, "[|]"))
  b <- b[-1]
  s <- 0
  for(i in 1:length(b)){
    s<- s+ sum(!grepl(substr(b[i], 2,10 ), b[-i]))
  }
  mouse_splice_non_REDUND$Mesma.Isoforma[x] <- s
}
# 3. Filtra o data.frame pelas linhas em que deu mais do que 1 match nas isoformas
mouse_df_spliceprot_complete <- mouse_splice_non_REDUND[mouse_splice_non_REDUND$Mesma.Isoforma==0,]

mouse_spliceprot <- subset(mouse_df_spliceprot_complete, mouse_df_spliceprot_complete$`Protein count`==1)

#### rat ####
rat <- read.table("../Clean_PXD016793-splice-leticia-181121_peptides.txt", sep = "\t")
colnames(rat) <- c("Sequence", "Spectrum count", "Protein count", "The Protein Loci")
rat$Sequence <- gsub("\\[15.9949]","", rat$Sequence)
rat$Sequence <- gsub("\\[]","", rat$Sequence)

rat_splice_non_REDUND <- rat %>% group_by(Sequence, 
                                               `The Protein Loci`,
                                               `Protein count`) %>% 
  summarize(`Spectrum count` = sum(`Spectrum count`))

rat_splice_non_REDUND $Mesma.Isoforma <- rep(0, nrow(rat_splice_non_REDUND ) )
# 2. descobre se tem mesma isoforma ou não
for(x in 1:nrow(rat_splice_non_REDUND )){
  print(paste0("Linha ",x," do data.frame"))
  b <- as.character(rat_splice_non_REDUND $`The Protein Loci`[x])
  b <- rat_splice_non_REDUND$`The Protein Loci`[x]
  b <- unlist(strsplit(b, "[|]"))
  b <- b[-1]
  s <- 0
  for(i in 1:length(b)){
    s<- s+ sum(!grepl(substr(b[i], 2,10 ), b[-i]))
  }
  rat_splice_non_REDUND $Mesma.Isoforma[x] <- s
}
# 3. Filtra o data.frame pelas linhas em que deu mais do que 1 match nas isoformas
rat_df_spliceprot_complete <- rat_splice_non_REDUND[rat_splice_non_REDUND$Mesma.Isoforma==0,]

rat_spliceprot <- subset(rat_df_spliceprot_complete, rat_df_spliceprot_complete$`Protein count`==1)

#### venn ####

mycol <- c("skyblue", "orange", "mediumorchid")

venn.diagram(
  x = list(human_spliceprot$Sequence, 
           mouse_spliceprot$Sequence, 
           rat_spliceprot$Sequence),
  category.names = c("Human" , "Mouse", "Rat"),
  #category.names = c("" , "", ""),
  filename = '../Figuras/mouse-venn-human-mouse-rat.png',
  output=TRUE,
  #main = "HUMAN",
  #main = "OpenProt ALL_DATABASE",
  #main.cex = 0.5,
  
  # Output features
  imagetype="png" ,
  height = 400 , 
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

#### comparando resultados SpliceProt Digerido ####
#### human ####

human <- read.table("B:/Doutorado/PatternLab/Outputs_Certos_Maio/SpliceProt_Digerido/PXD008720-spliceprot-DIGERIDO-SEMNMD-17052022/old/peptides-PXD008720-spliceprot-DIGERIDO-SEMNMD-17052022-old.txt", sep = "\t")
colnames(human) <- c("Sequence", "Spectrum count", "Protein count", "The Protein Loci")
human <- human[!grepl("contaminant",human$`The Protein Loci`),]
human <- human[!grepl("Reverse",human$`The Protein Loci`),]

#### novo filtro por PeptideScan >2 ####

Scans <- read_excel("B:/Doutorado/PatternLab/Outputs_Certos_Maio/SpliceProt_Digerido/PXD008720-spliceprot-DIGERIDO-SEMNMD-17052022/old/Scans-PXD008720-spliceprot-DIGERIDO-SEMNMD-17052022-old.xlsx", 
                    col_types = c("numeric", "text", "text", 
                                  "text", "text", "text", "text", "text", 
                                  "numeric", "numeric", "text", "numeric", 
                                  "numeric", "numeric", "text", "text", 
                                  "numeric", "text", "numeric", "text", 
                                  "text"))

#View(Scans)

teste2 <- Scans %>% filter(Scans$PrimaryScore >= 2.5 & Scans$DeltaCN > 0.05)

common <- intersect(human$Sequence,teste2$PeptideSequence)
common <- data.frame(common)
colnames(common) <- c('Sequence')

final <- sqldf('
select  *
from    human two
        inner join common one
          on one.Sequence = two.Sequence')

a <- final[complete.cases(final), ]

human <- data.frame(a$Sequence,a$`Spectrum count`,a$`Protein count`,a$`The Protein Loci`)
colnames(human) <- c("Sequence", "Spectrum count", "Protein count", "The Protein Loci")
human$Sequence <- gsub("\\[15.9949]","", human$Sequence)
human$Sequence <- gsub("\\[]","", human$Sequence)

human_non_REDUND <- human %>% group_by(Sequence, 
                                               `The Protein Loci`,
                                               `Protein count`) %>% 
  summarize(`Spectrum count` = sum(`Spectrum count`))

human_non_REDUND$Mesma.Isoforma <- rep(0, nrow(human_non_REDUND))
# 2. descobre se tem mesma isoforma ou nÃÂ£o
for(x in 1:nrow(human_non_REDUND)){
  print(paste0("Linha ",x," do data.frame"))
  b <- as.character(human_non_REDUND$`The Protein Loci`[x])
  b <- human_non_REDUND$`The Protein Loci`[x]
  b <- unlist(strsplit(b, "[|]"))
  b <- b[-1]
  s <- 0
  for(i in 1:length(b)){
    s<- s+ sum(!grepl(substr(b[i], 2,10 ), b[-i]))
  }
  human_non_REDUND$Mesma.Isoforma[x] <- s
}
# 3. Filtra o data.frame pelas linhas em que deu mais do que 1 match nas isoformas
human_df_spliceprot_complete <- human_non_REDUND[human_non_REDUND$Mesma.Isoforma==0,]

human_spliceprot <- subset(human_df_spliceprot_complete, human_df_spliceprot_complete$`Protein count`==1)

#### mouse ####
mouse <-read.table("B:/Doutorado/PatternLab/Outputs_Certos_Maio/SpliceProt_Digerido/PXD020656-spliceprot-DIGERIDO-SEMNMD-18052022/peptides-PXD020656-spliceprot-DIGERIDO-SEMNMD-18052022.txt", sep = "\t")
colnames(mouse) <- c("Sequence", "Spectrum count", "Protein count", "The Protein Loci")
mouse <- mouse[!grepl("contaminant",mouse$`The Protein Loci`),]
mouse <- mouse[!grepl("Reverse",mouse$`The Protein Loci`),]

#### novo filtro por PeptideScan >2 ####

Scans <- read_excel("B:/Doutorado/PatternLab/Outputs_Certos_Maio/SpliceProt_Digerido/PXD020656-spliceprot-DIGERIDO-SEMNMD-18052022/Scans-PXD020656-spliceprot-DIGERIDO-SEMNMD-18052022.xlsx", 
                    col_types = c("numeric", "text", "text", 
                                  "text", "text", "text", "text", "text", 
                                  "numeric", "numeric", "text", "numeric", 
                                  "numeric", "numeric", "text", "text", 
                                  "numeric", "text", "numeric", "text", 
                                  "text"))

#View(Scans)

teste2 <- Scans %>% filter(Scans$PrimaryScore >= 2.5 & Scans$DeltaCN > 0.05)

common <- intersect(mouse$Sequence,teste2$PeptideSequence)
common <- data.frame(common)
colnames(common) <- c('Sequence')

final <- sqldf('
select  *
from    mouse two
        inner join common one
          on one.Sequence = two.Sequence')

a <- final[complete.cases(final), ]

mouse <- data.frame(a$Sequence,a$`Spectrum count`,a$`Protein count`,a$`The Protein Loci`)
colnames(mouse) <- c("Sequence", "Spectrum count", "Protein count", "The Protein Loci")
mouse$Sequence <- gsub("\\[15.9949]","", mouse$Sequence)
mouse$Sequence <- gsub("\\[]","", mouse$Sequence)

mouse_non_REDUND <- mouse %>% group_by(Sequence, 
                                               `The Protein Loci`,
                                               `Protein count`) %>% 
  summarize(`Spectrum count` = sum(`Spectrum count`))

mouse_non_REDUND$Mesma.Isoforma <- rep(0, nrow(mouse_non_REDUND))
# 2. descobre se tem mesma isoforma ou nÃ£o
for(x in 1:nrow(mouse_non_REDUND)){
  print(paste0("Linha ",x," do data.frame"))
  b <- as.character(mouse_non_REDUND$`The Protein Loci`[x])
  b <- mouse_non_REDUND$`The Protein Loci`[x]
  b <- unlist(strsplit(b, "[|]"))
  b <- b[-1]
  s <- 0
  for(i in 1:length(b)){
    s<- s+ sum(!grepl(substr(b[i], 2,10 ), b[-i]))
  }
  mouse_non_REDUND$Mesma.Isoforma[x] <- s
}
# 3. Filtra o data.frame pelas linhas em que deu mais do que 1 match nas isoformas
mouse_df_spliceprot_complete <- mouse_non_REDUND[mouse_non_REDUND$Mesma.Isoforma==0,]

mouse_spliceprot <- subset(mouse_df_spliceprot_complete, mouse_df_spliceprot_complete$`Protein count`==1)

#### rat ####

#### lendo output PatternLab ####
rat <-read.table("B:/Doutorado/PatternLab/Outputs_Certos_Maio/SpliceProt_Digerido/PXD016793-spliceprot-DIGERIDO-SEMNMD-17052022/peptides-PXD016793-spliceprot-DIGERIDO-SEMNMD-17052022.txt", sep = "\t")
colnames(rat) <- c("Sequence", "Spectrum count", "Protein count", "The Protein Loci")
rat <- rat[!grepl("contaminant",rat$`The Protein Loci`),]
rat <- rat[!grepl("Reverse",rat$`The Protein Loci`),]

#### novo filtro por PeptideScan >2 ####

Scans <- read_excel("B:/Doutorado/PatternLab/Outputs_Certos_Maio/SpliceProt_Digerido/PXD016793-spliceprot-DIGERIDO-SEMNMD-17052022/Scans-PXD016793-spliceprot-DIGERIDO-SEMNMD-17052022.xlsx", 
                    col_types = c("numeric", "text", "numeric", 
                                  "text", "text", "text", "text", "text", 
                                  "numeric", "numeric", "text", "text", 
                                  "numeric", "numeric", "numeric", 
                                  "text", "text", "numeric", "text", 
                                  "numeric", "text", "text"))

#View(Scans)

teste2 <- Scans %>% filter(Scans$PrimaryScore >= 2.5 & Scans$DeltaCN > 0.05)

common <- intersect(rat$Sequence,teste2$PeptideSequence)
common <- data.frame(common)
colnames(common) <- c('Sequence')

final <- sqldf('
select  *
from    rat two
        inner join common one
          on one.Sequence = two.Sequence')

a <- final[complete.cases(final), ]

rat <- data.frame(a$Sequence,a$`Spectrum count`,a$`Protein count`,a$`The Protein Loci`)
colnames(rat) <- c("Sequence", "Spectrum count", "Protein count", "The Protein Loci")
rat$Sequence <- gsub("\\[15.9949]","", rat$Sequence)
rat$Sequence <- gsub("\\[]","", rat$Sequence)

rat_non_REDUND <- rat %>% group_by(Sequence, 
                                               `The Protein Loci`,
                                               `Protein count`) %>% 
  summarize(`Spectrum count` = sum(`Spectrum count`))

rat_non_REDUND$Mesma.Isoforma <- rep(0, nrow(rat_non_REDUND))
# 2. descobre se tem mesma isoforma ou nÃ£o
for(x in 1:nrow(rat_non_REDUND)){
  print(paste0("Linha ",x," do data.frame"))
  b <- as.character(rat_non_REDUND$`The Protein Loci`[x])
  b <- rat_non_REDUND$`The Protein Loci`[x]
  b <- unlist(strsplit(b, "[|]"))
  b <- b[-1]
  s <- 0
  for(i in 1:length(b)){
    s<- s+ sum(!grepl(substr(b[i], 2,10 ), b[-i]))
  }
  rat_non_REDUND$Mesma.Isoforma[x] <- s
}
# 3. Filtra o data.frame pelas linhas em que deu mais do que 1 match nas isoformas
rat_df_spliceprot_complete <- rat_non_REDUND[rat_non_REDUND$Mesma.Isoforma==0,]

rat_spliceprot <- subset(rat_df_spliceprot_complete, rat_df_spliceprot_complete$`Protein count`==1)

#### venn ####

mycol <- c("skyblue", "orange", "mediumorchid")

venn.diagram(
  x = list(human_spliceprot$Sequence, 
           mouse_spliceprot$Sequence, 
           rat_spliceprot$Sequence),
  category.names = c("Human" , "Mouse", "Rat"),
  #category.names = c("" , "", ""),
  filename = 'B:/Doutorado/PatternLab/Outputs_Certos_Maio/SpliceProt_Digerido/Figuras/mouse-venn-human-mouse-rat.png',
  output=TRUE,
  #main = "HUMAN",
  #main = "OpenProt ALL_DATABASE",
  #main.cex = 0.5,
  
  # Output features
  imagetype="png" ,
  height = 400 , 
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