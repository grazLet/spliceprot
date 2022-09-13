library(readr)
library(dplyr)

####lendo arquivo hsap vs mmus ####
hsap_vs_mmus <- read_table2("B:/Doutorado/ortologos_flavia/hsap-vs-mmus/spliceprot_FinalTab_App1.txt", 
                            col_types = cols(similarity = col_skip(), 
                                             similarityRatio = col_number(), identity = col_skip(), 
                                             identityRatio = col_number(), qCover = col_number(), 
                                             tCover = col_number()))
#identicos
##needle ##
grupo1 <- hsap_vs_mmus %>% filter(hsap_vs_mmus$identityRatio >= 1)

grupo1a <- hsap_vs_mmus %>% filter(hsap_vs_mmus$qCover == 100
                                   & hsap_vs_mmus$tCover == 100)

## rbh ## 
# identidade menor que 100%
grupo2 <- hsap_vs_mmus %>% filter(hsap_vs_mmus$identityRatio < 1)

grupo2a <- hsap_vs_mmus %>% filter(hsap_vs_mmus$qCover < 100
                        & hsap_vs_mmus$tCover < 100)


#### lendo arquivo hsap vs rnov ####
hsap_vs_rnov <- read_table2("B:/Doutorado/ortologos_flavia/hsap-vs-rnov/spliceprot_FinalTab_App1.txt", 
                            col_types = cols(similarity = col_skip(), 
                                             similarityRatio = col_number(), identity = col_skip(), 
                                             identityRatio = col_number(), qCover = col_number(), 
                                             tCover = col_number()))

#identicos
grupo3 <- hsap_vs_rnov %>% filter(hsap_vs_rnov$identityRatio >= 1)

grupo3a <- hsap_vs_rnov %>% filter(hsap_vs_rnov$qCover == 100
                                   & hsap_vs_rnov$tCover == 100)

# identidade menor que 100%
grupo4 <- hsap_vs_rnov %>% filter(hsap_vs_rnov$identityRatio < 1)

grupo4a <- hsap_vs_rnov %>% filter(hsap_vs_rnov$qCover < 100
                                   & hsap_vs_rnov$tCover < 100)

#### lendo arquivo rnov_vs_mmus ####
rnov_vs_mmus <- read_table2("B:/Doutorado/ortologos_flavia/rnov-vs-mmus/analysis/finalTab/spliceprot_FinalTab_App1.txt", 
                            col_types = cols(similarity = col_skip(), 
                                             similarityRatio = col_number(), identity = col_skip(), 
                                             identityRatio = col_number(), qCover = col_number(), 
                                             tCover = col_number()))


#identicos
grupo5 <- rnov_vs_mmus %>% filter(rnov_vs_mmus$identityRatio >= 1)

grupo5a <- rnov_vs_mmus %>% filter(rnov_vs_mmus$qCover == 100
                                   & rnov_vs_mmus$tCover == 100)

# identidade menor que 100%
grupo6 <- rnov_vs_mmus %>% filter(rnov_vs_mmus$identityRatio < 1)

grupo6a <- rnov_vs_mmus %>% filter(rnov_vs_mmus$qCover < 100
                                   & rnov_vs_mmus$tCover < 100)

#### needle vs rbh ####

#### humano genes ####
human_genes_needle <- data.frame(hsap_vs_mmus$hsTranscriptID)
human_genes_needle <- gsub("\\-.*","",human_genes_needle$hsap_vs_mmus.hsTranscriptID)
unique_human_needle <- data.frame(unique(human_genes_needle))

human_genes_rbh <- data.frame(hsap_vs_mmus$Query)
human_genes_rbh <- gsub("\\-.*","",human_genes_rbh$hsap_vs_mmus.Query)
unique_human_rbh <- data.frame(unique(human_genes_rbh))

common_genes_human <- data.frame(intersect(unique_human_needle$unique.human_genes_needle.,unique_human_rbh$unique.human_genes_rbh.))

#### humano proteínas hsap x mmus ####
# needle 
human_proteins_needle <- data.frame(hsap_vs_mmus$hsTranscriptID)
unique_human_proteins_needle <- data.frame(unique(human_proteins_needle))

#rbh 
human_proteins_rbh <- data.frame(hsap_vs_mmus$Query)
unique_human_proteins_rbh <- data.frame(unique(human_proteins_rbh))

common_proteins_human <- data.frame(intersect(unique_human_proteins_needle$hsap_vs_mmus.hsTranscriptID, 
                                   unique_human_proteins_rbh$hsap_vs_mmus.Query)) 

#### humano proteínas hsap x rnov ####

human_proteins_needle2 <- data.frame(hsap_vs_rnov$hsTranscriptID)
unique_human_proteins_needle2 <- data.frame(unique(human_proteins_needle2))

#rbh 
human_proteins_rbh2 <- data.frame(hsap_vs_rnov$Query)
unique_human_proteins_rbh2 <- data.frame(unique(human_proteins_rbh2))

common_proteins_human2 <- data.frame(intersect(unique_human_proteins_needle2$hsap_vs_rnov.hsTranscriptID, 
                                              unique_human_proteins_rbh2$hsap_vs_rnov.Query)) 


#### camundongo ####
mouse_genes <- data.frame(hsap_vs_mmus$mmTranscriptID)
mouse_genes <- gsub("\\-.*","",mouse_genes$hsap_vs_mmus.hsTranscriptID)
unique_mouse <- data.frame(unique(mouse_genes))

#### rato ####
rat_genes <- data.frame(hsap_vs_rnov$rnTranscriptID)
rat_genes <- gsub("\\-.*","",rat_genes$rnov_vs_mmus.rnTranscriptID)
unique_rat <- data.frame(unique(rat_genes))
