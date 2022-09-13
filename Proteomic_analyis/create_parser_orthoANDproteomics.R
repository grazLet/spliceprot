#datasets
#proteotipicos_spliceprot_final
#dicionario

proteotipicos_spliceprot_final <- read_table2("B:/Doutorado/PatternLab/Outputs_Certos_Maio/SpliceProt_Completo/PXD020656--spliceprot-COMPLETO-SEMNMD-leticia-11052022/proteotipicos-PXD020656-SPLICE-COMPLETO-final.csv", 
                                              col_names = FALSE, col_types = cols(X1 = col_character(), 
                                                                                  X3 = col_number()))
colnames(proteotipicos_spliceprot_final) <- c('TheProteinLoci','Sequence','SpectrumCount')

data <- sqldf('
 select  *
 from    proteotipicos_PXD016793_testerato_final two
         join dicionario one
           on one.HeaderNumerica = two.TheProteinLoci')

data <- data.frame(data$TheProteinLoci,data$Sequence,data$SpectrumCount,data$HeaderOriginal)

colnames(data) <- c('TheProteinLoci','Sequence','SpectrumCount','HeaderOriginal')

write.csv(data,"../../20171209_iBAQ_liver_1000ng_05/../../ortologos_flavia/resultados_finais/PXD016793_lista_proteotipicos_headeroriginal.csv")

#### tentando com as exclusivas do spliceprot ####
exclusivas_spliceprot$`The Protein Loci` <- sub('\\| ','',exclusivas_spliceprot$`The Protein Loci`)

exclusivas_spliceprot$`The Protein Loci` <- sub('\\s+','',exclusivas_spliceprot$`The Protein Loci`)

colnames(exclusivas_spliceprot) <- c('Sequence','TheProteinLoci','ProteinCount','SpectrumCount','MesmaIsoforma')

colnames(dicionario) <- c('HeaderOriginal','HeaderNumerica')

data <- sqldf('
 select  *
 from    exclusivas_spliceprot two
         join dicionario one
           on one.HeaderNumerica = two.TheProteinLoci')

data <- data.frame(data$TheProteinLoci,data$Sequence,data$SpectrumCount,data$HeaderOriginal)

colnames(data) <- c('TheProteinLoci','Sequence','SpectrumCount','HeaderOriginal')

write.csv(data,"../../Doutorado/ortologos_flavia/resultados_finais/PXD016793_lista_exclusivos_spliceprot.csv")
