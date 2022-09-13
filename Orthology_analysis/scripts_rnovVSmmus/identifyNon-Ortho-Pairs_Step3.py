import sys
import os
from Bio.Emboss.Applications import NeedleCommandline
from Bio import AlignIO
import subprocess
import pickle

##### ETAPA 4 #####
# entrada: diretorio q contem os arquivos fastas 
# com as proteoformas de Mmus e Rnov para cada gene
# le automaticamente o dicionario gerado pelo script
# findOrthoTranscripts_Translated-data_3.py
# MmRn-ProteoOrtho_v1.pickle guardado no diretorio ortho-pairs_app1/
# linha comando python scripts/identifyNon-Ortho-Pairs_approach2_4.py MmusRnov_OrthoProt_r100/
# gera nova versao dos arquivos contidos no dir MmusRnov_OrthoProt_r100/
# em q foram removidas as proteoformas cujos pares ja foram encontrados
# na abordagem 1 >> needle, similaridade >= 0.8; gap =< 0.01 e highest score
# no dir FASTA/ serao criados os arquivos que servirao de entrada na abordagem2

# Funcao
def multiFasta_printDict_function(multiFasta,dbName):
    mmDict={}
    rnDict={}
    for fasta in multiFasta:
        header,seq = fasta.split("\n",1)
        seq = seq.replace("\n","")
        header_itens = header.split("|")
        sp = header_itens.pop(0)
        
        if "spliceprot" in dbName:
            transEnsL = header_itens.pop(-1)
            transcriptID = "-".join(header_itens)#header_itens[0]+"_"+header_itens[1]
        # if "ccds" in dbName:
        #     transcriptID = "-".join(header_itens)#header_itens[1]+"_"+header_itens[0]
        if "refseq" in dbName:
            transcriptID = "-".join(header_itens)#header_itens[0]+"_"+header_itens[1]
        if "ensembl" in dbName:
            transcriptID = "-".join(header_itens)#header_itens[0]+"_"+header_itens[1]

        if "Mmus" in sp:
            fasta2 = ">Mmus_"+transcriptID+"\n"+seq+"\n" #">Mmus_"+transcriptID+"\n"+seq+"\n"
            mmDict[transcriptID] = {"fasta":fasta2}
                        
        if "Rnov" in sp:
            fasta2 = ">Rnov_"+transcriptID+"\n"+seq+"\n" #">Rnov_"+transcriptID+"\n"+seq+"\n"
            rnDict[transcriptID] = {"fasta":fasta2}

    return (mmDict,rnDict) 


folder = sys.argv[1]
folder2 = folder.replace("MmusRnov_OrthoProt/","")
folder3 = folder.replace("needle/MmusRnov_OrthoProt/","")
folder_itens = folder.split("/")
database = folder_itens[0]

PickleFile = open(folder2+"ortho-pairs_app1/MmRn-ProteoOrtho_v1.pickle","rb")
ProtOrtDict = pickle.load(PickleFile)
# ProtoOrtDict contem Mm-Proteoforms como key
# e Rn-Proteoforms como value
# key, value  = pares de proteoformas ortologas
#protOrtMm = ProtOrtDict.keys()
#protOrtRn = ProtOrtDict.values()
protOrtMm = []
protOrtRn = []
for key,value in ProtOrtDict.items():
    if key not in protOrtMm:
        protOrtMm.append(key)
    if type(value) is list:
        value_list = value
        for value2 in value_list:
            if value2 not in protOrtRn:
                protOrtRn.append(value2)
    if type(value) is not list:
        value2 = value
        if value2 not in protOrtRn:
            protOrtRn.append(value2)

# passar diretorio onde estao os dados q devem ser analisados
# dados = arquivos fasta com os transcritos dos genes ortologos
# das duas ou tres especies analisadas
# rodar ls comando


out = subprocess.Popen(['ls', folder],
           stdout = subprocess.PIPE,
           stderr = subprocess.STDOUT, encoding='utf-8')
stdout2,stderr2 = out.communicate()

fileList = stdout2.split()
#os.mkdir(folder2+"non-ortho-app2/")
#os.mkdir(folder2+"non-ortho-app2/FASTA/")
# criar diretorio para rodar o diamond

#os.mkdir(folder3+"diamond_needle")
#os.mkdir(folder3+"diamond_needle/data")
# Criar arquivos com todas as isoformas de todos os genes
# para cada especie
# dado de entrada do orthofinder
FileMmusNonOrtho2 = open(folder3+"diamond_needle/data/mmus2.faa","w")
FileRnovNonOrtho2 = open(folder3+"diamond_needle/data/rnov2.faa","w")

for files2 in fileList:
        files = folder+files2
        multiFasta2 = open(files).read().split(">")[1:]
        fileOutName = folder2+"non-ortho-app2/FASTA/"+files2.replace("_ortho_proteins.fasta","_non-ortho-app1.fa")
        fileOut = open(fileOutName,"w")
        multiFasta3 = multiFasta_printDict_function(multiFasta2,database)
        mmD = multiFasta3[0]
        rnD = multiFasta3[1]

        for mm1, mm2 in mmD.items():
                if mm1 not in protOrtMm:
                        fileOut.write(mm2["fasta"])
                        #FileMmusNonOrtho.write(mm2["fasta"])
        for rn1,rn2 in rnD.items():
                if rn1 not in protOrtRn:
                        fileOut.write(rn2["fasta"])
                        #FileRnovNonOrtho.write(rn2["fasta"])
    # separar seqs q nao estao presentes dict ProtOrtDict (pares de proteoformas)
    # criar novos arquivos
    # para fazer nova abordagem --> approach 2
    

fileOut.close()

#FileMmusNonOrtho.close()
#FileRnovNonOrtho.close()

filesFastaDir = os.listdir(folder2+"non-ortho-app2/FASTA/")
for eachFile in filesFastaDir:
    pathFile = "%s%s" % (folder2+"non-ortho-app2/FASTA/",eachFile)
    readPathFile = open(pathFile).read()
    if os.path.getsize(pathFile) == 0:
        os.remove(pathFile)
    
    else:
        if "Rnov" not in readPathFile:
            os.remove(pathFile)
        if "Mmus" not in readPathFile:
            os.remove(pathFile)

filesFastaDir2 = os.listdir(folder2+"non-ortho-app2/FASTA/")
for eachGeneFile in filesFastaDir2:
        pathFile2 = "%s%s" % (folder2+"non-ortho-app2/FASTA/",eachGeneFile)
        readFastaFile = open(pathFile2).read().split(">")[1:]
        for FASTA1 in readFastaFile:
                if "Mmus" in FASTA1:
                        #FASTA1 = FASTA1.replace("Mmus_","")
                        FASTA2 = "%s%s" % (">",FASTA1)
                        FileMmusNonOrtho2.write(FASTA2)
                if "Rnov" in FASTA1:
                        #FASTA1 = FASTA1.replace("Rnov_","")
                        FASTA2 = "%s%s" % (">",FASTA1)
                        FileRnovNonOrtho2.write(FASTA2)


FileMmusNonOrtho2.close()
FileRnovNonOrtho2.close()