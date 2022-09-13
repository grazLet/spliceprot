import sys
import os

featuresList = ["CDS","exon","mRNA","gene","three_prime_UTR","five_prime_UTR"]
 
hsGFFfile = open("/home/flavia.freitas/disk1/flavia.freitas/data/ensembl/hsap/Homo_sapiens.GRCh38.100.gff3").read().split("###")[201:-1]
mmGFFfile = open("/home/flavia.freitas/disk1/flavia.freitas/data/ensembl/hsap/Homo_sapiens.GRCh38.100.gff3").read().split("\n")[:-1]

def parseGFF(GFFfile):
    featuresList = ["CDS","exon","mRNA","gene","three_prime_UTR","five_prime_UTR"]
    for entry in GFFfile:
        lines = entry.split("\n")[:-1]
        for line in lines:
            col = line.split()
            feature = col[2].strip()
            if feature in featuresList:
                if feature == "gene":
                    loc_feature = GFFfile.index(line)
                    info = col[8].split(";")
                    geneID = info[0].replace("ID=gene:","")
                    return(">>>",entry)
                #if feature != "gene" and geneID in line:
                    if feature == "mRNA":
                        trIDinfo = col[8].split(";")
                        trID = trIDinfo[0].replace("ID=transcript:","")
                    if geneID in gffD.keys():
                        if trID in gffD[geneID][trID]:
                            gffD[geneID][trID].append(col)
                        if trID not in gffD[geneID][trID]:
                            gffD[geneID][trID] = [col]
                    if geneID not in gffD.keys():
                        gffD[geneID][trID] = [col]
 #   return(gffD)

hsD = parseGFF(hsGFFfile)
#mmD = parseGFF(mmGFFfile)
print(hsD)
TabFinal = open("/home/flavia.freitas/disk1/flavia.freitas/results/proteoformas/analysis/finalTab/ensembl_FinalTab_App1.txt").read().split("\n")[:-1]
header = TabFinal.pop(0)

def parseTabFinal(TabFinalFile):
    for line in TabFinalFile:
        col = line.split()
        hsG = col[0]
        mmG = col[1]
        method = col[-1]
        if "needle" in method:
            hsT = col[2]
            mmT = col[3]
        if "rbh" in method:
            hsT = col[11]
            mmT = col[12]
        
         

