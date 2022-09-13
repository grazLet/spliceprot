import sys
import os
import itertools
import pickle

TableOrtho = open("/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/data/ortho_ensembl_101/Table_ID_Hsap_Mmus_Rnov_Ortho.tsv").read().split("\n")[1:-1] #leticia
orthodict2 = {}
OrthoDict={}
for line in TableOrtho:
    col = line.split()
    #  considerar apenas ortologo 1:1
    mmusID2 = col[1]
    #mmusID2 = col[1]
    rnovID2 = col[2]
#    pair = mmusID+"-"+mmusID
    if "," not in rnovID2:
        if rnovID2 in orthodict2.keys():
            orthodict2[rnovID2].append(mmusID2)
        if rnovID2 not in orthodict2.keys():
            orthodict2[rnovID2] = [mmusID2]

OrthoDictPrint = {}
for rnovID,mmusID3 in orthodict2.items():
    if len(mmusID3) == 1:
        mmusID = "".join(mmusID3)
        #pair = mmusID+"-"+mmusID
        pair = mmusID+"-"+rnovID
        OrthoDict[pair] = {"mmensT":[], "rnensT":[],
                            "mmensP":[], "rnensP":[],
                            #"mmccds":[], "rnccds":[],
                            "mmSp":[], "rnSp":[],
                            "mmSpT":[], "rnSpT":[],
                            "mmrefseq":[],"rnrefseq":[]} 
        OrthoDictPrint[mmusID] = rnovID

#pickleFileName = "/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/data/ortho_ensembl_101/Table_ID_Hsap_Rnov_Ortho_1-1.pickle"
pickleFileName = "/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/data/ortho_ensembl_101/Table_ID_Rnov_Mmus_Ortho_1-1.pickle"
pickleFileout =  open(pickleFileName, "wb")
pickle.dump(OrthoDictPrint,pickleFileout)
pickleFileout.close()

#hsEnsemblFileName = "/home/flavia.freitas/disk1/flavia.freitas/data/ensembl/Hsap_Gene_Proteins_Seq.tsv"
#mmEnsemblFileName = "/home/flavia.freitas/disk1/flavia.freitas/data/ensembl/Mmus_Gene_Proteins_Seq.tsv"

#hsEnsemblFileName = "/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/data/ensembl/hsap/Homo_sapiens.GRCh38.pep.all.fa"
mmEnsemblFileName = "/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/data/ensembl/mmus/Mus_musculus.GRCm38.pep.all.fa"
rnEnsemblFileName = "/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/data/ensembl/rnov/Rattus_norvegicus.Rnor_6.0.pep.all.fa"

ensemblFiles = [mmEnsemblFileName,rnEnsemblFileName]
ensDict = {}
EnsProtDict = {}
trEnsDict = {}
for ensFileName in ensemblFiles:
    ensFile = open(ensFileName).read().split(">")[1:]
    for fasta in ensFile:
        header,seq = fasta.split("\n",1)
        seq = seq.replace("\n","")
        header_itens = header.split(" ")
        protID,lixo = header_itens[0].split(".")
        geneID,lixo = header_itens[3].replace("gene:","").split(".")
        trID,lixo = header_itens[4].replace("transcript:","").split(".")
        trEnsDict[trID] = geneID
        if len(seq) >= 20:
            EnsProtDict[protID] = {"seq":seq,"trID":trID, "geneID":geneID}
        if geneID in ensDict.keys():
            ensDict[geneID]["TrID"].append(trID)
            ensDict[geneID]["ProtID"].append(protID)
        if geneID not in ensDict.keys():
            ensDict[geneID] = {"TrID": [trID],"ProtID":[protID]}
        
#hsGene2Ens = "/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/data/gene_info/hsap_gene2ensembl.txt"
mmGene2Ens = "/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/data/gene_info/mmus_gene2ensembl.txt"
rnGene2Ens = "/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/data/gene_info/rnov_gene2ensembl.txt"
#hsGene2RefSeq = "/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/data/gene_info/hsap_gene2refseq.txt"
mmGene2RefSeq = "/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/data/gene_info/mmus_gene2refseq.txt"
rnGene2RefSeq = "/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/data/gene_info/rnov_gene2refseq.txt" 

Gene2EnsFiles = [mmGene2Ens,rnGene2Ens]
Gene2RefSeqFiles = [mmGene2RefSeq,rnGene2RefSeq]
refSeqDict = {}
geneidD = {}

for Gene2EnsFileName in Gene2EnsFiles:
    Gene2EnsFile = open(Gene2EnsFileName).read().split("\n")[1:-1]
    #alimentar dicitonario com geneID e ensemblID
    for line2 in Gene2EnsFile:
        col2 = line2.split()
        gID3 = col2[1].strip()
        ensgID = col2[2].strip()
        if gID3 not in geneidD.keys():
            geneidD[gID3] = ensgID


for Gene2RefSeqFileName in Gene2RefSeqFiles:
    Gene2RefSeqFile = open(Gene2RefSeqFileName).read().split("\n")[1:-1]
    for line in Gene2RefSeqFile:
        col = line.split()
        gID2 = col[1].strip()
        if gID2 in geneidD.keys():
            gID = geneidD[gID2]
            rnaID = col[3]
            prID = col[5]
            if rnaID != "-" and prID != "-":
                refseqTrID,lixo = rnaID.split(".")
                refseqPrID,lixo = prID.split(".")

                if gID in refSeqDict.keys():
                    refSeqDict[gID]["refseqTrID"].append(refseqTrID)
                    refSeqDict[gID]["refseqPrID"].append(refseqPrID)
                
                if gID not in refSeqDict.keys():
                    refSeqDict[gID] = {"refseqTrID":[refseqTrID],
                                    "refseqPrID":[refseqPrID]}
            
#hsCCDSFileName = "/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/data/ccds/Hs.CCDS2Sequence.current.txt"
#rnCCDSFileName = "/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/ccds/Rn.CCDS2Sequence.current.txt"

""" CCDSFileNames = [hsCCDSFileName,rnCCDSFileName]#
ccdsDict = {}
for CCDSFileName in CCDSFileNames:
    CCDSFile = open(CCDSFileName).read().split("\n")[1:-1]
    for line3 in CCDSFile:
        col3 = line3.split()
        if col3[3] == "EBI":
            ccdsID,lixo5 = col3[0].split(".")
            ensTr2,lixo6 = col3[4].split(".")
            if ccdsID in ccdsDict.keys():
                ccdsDict[ccdsID].append(ensTr2)
            if ccdsID not in ccdsDict.keys():
                ccdsDict[ccdsID] = [ensTr2]

ccdsGeneDict = {}
for ccds,trList in ccdsDict.items():
    for trID in trList:
        if trID in trEnsDict.keys():
            geneID = trEnsDict[trID]
            if geneID in ccdsGeneDict.keys():
                ccdsGeneDict[geneID].append(ccds)
            if geneID not in ccdsGeneDict.keys():
                ccdsGeneDict[geneID] = [ccds]
 """

#hsSpFASTA = open("/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/data/spliceprot_r100_v2/Hs.traducao_CORRIGIDA_release100.fasta").read().split(">")[1:]
rnSpFASTA = open("/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/data/spliceprot_r100_v2/Rn.traducao_CORRIGIDA_release100.fasta").read().split(">")[1:]
mmSpFASTA = open("/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/data/spliceprot_r100_v2/Mm.traducao_CORRIGIDA_release100.fasta").read().split(">")[1:]

MMspProtDict = {}
RNspProtDict = {}

for mmfasta in mmSpFASTA:
    mmHeader1,mmSeq = mmfasta.split("\n",1)
    MmSeq = mmSeq.strip()
    mmHeader2 = mmHeader1.split("_",1)
    mmHeader3 = mmHeader2[0].replace("Mm.","")
    #correcao do ID do fasta
    # modelo ENSMUSG00000100280
    # len("ENSMUSG00000100280") = 18
    mmIdLen = len(mmHeader3)
    mmAdd = 11 - mmIdLen
    aa = []
    for a1 in range(mmAdd):aa.append('0')
    mmAdd2 = "".join(aa)
    MmgeneID = "ENSMUSG"+mmAdd2+mmHeader3
    MmProtID2,MmEnTrID2 = mmHeader2[1].split(" ")
    MmProtID4 = MmProtID2.replace(":","").replace("__","_") #MmHeader2[1].strip()
    MmProtID = MmgeneID+"-"+MmProtID4
    MmEnTrID3 = MmEnTrID2.replace("[","").split("]")[:-1]
    MmEnTrID = "_".join(MmEnTrID3)
    if MmgeneID in MMspProtDict.keys():
        MMspProtDict[MmgeneID]["MmProtID"].append(MmProtID)
        MMspProtDict[MmgeneID]["MmSeq"].append(MmSeq)
        MMspProtDict[MmgeneID]["MmensTr"].append(MmEnTrID)

    if MmgeneID not in MMspProtDict.keys():
        MMspProtDict[MmgeneID] = {"MmProtID":[MmProtID],"MmSeq":[MmSeq],"MmensTr":[MmEnTrID]}

for rnfasta in rnSpFASTA:
    rnHeader1,rnSeq = rnfasta.split("\n",1)
    RnSeq = rnSeq.strip()
    rnHeader2 = rnHeader1.split("_",1)
    rnHeader3 = rnHeader2[0].replace("Rn.","")
    #correcao do ID do fasta
    # modelo ENSRNOG00000008786 
    # len("ENSRNOG00000008786") = 18
    rnIdLen = len(rnHeader3)
    rnAdd = 11 - rnIdLen
    bb = []
    for b1 in range(rnAdd):
        bb.append('0')
    rnAdd2 = "".join(bb)
    RngeneID = "ENSRNOG"+rnAdd2+rnHeader3
    RnProtID2,RnEnTrID2 = rnHeader2[1].split(" ")
    RnProtID4 = RnProtID2.replace(":","").replace("__","_") #rnHeader2[1].strip()
    RnProtID = RngeneID+"-"+RnProtID4
    RnEnTrID3 = RnEnTrID2.replace("[","").split("]")[:-1]
    RnEnTrID = "_".join(RnEnTrID3)
    if RngeneID in RNspProtDict.keys():
        RNspProtDict[RngeneID]["RnProtID"].append(RnProtID)
        RNspProtDict[RngeneID]["RnSeq"].append(RnSeq)
        RNspProtDict[RngeneID]["RnensTr"].append(RnEnTrID)
    if RngeneID not in RNspProtDict.keys():
        RNspProtDict[RngeneID] = {"RnProtID": [RnProtID], 
                                  "RnSeq": [RnSeq],
                                  "RnensTr": [RnEnTrID]}                                  


for ortho,info in OrthoDict.items():
    mmID, rnID = ortho.split("-")
    if mmID in ensDict.keys():
        info["mmensT"] = ensDict[mmID]["TrID"]
        info["mmensP"] = ensDict[mmID]["ProtID"]
    if mmID in refSeqDict.keys():
        info["mmrefseq"] = refSeqDict[mmID]["refseqTrID"]
#        info["mmensT"] = geneDict[mmID]["ensTrID"]
 #       info["mmensP"] = geneDict[mmID]["ensPrID"]
    """ if mmID in ccdsGeneDict.keys():
        info["mmccds"] = ccdsGeneDict[mmID] """
    if mmID in MMspProtDict.keys():
        info["mmSp"] = MMspProtDict[mmID]["MmProtID"]
        info["mmSpT"] = MMspProtDict[mmID]["MmensTr"]
    if rnID in ensDict.keys():
        info["rnensT"] = ensDict[rnID]["TrID"]
        info["rnensP"] = ensDict[rnID]["ProtID"]
    if rnID in refSeqDict.keys():
        info["rnrefseq"] = refSeqDict[rnID]["refseqTrID"]
    #    info["rnensT"] = geneDict[rnID]["ensTrID"]
     #   info["rnensP"] = geneDict[rnID]["ensPrID"]
    #if rnID in ccdsGeneDict.keys():
        #info["rnccds"] = ccdsGeneDict[rnID]
    if rnID in RNspProtDict.keys():
        info["rnSp"] = RNspProtDict[rnID]["RnProtID"]
        info["rnSpT"] = RNspProtDict[rnID]["RnensTr"]


#for a,b in OrthoDict.items():
 #   print(a,b)

#hsCCDSFASTAfileName = "/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/data/ccds/Hsap_CCDS_protein.current.faa"
#rnCCDSFASTAfileName = "/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/data/ccds/Rnov_CCDS_protein.current.faa"
#ccdsFASTAfileNames = [hsCCDSFASTAfileName,rnCCDSFASTAfileName]

#hsCCDSFasta = open(hsCCDSFASTAfileName).read().split(">")[1:]
#mmCCDSFasta = open(mmCCDSFASTAfileName).read().split(">")[1:]
#ccdsFastaD = {}
""" for ccdsFastaFileName in ccdsFASTAfileNames:
    CCDSFastaFile = open(ccdsFastaFileName).read().split(">")[1:]   
    for ccdsFasta in CCDSFastaFile:
        h,s = ccdsFasta.split("\n",1)
        seq = s.replace("\n","")
        h2 = h.split("|")
        h3,lixo = h2[0].split(".")
        if len(seq) >= 20:
            ccdsFastaD[h3] = seq """

#os.mkdir("ccds/")
#os.mkdir("ccds/needle/") 
#os.mkdir("ccds/needle/HsapRnov_OrthoProt")
#os.mkdir("ccds/needle/HsapRnov_FASTA")
#hsFileOut = open("ccds/needle/HsapRnov_FASTA/hsap.faa","w")
#rnFileOut = open("ccds/needle/HsapRnov_FASTA/rnov.faa","w") 

""" for ortho2,infoID in OrthoDict.items():
    hsap,rnov = ortho2.split("-")
    mouse = ""
    mouse = ""
    rat = ""
    hsPrint = []
    rnPrint = []
    if infoID["hsccds"] != [] and infoID["rnccds"] != []:
        hsCCDSlist = set(infoID["hsccds"])
        rnCCDSlist = set(infoID["rnccds"])
        for ccds1 in hsCCDSlist:
            if ccds1 in ccdsFastaD.keys():
                sequence1 = ccdsFastaD[ccds1]
                header = "|".join(["Hsap",hsap,ccds1])
                header1 = ">"+header+"\n"
                mouse = hsap
                hsPrint.append((header1,sequence1))
        for ccds2 in rnCCDSlist:
            if ccds2 in ccdsFastaD.keys():
                sequence2 = ccdsFastaD[ccds2]
                header2 = "|".join(["Rnov",rnov,ccds2])
                header12 = ">"+header2+"\n"
                rat = rnov
                rnPrint.append((header12,sequence2))
        if mouse != "" and rat != "":
            fileOutName = "ccds/needle/HsapRnov_OrthoProt/"+mouse+"_"+rat+"_ccds_ortho_proteins.fa"
            fileOut = open(fileOutName,"w")
            
            for i in range(len(hsPrint)):
                 HeaderHs = hsPrint[i][0]
                 SequenceHs = hsPrint[i][1]
                 fileOut.write(HeaderHs)
                 fileOut.write(SequenceHs)
                 fileOut.write("\n")
                 hsFileOut.write(HeaderHs)
                 hsFileOut.write(SequenceHs)
                 hsFileOut.write("\n")

            for j in range(len(rnPrint)):
                 HeaderRn = rnPrint[j][0]
                 SequenceRn = rnPrint[j][1]
                 fileOut.write(HeaderRn)
                 fileOut.write(SequenceRn)
                 fileOut.write("\n")
                 rnFileOut.write(HeaderRn)
                 rnFileOut.write(SequenceRn)
                 rnFileOut.write("\n")
        fileOut.close()
hsFileOut.close()
rnFileOut.close() """


#os.mkdir("refseq/")
#os.mkdir("refseq/needle/")
os.mkdir("refseq/needle/MmusRnov_OrthoProt")
os.mkdir("refseq/needle/MmusRnov_FASTA")
#hsFileOut2 = open("/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/final_pipeline/refseq/needle/HsapRnov_FASTA/hsap.faa","w")
rnFileOut2 = open("/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/final_pipeline/refseq/needle/MmusRnov_FASTA/rnov.faa","w")
mmFileOut2 = open("/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/final_pipeline/refseq/needle/MmusRnov_FASTA/mmus.faa","w")

refseqFASTAfileNames = ["/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/data/refseq/Mmus_refseq_prot.faa","/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/data/refseq/Rnov_refseq_prot.faa"]

refseqFastaD = {}
for refseqFastaFileName in refseqFASTAfileNames:
    refseqFastaFile = open(refseqFastaFileName).read().split("\n>")[1:]

    for refseqFasta in refseqFastaFile:
        h,s = refseqFasta.split("\n",1)
        seq = s.replace("\n","")
        h2 = h.split(" ")
        h3,lixo = h2[0].split(".")
        if len(seq) >= 20:
            refseqFastaD[h3] = seq

for ortho2,infoID in OrthoDict.items():
    mmus,rnov = ortho2.split("-")
    mouse = ""
    rat = ""
    mmPrint = []
    rnPrint = []
    if mmus in refSeqDict.keys() and rnov in refSeqDict.keys():
        mmTrList = refSeqDict[mmus]["refseqTrID"]
        mmPrList = refSeqDict[mmus]["refseqPrID"]
        rnTrList = refSeqDict[rnov]["refseqTrID"]
        rnPrList = refSeqDict[rnov]["refseqPrID"]

    if infoID["mmrefseq"] != [] and infoID["rnrefseq"] != []:
        mmREFSEQlist = set(infoID["mmrefseq"])
        rnREFSEQlist = set(infoID["rnrefseq"])

        for refseq1 in mmREFSEQlist:
            if refseq1 in mmTrList:
                rstrIndex = mmTrList.index(refseq1)
                refseqPr1 = mmPrList[rstrIndex]
                if refseqPr1 in refseqFastaD.keys():
                    sequence1 = refseqFastaD[refseqPr1]
                    header = "|".join(["Mmus",mmus,refseqPr1,refseq1])
                    header1 = ">"+header+"\n"
                    mouse = mmus 
                    mmPrint.append((header1,sequence1))
        for refseq2 in rnREFSEQlist:
            if refseq2 in rnTrList:
                rstrIndex2 = rnTrList.index(refseq2)
                refseqPr2 = rnPrList[rstrIndex2]
                if refseqPr2 in refseqFastaD.keys():
                    sequence2 = refseqFastaD[refseqPr2]
                    header2 = "|".join(["Rnov",rnov,refseqPr2,refseq2])
                    header12 = ">"+header2+"\n"
                    rat = rnov
                    rnPrint.append((header12,sequence2))
        if mouse != "" and rat != "":
            fileOutName = "/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/final_pipeline/refseq/needle/MmusRnov_OrthoProt/"+mouse+"_"+rat+"_refseq_ortho_proteins.fa"
            fileOut = open(fileOutName,"w")
            for i in range(len(mmPrint)):
                HeaderMm = mmPrint[i][0]
                SequenceMm = mmPrint[i][1]
                fileOut.write(HeaderMm)
                fileOut.write(SequenceMm)
                fileOut.write("\n")
                mmFileOut2.write(HeaderMm)
                mmFileOut2.write(SequenceMm)
                mmFileOut2.write("\n")

            for j in range(len(rnPrint)):
                HeaderRn = rnPrint[j][0]
                SequenceRn = rnPrint[j][1]
                fileOut.write(HeaderRn)
                fileOut.write(SequenceRn)
                fileOut.write("\n")
                rnFileOut2.write(HeaderRn)
                rnFileOut2.write(SequenceRn)
                rnFileOut2.write("\n")
        fileOut.close()
mmFileOut2.close()
rnFileOut2.close() 



#os.mkdir("spliceprot/")
#os.mkdir("spliceprot/needle/")
os.mkdir("/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/final_pipeline/spliceprot/needle/MmusRnov_OrthoProt")
os.mkdir("/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/final_pipeline/spliceprot/needle/MmusRnov_FASTA")
#hsFileOut3 = open("/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/final_pipeline/spliceprot/needle/MmusRnov_FASTA/hsap.faa","w")
mmFileOut3 = open("/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/final_pipeline/spliceprot/needle/MmusRnov_FASTA/mmus.faa","w")
rnFileOut3 = open("/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/final_pipeline/spliceprot/needle/MmusRnov_FASTA/rnov.faa","w")

SpFASTAFileNames = ["/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/data/spliceprot_r100_v2/Mm.traducao_CORRIGIDA_release100.fasta","/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/data/spliceprot_r100_v2/Rn.traducao_CORRIGIDA_release100.fasta"]
SpFastaD = {}
for SpFASTAFilename in SpFASTAFileNames:
    SpFastaFile = open(SpFASTAFilename).read().split(">")[1:]
    for spFASTA in SpFastaFile:
        Header1,Seq = spFASTA.split("\n",1)
        Seq = Seq.strip()
        Header2 = Header1.split("_",1)
        if len(Seq) >= 20:
            if "Mm." in Header2[0]:
                Header3 = Header2[0].replace("Mm.","")
                #correcao do ID do fasta
                # modelo ENSMUSG00000100280
                # len("ENSMUSG00000100280") = 18
                IdLen = len(Header3)
                Add = 11 - IdLen
                aa = []
                for a1 in range(Add):
                    aa.append('0')
                Add2 = "".join(aa)
                geneID = "ENSMUSG"+Add2+Header3
            if "Rn." in Header2[0]:
                Header3 = Header2[0].replace("Rn.","")
                #correcao do ID do fasta
                # modelo ENSRNOG00000008786
                # len("ENSRNOG00000008786") = 18
                IdLen = len(Header3)
                Add = 11 - IdLen
                aa = []
                for a1 in range(Add):
                    aa.append('0')
                Add2 = "".join(aa)
                geneID = "ENSRNOG"+Add2+Header3

            ProtID2,EnTrID2 = Header2[1].split(" ")
            ProtID3 = ProtID2.replace(":","").replace("__","_") 
            ProtID = geneID+"-"+ProtID3 
            EnTrID3 = EnTrID2.replace("[","").split("]")[:-1]
            EnTrID = "_".join(EnTrID3)
            SpFastaD[ProtID] = Seq 

for ortho2,infoID in OrthoDict.items():
    mmus,rnov = ortho2.split("-")
    mouse = ""
    rat = ""
    mmPrint = []
    rnPrint = []
    if infoID["mmSp"] != [] and infoID["rnSp"] != []:
        mmSPlist = set(infoID["mmSp"])
        rnSPlist = set(infoID["rnSp"])
        mmSpliceProtL = infoID["mmSpT"]
        rnSpliceProtL = infoID["rnSpT"]

        for Sp1 in mmSPlist:
            sp1index = infoID["mmSp"].index(Sp1)
            sp1T = mmSpliceProtL[sp1index]

            if Sp1 in SpFastaD.keys():
                sequence1 = SpFastaD[Sp1]
                header = "|".join(["Mmus",Sp1,sp1T])
                header1 = ">"+header+"\n"
                mouse = mmus 
                mmPrint.append((header1,sequence1))
        for Sp2 in rnSPlist:
            sp2index = infoID["rnSp"].index(Sp2)
            sp2T = rnSpliceProtL[sp2index]

            if Sp2 in SpFastaD.keys():
                sequence2 = SpFastaD[Sp2]
                header2 = "|".join(["Rnov",Sp2,sp2T])
                header12 = ">"+header2+"\n"
                rat = rnov
                rnPrint.append((header12,sequence2))
        if mouse != "" and rat != "":
            fileOutName = "/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/final_pipeline/spliceprot/needle/MmusRnov_OrthoProt/"+mouse+"_"+rat+"_ortho_proteins.fa"
            fileOut = open(fileOutName,"w")
            for i in range(len(mmPrint)):
                HeaderMm = mmPrint[i][0]
                SequenceMm = mmPrint[i][1]
                fileOut.write(HeaderMm)
                fileOut.write(SequenceMm)
                fileOut.write("\n")
                mmFileOut3.write(HeaderMm)
                mmFileOut3.write(SequenceMm)
                mmFileOut3.write("\n")

            for j in range(len(rnPrint)):
                HeaderRn = rnPrint[j][0]
                SequenceRn = rnPrint[j][1]
                fileOut.write(HeaderRn)
                fileOut.write(SequenceRn)
                fileOut.write("\n")
                rnFileOut3.write(HeaderRn)
                rnFileOut3.write(SequenceRn)
                rnFileOut3.write("\n")
        fileOut.close()
mmFileOut3.close()
rnFileOut3.close()

#os.mkdir("ensembl/")
#os.mkdir("ensembl/needle/")
#os.mkdir("ensembl/needle/MmusRnov_OrthoProt")
#os.mkdir("ensembl/needle/MmusRnov_FASTA")
#hsFileOut4 = open("/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/final_pipeline/ensembl/needle/HsapRnov_FASTA/hsap.faa","w")
mmFileOut4 = open("/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/final_pipeline/ensembl/needle/MmusRnov_FASTA/mmus.faa","w")
rnFileOut4 = open("/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/final_pipeline/ensembl/needle/MmusRnov_FASTA/rnov.faa","w")

for ortho2,infoID in OrthoDict.items():
    mmus,rnov = ortho2.split("-")
    mouse = ""
    rat = ""
    mmPrint = []
    rnPrint = []
    if infoID["mmensP"] != [] and infoID["rnensP"] != []:
        mmEnslist = set(infoID["mmensP"])
        rnEnslist = set(infoID["rnensP"])
        for Ens1 in mmEnslist:
            if Ens1 in EnsProtDict.keys():
                sequence1 = EnsProtDict[Ens1]["seq"]
                TrID = EnsProtDict[Ens1]["trID"]
                header = "|".join(["Mmus",mmus,TrID,Ens1])
                header1 = ">"+header+"\n"
                mouse = mmus
                mmPrint.append((header1,sequence1))
        for Ens2 in rnEnslist:
            if Ens2 in EnsProtDict.keys():
                sequence2 = EnsProtDict[Ens2]["seq"]
                TrID2 = EnsProtDict[Ens2]["trID"]
                header2 = "|".join(["Rnov",rnov,TrID2,Ens2])
                header12 = ">"+header2+"\n"
                rat = rnov
                rnPrint.append((header12,sequence2))
        if mouse != "" and rat != "":
            fileOutName = "/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/final_pipeline/ensembl/needle/MmusRnov_OrthoProt/"+mouse+"_"+rat+"_ortho_proteins.fa"
            fileOut = open(fileOutName,"w")
            for i in range(len(mmPrint)):
                HeaderMm = mmPrint[i][0]
                SequenceMm = mmPrint[i][1]
                fileOut.write(HeaderMm)
                fileOut.write(SequenceMm)
                fileOut.write("\n")
                mmFileOut4.write(HeaderMm)
                mmFileOut4.write(SequenceMm)
                mmFileOut4.write("\n")
            for j in range(len(rnPrint)):
                HeaderRn = rnPrint[j][0]
                SequenceRn = rnPrint[j][1]
                fileOut.write(HeaderRn)
                fileOut.write(SequenceRn)
                fileOut.write("\n")
                rnFileOut4.write(HeaderRn)
                rnFileOut4.write(SequenceRn)
                rnFileOut4.write("\n")


        fileOut.close()
mmFileOut4.close()
rnFileOut4.close()