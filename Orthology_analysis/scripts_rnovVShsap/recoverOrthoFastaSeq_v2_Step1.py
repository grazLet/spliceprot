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
    hsapID2 = col[0]
    #mmusID2 = col[1]
    rnovID2 = col[2]
#    pair = hsapID+"-"+mmusID
    if "," not in rnovID2:
        if rnovID2 in orthodict2.keys():
            orthodict2[rnovID2].append(hsapID2)
        if rnovID2 not in orthodict2.keys():
            orthodict2[rnovID2] = [hsapID2]

OrthoDictPrint = {}
for rnovID,hsapID3 in orthodict2.items():
    if len(hsapID3) == 1:
        hsapID = "".join(hsapID3)
        #pair = hsapID+"-"+mmusID
        pair = hsapID+"-"+rnovID
        OrthoDict[pair] = {"hsensT":[], "rnensT":[],
                            "hsensP":[], "rnensP":[],
                            #"hsccds":[], "rnccds":[],
                            "hsSp":[], "rnSp":[],
                            "hsSpT":[], "rnSpT":[],
                            "hsrefseq":[],"rnrefseq":[]} 
        OrthoDictPrint[hsapID] = rnovID

pickleFileName = "/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/data/ortho_ensembl_101/Table_ID_Hsap_Rnov_Ortho_1-1.pickle"
pickleFileout =  open(pickleFileName, "wb")
pickle.dump(OrthoDictPrint,pickleFileout)
pickleFileout.close()

#hsEnsemblFileName = "/home/flavia.freitas/disk1/flavia.freitas/data/ensembl/Hsap_Gene_Proteins_Seq.tsv"
#mmEnsemblFileName = "/home/flavia.freitas/disk1/flavia.freitas/data/ensembl/Mmus_Gene_Proteins_Seq.tsv"

hsEnsemblFileName = "/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/data/ensembl/hsap/Homo_sapiens.GRCh38.pep.all.fa"
rnEnsemblFileName = "/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/data/ensembl/rnov/Rattus_norvegicus.Rnor_6.0.pep.all.fa"

ensemblFiles = [hsEnsemblFileName,rnEnsemblFileName]
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
        
hsGene2Ens = "/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/data/gene_info/hsap_gene2ensembl.txt"
rnGene2Ens = "/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/data/gene_info/rnov_gene2ensembl.txt"
hsGene2RefSeq = "/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/data/gene_info/hsap_gene2refseq.txt"
rnGene2RefSeq = "/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/data/gene_info/rnov_gene2refseq.txt" 

Gene2EnsFiles = [hsGene2Ens,rnGene2Ens]
Gene2RefSeqFiles = [hsGene2RefSeq,rnGene2RefSeq]
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

hsSpFASTA = open("/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/data/spliceprot_r100_v2/Hs.traducao_CORRIGIDA_release100.fasta").read().split(">")[1:]
rnSpFASTA = open("/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/data/spliceprot_r100_v2/Rn.traducao_CORRIGIDA_release100.fasta").read().split(">")[1:]

HSspProtDict = {}
RNspProtDict = {}

for hsfasta in hsSpFASTA:
    hsHeader1,hsSeq = hsfasta.split("\n",1)
    HsSeq = hsSeq.strip()
    hsHeader2 = hsHeader1.split("_",1)
    hsHeader3 = hsHeader2[0].replace("Hs.","")
    #correcao do ID do fasta
    # modelo ENSG00000100280
    # len("ENSG00000100280") = 15
    hsIdLen = len(hsHeader3)
    hsAdd = 11 - hsIdLen
    aa = []
    for a1 in range(hsAdd):aa.append('0')
    hsAdd2 = "".join(aa)
    HsgeneID = "ENSG"+hsAdd2+hsHeader3
    HsProtID2,HsEnTrID2 = hsHeader2[1].split(" ")
    HsProtID4 = HsProtID2.replace(":","").replace("__","_") #hsHeader2[1].strip()
    HsProtID = HsgeneID+"-"+HsProtID4
    HsEnTrID3 = HsEnTrID2.replace("[","").split("]")[:-1]
    HsEnTrID = "_".join(HsEnTrID3)
    if HsgeneID in HSspProtDict.keys():
        HSspProtDict[HsgeneID]["HsProtID"].append(HsProtID)
        HSspProtDict[HsgeneID]["HsSeq"].append(HsSeq)
        HSspProtDict[HsgeneID]["HsensTr"].append(HsEnTrID)

    if HsgeneID not in HSspProtDict.keys():
        HSspProtDict[HsgeneID] = {"HsProtID":[HsProtID],"HsSeq":[HsSeq],"HsensTr":[HsEnTrID]}

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
    hsID, rnID = ortho.split("-")
    if hsID in ensDict.keys():
        info["hsensT"] = ensDict[hsID]["TrID"]
        info["hsensP"] = ensDict[hsID]["ProtID"]
    if hsID in refSeqDict.keys():
        info["hsrefseq"] = refSeqDict[hsID]["refseqTrID"]
#        info["hsensT"] = geneDict[hsID]["ensTrID"]
 #       info["hsensP"] = geneDict[hsID]["ensPrID"]
    """ if hsID in ccdsGeneDict.keys():
        info["hsccds"] = ccdsGeneDict[hsID] """
    if hsID in HSspProtDict.keys():
        info["hsSp"] = HSspProtDict[hsID]["HsProtID"]
        info["hsSpT"] = HSspProtDict[hsID]["HsensTr"]
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
    human = ""
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
                human = hsap
                hsPrint.append((header1,sequence1))
        for ccds2 in rnCCDSlist:
            if ccds2 in ccdsFastaD.keys():
                sequence2 = ccdsFastaD[ccds2]
                header2 = "|".join(["Rnov",rnov,ccds2])
                header12 = ">"+header2+"\n"
                rat = rnov
                rnPrint.append((header12,sequence2))
        if human != "" and rat != "":
            fileOutName = "ccds/needle/HsapRnov_OrthoProt/"+human+"_"+rat+"_ccds_ortho_proteins.fa"
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
#os.mkdir("refseq/needle/HsapRnov_OrthoProt")
#os.mkdir("refseq/needle/HsapRnov_FASTA")
hsFileOut2 = open("/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/final_pipeline/refseq/needle/HsapRnov_FASTA/hsap.faa","w")
rnFileOut2 = open("/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/final_pipeline/refseq/needle/HsapRnov_FASTA/rnov.faa","w")

refseqFASTAfileNames = ["/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/data/refseq/Hsap_refseq_prot.faa","/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/data/refseq/Rnov_refseq_prot.faa"]

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
    hsap,rnov = ortho2.split("-")
    human = ""
    rat = ""
    hsPrint = []
    rnPrint = []
    if hsap in refSeqDict.keys() and rnov in refSeqDict.keys():
        hsTrList = refSeqDict[hsap]["refseqTrID"]
        hsPrList = refSeqDict[hsap]["refseqPrID"]
        rnTrList = refSeqDict[rnov]["refseqTrID"]
        rnPrList = refSeqDict[rnov]["refseqPrID"]

    if infoID["hsrefseq"] != [] and infoID["rnrefseq"] != []:
        hsREFSEQlist = set(infoID["hsrefseq"])
        rnREFSEQlist = set(infoID["rnrefseq"])

        for refseq1 in hsREFSEQlist:
            if refseq1 in hsTrList:
                rstrIndex = hsTrList.index(refseq1)
                refseqPr1 = hsPrList[rstrIndex]
                if refseqPr1 in refseqFastaD.keys():
                    sequence1 = refseqFastaD[refseqPr1]
                    header = "|".join(["Hsap",hsap,refseqPr1,refseq1])
                    header1 = ">"+header+"\n"
                    human = hsap
                    hsPrint.append((header1,sequence1))
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
        if human != "" and rat != "":
            fileOutName = "/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/final_pipeline/refseq/needle/HsapRnov_OrthoProt/"+human+"_"+rat+"_refseq_ortho_proteins.fa"
            fileOut = open(fileOutName,"w")
            for i in range(len(hsPrint)):
                HeaderHs = hsPrint[i][0]
                SequenceHs = hsPrint[i][1]
                fileOut.write(HeaderHs)
                fileOut.write(SequenceHs)
                fileOut.write("\n")
                hsFileOut2.write(HeaderHs)
                hsFileOut2.write(SequenceHs)
                hsFileOut2.write("\n")

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
hsFileOut2.close()
rnFileOut2.close() 



#os.mkdir("spliceprot/")
#os.mkdir("spliceprot/needle/")
#os.mkdir("/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/final_pipeline/spliceprot/needle/HsapRnov_OrthoProt")
#os.mkdir("/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/final_pipeline/spliceprot/needle/HsapRnov_FASTA")
hsFileOut3 = open("/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/final_pipeline/spliceprot/needle/HsapRnov_FASTA/hsap.faa","w")
rnFileOut3 = open("/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/final_pipeline/spliceprot/needle/HsapRnov_FASTA/rnov.faa","w")

SpFASTAFileNames = ["/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/data/spliceprot_r100_v2/Hs.traducao_CORRIGIDA_release100.fasta","/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/data/spliceprot_r100_v2/Rn.traducao_CORRIGIDA_release100.fasta"]
SpFastaD = {}
for SpFASTAFilename in SpFASTAFileNames:
    SpFastaFile = open(SpFASTAFilename).read().split(">")[1:]
    for spFASTA in SpFastaFile:
        Header1,Seq = spFASTA.split("\n",1)
        Seq = Seq.strip()
        Header2 = Header1.split("_",1)
        if len(Seq) >= 20:
            if "Hs." in Header2[0]:
                Header3 = Header2[0].replace("Hs.","")
                #correcao do ID do fasta
                # modelo ENSG00000100280
                # len("ENSG00000100280") = 15
                IdLen = len(Header3)
                Add = 11 - IdLen
                aa = []
                for a1 in range(Add):
                    aa.append('0')
                Add2 = "".join(aa)
                geneID = "ENSG"+Add2+Header3
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
    hsap,rnov = ortho2.split("-")
    human = ""
    rat = ""
    hsPrint = []
    rnPrint = []
    if infoID["hsSp"] != [] and infoID["rnSp"] != []:
        hsSPlist = set(infoID["hsSp"])
        rnSPlist = set(infoID["rnSp"])
        hsSpliceProtL = infoID["hsSpT"]
        rnSpliceProtL = infoID["rnSpT"]

        for Sp1 in hsSPlist:
            sp1index = infoID["hsSp"].index(Sp1)
            sp1T = hsSpliceProtL[sp1index]

            if Sp1 in SpFastaD.keys():
                sequence1 = SpFastaD[Sp1]
                header = "|".join(["Hsap",Sp1,sp1T])
                header1 = ">"+header+"\n"
                human = hsap
                hsPrint.append((header1,sequence1))
        for Sp2 in rnSPlist:
            sp2index = infoID["rnSp"].index(Sp2)
            sp2T = rnSpliceProtL[sp2index]

            if Sp2 in SpFastaD.keys():
                sequence2 = SpFastaD[Sp2]
                header2 = "|".join(["Rnov",Sp2,sp2T])
                header12 = ">"+header2+"\n"
                rat = rnov
                rnPrint.append((header12,sequence2))
        if human != "" and rat != "":
            fileOutName = "/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/final_pipeline/spliceprot/needle/HsapRnov_OrthoProt/"+human+"_"+rat+"_ortho_proteins.fa"
            fileOut = open(fileOutName,"w")
            for i in range(len(hsPrint)):
                HeaderHs = hsPrint[i][0]
                SequenceHs = hsPrint[i][1]
                fileOut.write(HeaderHs)
                fileOut.write(SequenceHs)
                fileOut.write("\n")
                hsFileOut3.write(HeaderHs)
                hsFileOut3.write(SequenceHs)
                hsFileOut3.write("\n")

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
hsFileOut3.close()
rnFileOut3.close()

#os.mkdir("ensembl/")
#os.mkdir("ensembl/needle/")
#os.mkdir("ensembl/needle/HsapRnov_OrthoProt")
#os.mkdir("ensembl/needle/HsapRnov_FASTA")
hsFileOut4 = open("/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/final_pipeline/ensembl/needle/HsapRnov_FASTA/hsap.faa","w")
rnFileOut4 = open("/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/final_pipeline/ensembl/needle/HsapRnov_FASTA/rnov.faa","w")

for ortho2,infoID in OrthoDict.items():
    hsap,rnov = ortho2.split("-")
    human = ""
    rat = ""
    hsPrint = []
    rnPrint = []
    if infoID["hsensP"] != [] and infoID["rnensP"] != []:
        hsEnslist = set(infoID["hsensP"])
        rnEnslist = set(infoID["rnensP"])
        for Ens1 in hsEnslist:
            if Ens1 in EnsProtDict.keys():
                sequence1 = EnsProtDict[Ens1]["seq"]
                TrID = EnsProtDict[Ens1]["trID"]
                header = "|".join(["Hsap",hsap,TrID,Ens1])
                header1 = ">"+header+"\n"
                human = hsap
                hsPrint.append((header1,sequence1))
        for Ens2 in rnEnslist:
            if Ens2 in EnsProtDict.keys():
                sequence2 = EnsProtDict[Ens2]["seq"]
                TrID2 = EnsProtDict[Ens2]["trID"]
                header2 = "|".join(["Rnov",rnov,TrID2,Ens2])
                header12 = ">"+header2+"\n"
                rat = rnov
                rnPrint.append((header12,sequence2))
        if human != "" and rat != "":
            fileOutName = "/home/leticia.costa/disk1/leticia.costa/ortologos_flavia/final_pipeline/ensembl/needle/HsapRnov_OrthoProt/"+human+"_"+rat+"_ortho_proteins.fa"
            fileOut = open(fileOutName,"w")
            for i in range(len(hsPrint)):
                HeaderHs = hsPrint[i][0]
                SequenceHs = hsPrint[i][1]
                fileOut.write(HeaderHs)
                fileOut.write(SequenceHs)
                fileOut.write("\n")
                hsFileOut4.write(HeaderHs)
                hsFileOut4.write(SequenceHs)
                hsFileOut4.write("\n")
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
hsFileOut4.close()
rnFileOut4.close()

