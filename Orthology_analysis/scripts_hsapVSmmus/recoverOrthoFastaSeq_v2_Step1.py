import sys
import os
import itertools
import pickle


TableOrtho = open("/home/flavia.freitas/disk1/flavia.freitas/data/ortho_ensembl_101/Table_ID_Hsap_Mmus_Ortho.tsv").read().split("\n")[1:-1]
orthodict2 = {}
OrthoDict={}
for line in TableOrtho:
    col = line.split()
    #  considerar apenas ortologo 1:1
    hsapID2 = col[0]
    mmusID2 = col[1]
#    pair = hsapID+"-"+mmusID
    if "," not in mmusID2:
        if mmusID2 in orthodict2.keys():
            orthodict2[mmusID2].append(hsapID2)
        if mmusID2 not in orthodict2.keys():
            orthodict2[mmusID2] = [hsapID2]

OrthoDictPrint = {}
for mmusID,hsapID3 in orthodict2.items():
    if len(hsapID3) == 1:
        hsapID = "".join(hsapID3)
        pair = hsapID+"-"+mmusID
        OrthoDict[pair] = {"hsensT":[], "mmensT":[],
                            "hsensP":[], "mmensP":[],
                            "hsccds":[], "mmccds":[],
                            "hsSp":[], "mmSp":[],
                            "hsSpT":[], "mmSpT":[],
                            "hsrefseq":[],"mmrefseq":[]}
        OrthoDictPrint[hsapID] = mmusID 

pickleFileName = "/home/flavia.freitas/disk1/flavia.freitas/data/ortho_ensembl_101/Table_ID_Hsap_Mmus_Ortho_1-1.pickle"
pickleFileout =  open(pickleFileName, "wb")
pickle.dump(OrthoDictPrint,pickleFileout)
pickleFileout.close()

#hsEnsemblFileName = "/home/flavia.freitas/disk1/flavia.freitas/data/ensembl/Hsap_Gene_Proteins_Seq.tsv"
#mmEnsemblFileName = "/home/flavia.freitas/disk1/flavia.freitas/data/ensembl/Mmus_Gene_Proteins_Seq.tsv"

hsEnsemblFileName = "/home/flavia.freitas/disk1/flavia.freitas/data/ensembl/hsap/Homo_sapiens.GRCh38.pep.all.fa"
mmEnsemblFileName = "/home/flavia.freitas/disk1/flavia.freitas/data/ensembl/mmus/Mus_musculus.GRCm38.pep.all.fa"

ensemblFiles = [hsEnsemblFileName,mmEnsemblFileName]
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
        

hsGene2Ens = "/home/flavia.freitas/disk1/flavia.freitas/data/gene_info/hsap_gene2ensembl.txt"
mmGene2Ens = "/home/flavia.freitas/disk1/flavia.freitas/data/gene_info/mmus_gene2ensembl.txt"
hsGene2RefSeq = "/home/flavia.freitas/disk1/flavia.freitas/data/gene_info/hsap_gene2refseq.txt"
mmGene2RefSeq = "/home/flavia.freitas/disk1/flavia.freitas/data/gene_info/mmus_gene2refseq.txt"
Gene2EnsFiles = [hsGene2Ens,mmGene2Ens]
Gene2RefSeqFiles = [hsGene2RefSeq,mmGene2RefSeq]
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
            
hsCCDSFileName = "/home/flavia.freitas/disk1/flavia.freitas/data/ccds/Hs.CCDS2Sequence.current.txt"
mmCCDSFileName = "/home/flavia.freitas/disk1/flavia.freitas/data/ccds/Mm.CCDS2Sequence.current.txt"
CCDSFileNames = [hsCCDSFileName,mmCCDSFileName]
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


hsSpFASTA = open("/home/flavia.freitas/disk1/flavia.freitas/data/spliceprot_r100_v2/Hs.traducao_CORRIGIDA_release100.fasta").read().split(">")[1:]
mmSpFASTA = open("/home/flavia.freitas/disk1/flavia.freitas/data/spliceprot_r100_v2/Mm.traducao_CORRIGIDA_release100.fasta").read().split(">")[1:]
HSspProtDict = {}
MMspProtDict = {}

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

for mmfasta in mmSpFASTA:
    mmHeader1,mmSeq = mmfasta.split("\n",1)
    MmSeq = mmSeq.strip()
    mmHeader2 = mmHeader1.split("_",1)
    mmHeader3 = mmHeader2[0].replace("Mm.","")
    #correcao do ID do fasta
    # modelo ENSMUSG00000093804 
    # len("ENSMUSG00000093804") = 18
    mmIdLen = len(mmHeader3)
    mmAdd = 11 - mmIdLen
    bb = []
    for b1 in range(mmAdd):
        bb.append('0')
    mmAdd2 = "".join(bb)
    MmgeneID = "ENSMUSG"+mmAdd2+mmHeader3
    MmProtID2,MmEnTrID2 = mmHeader2[1].split(" ")
    MmProtID4 = MmProtID2.replace(":","").replace("__","_") #mmHeader2[1].strip()
    MmProtID = MmgeneID+"-"+MmProtID4
    MmEnTrID3 = MmEnTrID2.replace("[","").split("]")[:-1]
    MmEnTrID = "_".join(MmEnTrID3)
    if MmgeneID in MMspProtDict.keys():
        MMspProtDict[MmgeneID]["MmProtID"].append(MmProtID)
        MMspProtDict[MmgeneID]["MmSeq"].append(MmSeq)
        MMspProtDict[MmgeneID]["MmensTr"].append(MmEnTrID)
    if MmgeneID not in MMspProtDict.keys():
        MMspProtDict[MmgeneID] = {"MmProtID": [MmProtID], 
                                  "MmSeq": [MmSeq],
                                  "MmensTr": [MmEnTrID]}


for ortho,info in OrthoDict.items():
    hsID, mmID = ortho.split("-")
    if hsID in ensDict.keys():
        info["hsensT"] = ensDict[hsID]["TrID"]
        info["hsensP"] = ensDict[hsID]["ProtID"]
    if hsID in refSeqDict.keys():
        info["hsrefseq"] = refSeqDict[hsID]["refseqTrID"]
#        info["hsensT"] = geneDict[hsID]["ensTrID"]
 #       info["hsensP"] = geneDict[hsID]["ensPrID"]
    if hsID in ccdsGeneDict.keys():
        info["hsccds"] = ccdsGeneDict[hsID]
    if hsID in HSspProtDict.keys():
        info["hsSp"] = HSspProtDict[hsID]["HsProtID"]
        info["hsSpT"] = HSspProtDict[hsID]["HsensTr"]
    if mmID in ensDict.keys():
        info["mmensT"] = ensDict[mmID]["TrID"]
        info["mmensP"] = ensDict[mmID]["ProtID"]
    if mmID in refSeqDict.keys():
        info["mmrefseq"] = refSeqDict[mmID]["refseqTrID"]
    #    info["mmensT"] = geneDict[mmID]["ensTrID"]
     #   info["mmensP"] = geneDict[mmID]["ensPrID"]
    if mmID in ccdsGeneDict.keys():
        info["mmccds"] = ccdsGeneDict[mmID]
    if mmID in MMspProtDict.keys():
        info["mmSp"] = MMspProtDict[mmID]["MmProtID"]
        info["mmSpT"] = MMspProtDict[mmID]["MmensTr"]


#for a,b in OrthoDict.items():
 #   print(a,b)

hsCCDSFASTAfileName = "/home/flavia.freitas/disk1/flavia.freitas/data/ccds/Hsap_CCDS_protein.current.faa"
mmCCDSFASTAfileName = "/home/flavia.freitas/disk1/flavia.freitas/data/ccds/Mmus_CCDS_protein.current.faa"
ccdsFASTAfileNames = [hsCCDSFASTAfileName,mmCCDSFASTAfileName]

#hsCCDSFasta = open(hsCCDSFASTAfileName).read().split(">")[1:]
#mmCCDSFasta = open(mmCCDSFASTAfileName).read().split(">")[1:]
ccdsFastaD = {}
for ccdsFastaFileName in ccdsFASTAfileNames:
    CCDSFastaFile = open(ccdsFastaFileName).read().split(">")[1:]   
    for ccdsFasta in CCDSFastaFile:
        h,s = ccdsFasta.split("\n",1)
        seq = s.replace("\n","")
        h2 = h.split("|")
        h3,lixo = h2[0].split(".")
        if len(seq) >= 20:
            ccdsFastaD[h3] = seq

os.mkdir("ccds/")
os.mkdir("ccds/needle/") 
os.mkdir("ccds/needle/HsapMmus_OrthoProt")
os.mkdir("ccds/needle/HsapMmus_FASTA")
hsFileOut = open("ccds/needle/HsapMmus_FASTA/hsap.faa","w")
mmFileOut = open("ccds/needle/HsapMmus_FASTA/mmus.faa","w")

for ortho2,infoID in OrthoDict.items():
    hsap,mmus = ortho2.split("-")
    human = ""
    mouse = ""
    hsPrint = []
    mmPrint = []
    if infoID["hsccds"] != [] and infoID["mmccds"] != []:
        hsCCDSlist = set(infoID["hsccds"])
        mmCCDSlist = set(infoID["mmccds"])
        for ccds1 in hsCCDSlist:
            if ccds1 in ccdsFastaD.keys():
                sequence1 = ccdsFastaD[ccds1]
                header = "|".join(["Hsap",hsap,ccds1])
                header1 = ">"+header+"\n"
                human = hsap
                hsPrint.append((header1,sequence1))
        for ccds2 in mmCCDSlist:
            if ccds2 in ccdsFastaD.keys():
                sequence2 = ccdsFastaD[ccds2]
                header2 = "|".join(["Mmus",mmus,ccds2])
                header12 = ">"+header2+"\n"
                mouse = mmus
                mmPrint.append((header12,sequence2))
        if human != "" and mouse != "":
            fileOutName = "ccds/needle/HsapMmus_OrthoProt/"+human+"_"+mouse+"_ccds_ortho_proteins.fa"
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

            for j in range(len(mmPrint)):
                 HeaderMm = mmPrint[j][0]
                 SequenceMm = mmPrint[j][1]
                 fileOut.write(HeaderMm)
                 fileOut.write(SequenceMm)
                 fileOut.write("\n")
                 mmFileOut.write(HeaderMm)
                 mmFileOut.write(SequenceMm)
                 mmFileOut.write("\n")
        fileOut.close()
hsFileOut.close()
mmFileOut.close()


os.mkdir("refseq/")
os.mkdir("refseq/needle/")
os.mkdir("refseq/needle/HsapMmus_OrthoProt")
os.mkdir("refseq/needle/HsapMmus_FASTA")
hsFileOut2 = open("refseq/needle/HsapMmus_FASTA/hsap.faa","w")
mmFileOut2 = open("refseq/needle/HsapMmus_FASTA/mmus.faa","w")

refseqFASTAfileNames = ["/home/flavia.freitas/disk1/flavia.freitas/data/refseq/Hsap_refseq_prot.faa","/home/flavia.freitas/disk1/flavia.freitas/data/refseq/Mmus_refseq_prot.faa"]

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
    hsap,mmus = ortho2.split("-")
    human = ""
    mouse = ""
    hsPrint = []
    mmPrint = []
    if hsap in refSeqDict.keys() and mmus in refSeqDict.keys():
        hsTrList = refSeqDict[hsap]["refseqTrID"]
        hsPrList = refSeqDict[hsap]["refseqPrID"]
        mmTrList = refSeqDict[mmus]["refseqTrID"]
        mmPrList = refSeqDict[mmus]["refseqPrID"]

    if infoID["hsrefseq"] != [] and infoID["mmrefseq"] != []:
        hsREFSEQlist = set(infoID["hsrefseq"])
        mmREFSEQlist = set(infoID["mmrefseq"])
        
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
        for refseq2 in mmREFSEQlist:
            if refseq2 in mmTrList:
                rstrIndex2 = mmTrList.index(refseq2)
                refseqPr2 = mmPrList[rstrIndex2]
                if refseqPr2 in refseqFastaD.keys():
                    sequence2 = refseqFastaD[refseqPr2]
                    header2 = "|".join(["Mmus",mmus,refseqPr2,refseq2])
                    header12 = ">"+header2+"\n"
                    mouse = mmus
                    mmPrint.append((header12,sequence2))
        if human != "" and mouse != "":
            fileOutName = "refseq/needle/HsapMmus_OrthoProt/"+human+"_"+mouse+"_refseq_ortho_proteins.fa"
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

            for j in range(len(mmPrint)):
                HeaderMm = mmPrint[j][0]
                SequenceMm = mmPrint[j][1]
                fileOut.write(HeaderMm)
                fileOut.write(SequenceMm)
                fileOut.write("\n")
                mmFileOut2.write(HeaderMm)
                mmFileOut2.write(SequenceMm)
                mmFileOut2.write("\n")
        fileOut.close()
hsFileOut2.close()
mmFileOut2.close()



os.mkdir("spliceprot/")
os.mkdir("spliceprot/needle/")
os.mkdir("spliceprot/needle/HsapMmus_OrthoProt")
os.mkdir("spliceprot/needle/HsapMmus_FASTA")
hsFileOut3 = open("spliceprot/needle/HsapMmus_FASTA/hsap.faa","w")
mmFileOut3 = open("spliceprot/needle/HsapMmus_FASTA/mmus.faa","w")

SpFASTAFileNames = ["/home/flavia.freitas/disk1/flavia.freitas/data/spliceprot_r100_v2/Hs.traducao_CORRIGIDA_release100.fasta","/home/flavia.freitas/disk1/flavia.freitas/data/spliceprot_r100_v2/Mm.traducao_CORRIGIDA_release100.fasta"]
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
            if "Mm." in Header2[0]:
                Header3 = Header2[0].replace("Mm.","")
                #correcao do ID do fasta
                # modelo ENSMUSG00000093804
                # len("ENSMUSG00000093804") = 18
                IdLen = len(Header3)
                Add = 11 - IdLen
                aa = []
                for a1 in range(Add):
                    aa.append('0')
                Add2 = "".join(aa)
                geneID = "ENSMUSG"+Add2+Header3

            ProtID2,EnTrID2 = Header2[1].split(" ")
            ProtID3 = ProtID2.replace(":","").replace("__","_") 
            ProtID = geneID+"-"+ProtID3 
            EnTrID3 = EnTrID2.replace("[","").split("]")[:-1]
            EnTrID = "_".join(EnTrID3)
            SpFastaD[ProtID] = Seq 

for ortho2,infoID in OrthoDict.items():
    hsap,mmus = ortho2.split("-")
    human = ""
    mouse = ""
    hsPrint = []
    mmPrint = []
    if infoID["hsSp"] != [] and infoID["mmSp"] != []:
        hsSPlist = set(infoID["hsSp"])
        mmSPlist = set(infoID["mmSp"])
        hsSpliceProtL = infoID["hsSpT"]
        mmSpliceProtL = infoID["mmSpT"]

        for Sp1 in hsSPlist:
            sp1index = infoID["hsSp"].index(Sp1)
            sp1T = hsSpliceProtL[sp1index]

            if Sp1 in SpFastaD.keys():
                sequence1 = SpFastaD[Sp1]
                header = "|".join(["Hsap",Sp1,sp1T])
                header1 = ">"+header+"\n"
                human = hsap
                hsPrint.append((header1,sequence1))
        for Sp2 in mmSPlist:
            sp2index = infoID["mmSp"].index(Sp2)
            sp2T = mmSpliceProtL[sp2index]

            if Sp2 in SpFastaD.keys():
                sequence2 = SpFastaD[Sp2]
                header2 = "|".join(["Mmus",Sp2,sp2T])
                header12 = ">"+header2+"\n"
                mouse = mmus
                mmPrint.append((header12,sequence2))
        if human != "" and mouse != "":
            fileOutName = "spliceprot/needle/HsapMmus_OrthoProt/"+human+"_"+mouse+"_ortho_proteins.fa"
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

            for j in range(len(mmPrint)):
                HeaderMm = mmPrint[j][0]
                SequenceMm = mmPrint[j][1]
                fileOut.write(HeaderMm)
                fileOut.write(SequenceMm)
                fileOut.write("\n")
                mmFileOut3.write(HeaderMm)
                mmFileOut3.write(SequenceMm)
                mmFileOut3.write("\n")
        fileOut.close()
hsFileOut3.close()
mmFileOut3.close()

os.mkdir("ensembl/")
os.mkdir("ensembl/needle/")
os.mkdir("ensembl/needle/HsapMmus_OrthoProt")
os.mkdir("ensembl/needle/HsapMmus_FASTA")
hsFileOut4 = open("ensembl/needle/HsapMmus_FASTA/hsap.faa","w")
mmFileOut4 = open("ensembl/needle/HsapMmus_FASTA/mmus.faa","w")

for ortho2,infoID in OrthoDict.items():
    hsap,mmus = ortho2.split("-")
    human = ""
    mouse = ""
    hsPrint = []
    mmPrint = []
    if infoID["hsensP"] != [] and infoID["mmensP"] != []:
        hsEnslist = set(infoID["hsensP"])
        mmEnslist = set(infoID["mmensP"])
        for Ens1 in hsEnslist:
            if Ens1 in EnsProtDict.keys():
                sequence1 = EnsProtDict[Ens1]["seq"]
                TrID = EnsProtDict[Ens1]["trID"]
                header = "|".join(["Hsap",hsap,TrID,Ens1])
                header1 = ">"+header+"\n"
                human = hsap
                hsPrint.append((header1,sequence1))
        for Ens2 in mmEnslist:
            if Ens2 in EnsProtDict.keys():
                sequence2 = EnsProtDict[Ens2]["seq"]
                TrID2 = EnsProtDict[Ens2]["trID"]
                header2 = "|".join(["Mmus",mmus,TrID2,Ens2])
                header12 = ">"+header2+"\n"
                mouse = mmus
                mmPrint.append((header12,sequence2))
        if human != "" and mouse != "":
            fileOutName = "ensembl/needle/HsapMmus_OrthoProt/"+human+"_"+mouse+"_ortho_proteins.fa"
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
            for j in range(len(mmPrint)):
                HeaderMm = mmPrint[j][0]
                SequenceMm = mmPrint[j][1]
                fileOut.write(HeaderMm)
                fileOut.write(SequenceMm)
                fileOut.write("\n")
                mmFileOut4.write(HeaderMm)
                mmFileOut4.write(SequenceMm)
                mmFileOut4.write("\n")


        fileOut.close()
hsFileOut4.close()
mmFileOut4.close()

