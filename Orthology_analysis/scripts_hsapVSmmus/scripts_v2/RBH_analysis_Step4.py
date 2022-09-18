import sys
import os
import subprocess
import pickle
import pprint
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline
import itertools


# gerar lista e dicionário com pares
# dos genes ortologos de hsap e mmus
# usar essa base para filtrar os pares encontrados pelo RBH 
# antes de proceder com a analise dos resultados
# caso haja algum par de isoformas que não seja de genes ortologos

pickleOrthoFile = open("/home/flavia.freitas/disk1/flavia.freitas/data/ortho_ensembl_101/Table_ID_Hsap_Mmus_Ortho_1-1.pickle","rb")
OrthoDict = pickle.load(pickleOrthoFile)
orthoList = []

for hsapID,mmusID in OrthoDict.items():
    pair1 = hsapID+" "+mmusID
    orthoList.append(pair1)

RBHFilePath = sys.argv[1]
RBHFilePathInfo = RBHFilePath.split("/")
RBHFileName2 = RBHFilePath
RBHFolderPath = "/".join(RBHFilePathInfo[:-3])
RBHFolderPath = RBHFolderPath + "/"
# database options
# ccds, refseq, spliceprot, ensembl
database = RBHFilePathInfo[0]


os.mkdir(RBHFolderPath+"analysis_RBH")
def RBH_Filter(RBHFileName,sourceDB):
    
    RBHFile = open(RBHFileName).read().split("\n")[:-1]
    # database options
    # ccds, refseq, spliceprot, ensembl
    
    RBH_header_itens = RBHFile.pop(0).replace("#","").split()
    RBH_header = "\t".join(RBH_header_itens)
    RBH_NonMatchGene = open(RBHFolderPath+"analysis_RBH/RBH_nonMatch_OrthoGenes.tab.txt","w")
    RBHlist = []
    RBHdict = {}
    RBH_nonMatchList = []
    for linha in RBHFile:
        coluna = linha.split("\t")
        if "spliceprot" in sourceDB: 
            if "|" in coluna[0]:
                hsT1,hsEnsT = coluna[0].split("|")
            if "|" not in coluna[0]:
                hsT1 = coluna[0].replace("Hsap_","")
            if "|" in coluna[1]:
                mmT1,mmEnsT = coluna[1].split("|")
            if "|" not in coluna[1]:
                mmT1 = coluna[1].replace("Mmus_","")
            Tpair = hsT1+" "+mmT1
            RBHFileInfo = coluna[2:]
            hsGene,hsvar = hsT1.split("-",1)
            mmGene,mmvar = mmT1.split("-",1)
            hsGene = hsGene.strip()
            mmGene = mmGene.strip()
            pair2 = hsGene+" "+mmGene
        if "refseq" in sourceDB:
            infoHs = coluna[0].replace("|","-").replace("Hsap_","").split("-")
            infoMm = coluna[1].replace("|","-").replace("Mmus_","").split("-")
            hsGene = infoHs[0].strip()
            hsT1 = infoHs[2]
            hsProt1 = infoHs[1]
            mmGene = infoMm[0].strip()
            mmT1 = infoMm[2]
            mmProt1 = infoMm[1]
            Tpair = hsProt1+" "+mmProt1 #hsT1+" "+mmT1
            pair2 = hsGene+" "+mmGene
            RBHFileInfo = coluna[2:]
        if "ccds" in sourceDB:
            hsInfo = coluna[0].replace("|","-").replace("Hsap_","").split("-")
            mmInfo = coluna[1].replace("|","-").replace("Mmus_","").split("-")
            hsccdsID = hsInfo[1]#.split(".")
            mmccdsID = mmInfo[1]#.split(".")
            hsGene = hsInfo[0].strip()
            mmGene = mmInfo[0].strip()
            #hsT1 = "_".join(hsInfo[2:])
            #mmT1 = "_".join(mmInfo[2:])
            pair2 = hsGene+" "+mmGene
            Tpair = hsccdsID+" "+mmccdsID
            RBHFileInfo = coluna[2:]
        if "ensembl" in sourceDB:
            hsInfo = coluna[0].replace("|","-").replace("Hsap_","").split("-")
            mmInfo = coluna[1].replace("|","-").replace("Mmus_","").split("-")
            hsT1 = hsInfo[1]#.split(".")
            mmT1 = mmInfo[1]#.split(".")
            hsProt1 = hsInfo[2]
            mmProt1 = mmInfo[2]
            Tpair = hsProt1+" "+mmProt1 #hsT1+" "+mmT1
            hsGene = hsInfo[0].strip()
            mmGene = mmInfo[0].strip()
            pair2 = hsGene+" "+mmGene
            RBHFileInfo = coluna[2:]
        
        if pair2 not in orthoList:
            if pair2 not in RBH_nonMatchList:
                RBH_nonMatchList.append(pair2)
        if pair2 in orthoList:
            RBHdict[Tpair] = {"genePair":pair2, "RBHinfo":RBHFileInfo}

    for p2 in RBH_nonMatchList:
        RBH_NonMatchGene.write(p2)
        RBH_NonMatchGene.write("\n")

    RBH_NonMatchGene.close()
    return(RBHdict,RBH_header)


RBHFunctionResult = RBH_Filter(RBHFileName2,database)
RBHdict2 = RBHFunctionResult[0]
RBH_Header = RBHFunctionResult[1]
RBHHeader = "HsGeneID"+"\t"+"MmGeneID"+"\t"+RBH_Header
RBHdict3 ={}

#pprint.pprint(RBHdict2)
for key,value in RBHdict2.items():
    v = value["RBHinfo"]
    genePair = value["genePair"]
    hs1,mm1 = key.split(" ")
    evalue = float(v[0])
    bitscore = float(v[1])
    qStart = float(v[2])
    qEnd = float(v[3])
    qCover = float(v[4])
    tStart = float(v[5])
    tEnd = float(v[6])
    tCover = float(v[7])
    qualifier = v[8].strip()

    if qualifier == "RBH": 
        if qCover >= 60 and tCover >= 60:
            if hs1 in RBHdict3.keys():
                RBHdict3[hs1]["mm1"].append(mm1)
                RBHdict3[hs1]["evalue"].append(evalue)
                RBHdict3[hs1]["bitscore"].append(bitscore)
                RBHdict3[hs1]["qStart"].append(qStart)
                RBHdict3[hs1]["qEnd"].append(qEnd)
                RBHdict3[hs1]["qCover"].append(qCover)
                RBHdict3[hs1]["tStart"].append(tStart)
                RBHdict3[hs1]["tEnd"].append(tEnd)
                RBHdict3[hs1]["tCover"].append(tCover)
                RBHdict3[hs1]["qualifier"].append(qualifier)
                RBHdict3[hs1]["genePair"].append(genePair)

            if hs1 not in RBHdict3.keys():
                RBHdict3[hs1] = {"mm1":[mm1],
                        "evalue":[evalue],
                        "bitscore":[bitscore],
                        "qStart":[qStart],
                        "qEnd":[qEnd],
                        "qCover":[qCover],
                        "tStart":[tStart],
                        "tEnd":[tEnd],
                        "tCover":[tCover],
                        "qualifier":[qualifier],
                        "genePair":[genePair]}

#pprint.pprint(RBHdict3)
sortDict = {}
for hsT,Info in RBHdict3.items():
    mmTarget = Info["mm1"]
    qCoverList = Info["qCover"]
    tCoverList = Info["tCover"]
    evalueList = Info["evalue"]
    bitScoreList = Info["bitscore"]
    GenePairList = Info["genePair"]
    qualifierList = Info["qualifier"]
    tEndList = Info["tEnd"]
    tStartList = Info["tStart"]
    qStartList = Info["qStart"]
    qEndList = Info["qEnd"]
    indexList = []
    
    for qCIndex,qCvalue in enumerate(qCoverList):
        if qCoverList.count(max(qCoverList)) > 1:
            if float(qCvalue) == float(max(qCoverList)):
                    indexList.append(qCIndex)
        if qCoverList.count(max(qCoverList)) == 1:
            if float(qCvalue) == float(max(qCoverList)):# or \
                        #float(tCoverList[qCIndex]) == float(max(tCoverList)) or \
                        #float(evalueList[qCIndex]) == float(min(evalueList)):
                    indexList = [qCIndex]
    finalIndex = ""
    tCoverSubconj = []
    bitScoreSubconj = []
    evalueSubconj = []
    if len(indexList) == 1:
        finalIndex = indexList[0]
    if len(indexList) > 1:
        for qcindex in indexList:
            tCoverSubconj.append(tCoverList[qcindex])
            bitScoreSubconj.append(bitScoreList[qcindex])
            evalueSubconj.append(evalueList[qcindex])
        maxTCoverSubconj = max(tCoverSubconj)
        minEvalueSubconj = min(evalueSubconj)
        if tCoverSubconj.count(maxTCoverSubconj) == 1:
            finalIndex = tCoverList.index(maxTCoverSubconj)
        if tCoverSubconj.count(maxTCoverSubconj) > 1:
            finalIndex = evalueList.index(minEvalueSubconj)

    if type(finalIndex) != str:
        #mmT = mmTarget[finalIndex]
        if hsT in sortDict.keys():
            sortDict[hsT]["mm1"].append(mmTarget[finalIndex])
            sortDict[hsT]["qCover"].append(qCoverList[finalIndex])
            sortDict[hsT]["tCover"].append(tCoverList[finalIndex])
            sortDict[hsT]["evalue"].append(evalueList[finalIndex])
            sortDict[hsT]["bitscore"].append(bitScoreList[finalIndex])
            sortDict[hsT]["genePair"].append(GenePairList[finalIndex])
            sortDict[hsT]["qualifier"].append(qualifierList[finalIndex])
            sortDict[hsT]["tEnd"].append(tEndList[finalIndex])
            sortDict[hsT]["qEnd"].append(qEndList[finalIndex])
            sortDict[hsT]["qStart"].append(qStartList[finalIndex])
            sortDict[hsT]["tStart"].append(tStartList[finalIndex])
        if hsT not in sortDict.keys():
            sortDict[hsT] = {"mm1":[mmTarget[finalIndex]],
                                  "qCover": [qCoverList[finalIndex]],
                                  "tCover":[tCoverList[finalIndex]],
                                  "evalue":[evalueList[finalIndex]],
                                  "bitscore":[bitScoreList[finalIndex]],
                                  "genePair":[GenePairList[finalIndex]],
                                  "qualifier":[qualifierList[finalIndex]],
                                  "tEnd":[tEndList[finalIndex]],
                                  "qEnd":[qEndList[finalIndex]],
                                  "qStart":[qStartList[finalIndex]],
                                  "tStart":[tStartList[finalIndex]]}

#pprint.pprint(sortDict)
#pprint.pprint(sortDict)
sort2Dict = {}
mouseLista = []
for human,rbhInfo in sortDict.items():
    qCoverList2 = rbhInfo["qCover"]
    tCoverList2 = rbhInfo["tCover"]
    eValueList2 = rbhInfo["evalue"]
    bitScoreList2 = rbhInfo["bitscore"]
    mmList2 = rbhInfo["mm1"]
    GenePairList2 = rbhInfo["genePair"]
    qualifierList2 = rbhInfo["qualifier"]
    tEndList2 = rbhInfo["tEnd"]
    tStartList2 = rbhInfo["tStart"]
    qStartList2 = rbhInfo["qStart"]
    qEndList2 = rbhInfo["qEnd"]

    if len(qCoverList2) == 1:
        if qCoverList2[0] >= 60 and tCoverList2[0] >= 90 or\
                tCoverList2[0]>= 60 and qCoverList2[0] >= 90:
                mmID = mmList2[0]
                mouseLista.append(mmID)
                sort2Dict[human] = {"target":mmID,
                                "qCover": qCoverList2[0],
                                "tCover": tCoverList2[0],
                                "evalue": eValueList2[0],
                                "bitscore":bitScoreList2[0],
                                "genePair":GenePairList2[0],
                                "qualifier":qualifierList2[0],
                                "tEnd":tEndList2[0],
                                "qEnd":qEndList2[0],
                                "qStart":qStartList2[0],
                                "tStart":tStartList2[0]}


    if len(qCoverList2) > 1:
        if qCoverList2.count(max(qCoverList2)) == 1:
            FinalIndex = qCoverList2.index(max(qCoverList2))
            if qCoverList2[FinalIndex] >= 60 and tCoverList2[FinalIndex] >= 90 or\
                    tCoverList2[FinalIndex] >= 60 and qCoverList2[FinalIndex]  >= 90:
                mmID = mmList2[FinalIndex]
                mouseLista.append(mmID)
                sort2Dict[human] = {"target":mmID,
                                "qCover": qCoverList2[FinalIndex],
                                "tCover": tCoverList2[FinalIndex],
                                "evalue": eValueList2[FinalIndex],
                                "bitscore":bitScoreList2[FinalIndex],
                                "genePair":GenePairList2[FinalIndex],
                                "qualifier":qualifierList2[FinalIndex],
                                "tEnd":tEndList2[FinalIndex],
                                "qEnd":qEndList2[FinalIndex],
                                "qStart":qStartList2[FinalIndex],
                                "tStart":tStartList2[FinalIndex]}
        if qCoverList2.count(max(qCoverList2)) > 1:
            if tCoverList2.count(max(tCoverList2)) == 1:
                FinalIndex = tCoverList2.index(max(tCoverList2))
                if qCoverList2[FinalIndex] >= 60 and tCoverList2[FinalIndex] >= 90 or\
                        tCoverList2[FinalIndex] >= 60 and qCoverList2[FinalIndex] >= 90:
                    mmID = mmList2[FinalIndex]
                    mouseLista.append(mmID)
                    sort2Dict[human] = {"target":mmID,
                                    "qCover": qCoverList2[FinalIndex],
                                   "tCover": tCoverList2[FinalIndex],
                                   "evalue": eValueList2[FinalIndex],
                                   "bitscore":bitScoreList2[FinalIndex],
                                   "genePair":GenePairList2[FinalIndex],
                                   "qualifier":qualifierList2[FinalIndex],
                                   "tEnd":tEndList2[FinalIndex],
                                   "qEnd":qEndList2[FinalIndex],
                                   "qStart":qStartList2[FinalIndex],
                                   "tStart":tStartList2[FinalIndex]}

            if tCoverList2.count(max(tCoverList2)) > 1:
                FinalIndex = eValueList2.index(min(eValueList2))
                if qCoverList2[FinalIndex] >= 60 and tCoverList2[FinalIndex] >= 90 or \
                        tCoverList2[FinalIndex] >= 60 and qCoverList2[FinalIndex]  >= 90:

                    mmID = mmList2[FinalIndex]
                    mouseLista.append(mmID)
                    sort2Dict[human] = {"target":mmID,
                                    "qCover": qCoverList2[FinalIndex],
                                    "tCover": tCoverList2[FinalIndex],
                                    "evalue": eValueList2[FinalIndex],
                                    "bitscore":bitScoreList2[FinalIndex],
                                    "genePair":GenePairList2[FinalIndex],
                                    "qualifier":qualifierList2[FinalIndex],
                                    "tEnd":tEndList2[FinalIndex],
                                    "qEnd":qEndList2[FinalIndex],
                                    "qStart":qStartList2[FinalIndex],
                                    "tStart":tStartList2[FinalIndex]}

#pprint.pprint(sort2Dict)
#listaHs=[]
for hsTa,Infoa in RBHdict3.items():
    if hsTa not in sort2Dict.keys():
 #       listaHs.append(hsTa)
        mmTargeta = Infoa["mm1"]
        qCoverLista = Infoa["qCover"]
        tCoverLista = Infoa["tCover"]
        evalueLista = Infoa["evalue"]
        bitScoreLista = Infoa["bitscore"]
        GenePairLista = Infoa["genePair"]
        qualifierLista = Infoa["qualifier"]
        tEndLista = Infoa["tEnd"]
        tStartLista = Infoa["tStart"]
        qStartLista = Infoa["qStart"]
        qEndLista = Infoa["qEnd"]
        indexLista = []
        for tCIndex,tCvalue in enumerate(tCoverLista):
    #        if tCoverLista.count(max(tCoverLista)) > 1:
            if float(tCvalue) == float(max(tCoverLista)):
                #indexLista.append(tCIndex)
                finalIndexA = tCIndex
                #if tCoverLista.count(max(tCoverLista)) == 1:
                 #   if float(tCvalue) == float(max(tCoverLista)):
                  #      indexLista = [tCIndex]
#        finalIndexA = ""
#        qCoverSubconjA = []
 #       bitScoreSubconjA = []
  #      evalueSubconjA = []
   #     if len(indexLista) == 1:
    #        finalIndexA = indexLista[0]
     #   if len(indexLista) > 1:
      #      for tcindex in indexLista:
       #         qCoverSubconjA.append(qCoverLista[tcindex])
        #        bitScoreSubconjA.append(bitScoreLista[tcindex])
         #       evalueSubconjA.append(evalueLista[tcindex])
          #  maxQCoverSubconjA = max(qCoverSubconjA)
           # minEvalueSubconjA = min(evalueSubconjA)

            #if qCoverSubconjA.count(maxQCoverSubconjA) == 1:
             #   finalIndexA = qCoverLista.index(maxQCoverSubconjA)
            #if qCoverSubconjA.count(maxQCoverSubconjA) > 1:
             #   finalIndexA = evalueLista.index(minEvalueSubconjA)

        #if type(finalIndexA) != str:
        sort2Dict[hsTa] = {"target":mmTargeta[finalIndexA],
                "qCover": qCoverLista[finalIndexA],
                "tCover":tCoverLista[finalIndexA],
                "evalue":evalueLista[finalIndexA],
                "bitscore":bitScoreLista[finalIndexA],
                "genePair":GenePairLista[finalIndexA],
                "qualifier":qualifierLista[finalIndexA],
                "tEnd":tEndLista[finalIndexA],
                "qEnd":qEndLista[finalIndexA],
                "qStart":qStartLista[finalIndexA],
                "tStart":tStartLista[finalIndexA]}
#print("################")
#pprint.pprint(sort2Dict)
sort3Dict = {}
for human2,info2 in sort2Dict.items():
    qCover3 = info2["qCover"]
    tCover3 = info2["tCover"]
    eValue3 = info2["evalue"]
    bitScore3 = info2["bitscore"]
    mm3 = info2["target"]
    GenePair3 = info2["genePair"]
    qualifier3 = info2["qualifier"]
    tEnd3 = info2["tEnd"]
    tStart3 = info2["tStart"]
    qStart3 = info2["qStart"]
    qEnd3 = info2["qEnd"]
    
    if mm3 in sort3Dict.keys():
        qCover2 = qCover3
        qCover1 = sort3Dict[mm3]["qCover"]
        if qCover2 > qCover1:
            sort3Dict[mm3]["hsT"] = human2 
            sort3Dict[mm3]["qCover"] = qCover2 
            sort3Dict[mm3]["tCover"] = tCover3 
            sort3Dict[mm3]["evalue"] = eValue3 
            sort3Dict[mm3]["bitscore"] = bitScore3 
            sort3Dict[mm3]["genePair"] = GenePair3 
            sort3Dict[mm3]["qualifier"] = qualifier3 
            sort3Dict[mm3]["tEnd"] = tEnd3 
            sort3Dict[mm3]["qEnd"] = qEnd3 
            sort3Dict[mm3]["qStart"] = qStart3 
            sort3Dict[mm3]["tStart"] = tStart3
        if qCover2 == qCover1:
            tCover2 = tCover3
            tCover1 = sort3Dict[mm3]["tCover"]
            if tCover2 > tCover1:
                sort3Dict[mm3]["hsT"] = human2 
                sort3Dict[mm3]["qCover"] = qCover2 
                sort3Dict[mm3]["tCover"] = tCover2 
                sort3Dict[mm3]["evalue"] = eValue3 
                sort3Dict[mm3]["bitscore"] = bitScore3 
                sort3Dict[mm3]["genePair"] = GenePair3 
                sort3Dict[mm3]["qualifier"] = qualifier3 
                sort3Dict[mm3]["tEnd"] = tEnd3 
                sort3Dict[mm3]["qEnd"] = qEnd3 
                sort3Dict[mm3]["qStart"] = qStart3 
                sort3Dict[mm3]["tStart"] = tStart3


    if mm3 not in sort3Dict.keys():
        sort3Dict[mm3] = {"hsT": human2,
                "qCover": qCover3,
                "tCover": tCover3,
                "evalue": eValue3,
                "bitscore":bitScore3,
                "genePair":GenePair3,
                "qualifier":qualifier3,
                "tEnd":tEnd3,
                "qEnd":qEnd3,
                "qStart":qStart3,
                "tStart":tStart3}


#pprint.pprint(sort3Dict)
#print(">>")
#print("\n".join(listaHs))
RBHResultThresName = RBHFolderPath+"analysis_RBH/"+database+"_RBH_ortho_60-90_tab.txt"    
RBHResultThres = open(RBHResultThresName,"w")    
RBHResultThres.write(RBHHeader)
RBHResultThres.write("\n") 

for a2,b2 in sort3Dict.items():
    hsG,mmG = b2["genePair"].split(" ")
    mouseID = a2
    humanID = b2["hsT"]
    info3 = [hsG,mmG,humanID,mouseID,b2["evalue"], b2["bitscore"],b2["qStart"],b2["qEnd"],\
                b2["qCover"], b2["tStart"],b2["tEnd"],b2["tCover"],b2["qualifier"]]      
    info4 ="\t".join(map(str,info3))  
    RBHResultThres.write(info4)
    RBHResultThres.write("\n")

RBHResultThres.close()
