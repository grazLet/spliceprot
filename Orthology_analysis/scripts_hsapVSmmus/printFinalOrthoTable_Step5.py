import sys
import os
import itertools
import subprocess
import pickle
import pprint
########## Functions - organize datasets
def compareDataset2(dataset2FileName):
    dataset2File = open(dataset2FileName).read().split("\n")[1:-1]
    blanqList = []
    for line in dataset2File:
        hsCCDS,mmCCDS = line.split()
        hsCCDS = hsCCDS.strip()
        mmCCDS = mmCCDS.strip()
        info2 = hsCCDS+"\t"+mmCCDS
        if info2 not in blanqList:
            blanqList.append(info2)
    return(blanqList)

def compareDataset4(Dataset4FileName):
    Dataset4File = open(Dataset4FileName).read().split("\n")[:-1]
    dataset4List = []

    for line in Dataset4File:
        col = line.split()
        hsaplist = col[2].split(",")
        mmuslist = col[5].split(",")
        for hs,mm in itertools.product(hsaplist,mmuslist):
            pair = hs+"\t"+mm
            if pair not in dataset4List:
                dataset4List.append(pair)
    return(dataset4List)

def compareDataset3(Dataset3FileName):
    Dataset3File = open(Dataset3FileName).read().split("\n")[:-1]
    dataset3List = []

    for line in Dataset3File:
        col = line.split()
        geneName = col[0].strip()
        hsT = col[1].strip()
        mmT = col[2].strip()
        pair = hsT+"\t"+mmT
        if pair not in dataset3List:
            dataset3List.append(pair)
    return(dataset3List)

### Files from datasets -  for validation
fileDs2 = "/home/flavia.freitas/disk1/flavia.freitas/results/validation/dataset2_blanquart_ccds_ortho_pairs.txt"
fileDs4 = "/home/flavia.freitas/disk1/flavia.freitas/results/validation/dataset4_ensembl_ortho_pairs.txt"
fileDs3 = "/home/flavia.freitas/disk1/flavia.freitas/results/validation/dataset3_refseq_zambelli_ortho_pairs.txt"

# ccds
dataset2 = set(compareDataset2(fileDs2))
# refseq
dataset3 = set(compareDataset3(fileDs3))
# ensembl
dataset4 = set(compareDataset4(fileDs4))

FileListDs2 = open("analysis/validation/Dataset2_ccds_Blanquart_list.txt","w")
for p2 in list(dataset2):
    p2 = p2.replace("\t",",")+"\n"
    FileListDs2.write(p2)
FileListDs2.close()

FileListDs3 = open("analysis/validation/Dataset3_refseq_Zambelli_list.txt","w")
for p3 in list(dataset3):
    p3 = p3.replace("\t",",")+"\n"
    FileListDs3.write(p3)
FileListDs3.close()

FileListDs4 = open("analysis/validation/Dataset4_ensemb_R-Sole_list.txt","w")
for p4 in list(dataset4):
    p4 = p4.replace("\t",",")+"\n"
    FileListDs4.write(p4)
FileListDs4.close()


pickleOrthoFile = "/home/flavia.freitas/disk1/flavia.freitas/data/ortho_ensembl_101/Table_ID_Hsap_Mmus_Ortho_1-1.pickle"

def loadDict (pickleOrthoFile):
    pickleFile = open(pickleOrthoFile,"rb")
    OrthoDict = pickle.load(pickleFile)
    alignD = {}

    for hsapID,mmusID in OrthoDict.items():
        genePair = hsapID+"\t"+mmusID
        alignD[genePair] = {"needleTr":[],"needleInfo":[],"rbhTr":[],"rbhInfo":[],"method":[],"app2Tr":[],"app2Info":[]}

    return(alignD)

database_list = ["ccds","ensembl","spliceprot","refseq"]
#database_list = ["ccds","spliceprot","refseq"]


rbhNA = "\t".join(["NA"]*8)
rbhNA2 = "\t".join(["NA"]*2)
needleNA = "\t".join(["NA"]*7)
needleNA2 = "\t".join(["NA"]*2)

os.mkdir("analysis")
os.mkdir("analysis/finalTab")
os.mkdir("analysis/app1_rbh")
os.mkdir("analysis/app2_rbh")
os.mkdir("analysis/needle_resultados")

OutFileCompare4 = open("analysis/finalTab/App1_App2_dataframe_alldb.txt" ,"w")
x11 = "\t".join(["database","app1_needle","app1_rbh","app1","app2","app1_intersection_app2","app1_diff_app2", \
                    "app2_diff_app1","app1_symdiff_app2","needle_rbh2","rbh_rbh2","needle_dif_rbh2", \
                    "rbh2_dif_needle","rbh_dif_rbh2","rbh2_dif_rbh","\n"])

OutFileCompare4.write(x11)
OutFileCompare3 = open("analysis/finalTab/App1_App2_CompleteComparisonDB.txt" ,"w")
dbDict = {}
for database in database_list:
    FileOutFinalListNeedle = open("analysis/finalTab/"+database+"_FinalList-Needle.txt","w") 
    FileOutFinalListApp2_N = open("analysis/finalTab/"+database+"_FinalList_App2-Needle.txt","w")
    FileOutFinalListRBH = open("analysis/finalTab/"+database+"_FinalList-RBH.txt","w")
    FileOutFinalListApp2 = open("analysis/finalTab/"+database+"_FinalList-App2.txt","w")
    FileOutFinalListApp1 = open("analysis/finalTab/"+database+"_FinalList-App1.txt","w")
    FileOutRBH1 = open("analysis/RBH-Only_analise/"+database+"_FinalList-RBH2-rbh2xrb1.txt","w") 
    FileOutRBH2 = open("analysis/RBH-Only_analise/"+database+"_FinalList-RBH1-rbh2xrb1.txt","w")
    RBHOnlyFile = open("analysis/finalTab/"+database+"_RBHOnly.txt","w")
    a11 = "".join([">>> database:", database,"\n"])
    OutFileTab = open("analysis/finalTab/"+database+"_FinalTab_App1.txt" ,"w")
    
    # load ortho genes dict
    alignD = loadDict(pickleOrthoFile)

    # 1. NEEDLE
    NeedleOutFilePath = "analysis/needle_resultados/"
    
    # 1.1 needle DB
    NeedleDBFile = open(database+"/needle/Needle_db/Needle_complete_db_table.txt").read().split("\n")[:-1]
    NeedleDBOutFileName = open(NeedleOutFilePath+database+"_NeedleDB.txt","w")
    for line in NeedleDBFile:
        line = line.replace("#","")
        NeedleDBOutFileName.write(line)
        NeedleDBOutFileName.write("\n")
    NeedleDBOutFileName.close()

    # 1.2 needle filtered
    NeedleFilterFile = open(database+"/needle/align_results/OrthoAnalysis_table.txt").read().split("\n")[:-1]
    NeedleFilterOutFileName = open(NeedleOutFilePath+database+"_needle_filtered.txt","w")
    needleHeader = NeedleFilterFile[0]
    needleHeaderItens = needleHeader.split("\t")

    for line2 in NeedleFilterFile:
        NeedleFilterOutFileName.write(line2)
        NeedleFilterOutFileName.write("\n")

        if NeedleFilterFile.index(line2) != 0:
            col2 = line2.split("\t")
            if database=="spliceprot":
                hsTr = col2[0].strip()
                mmTr = col2[1].strip()
                hsgene,var1 = hsTr.split("-",1)
                mmgene,var2 = mmTr.split("-",1)
 
            if database!="spliceprot":
                hsInfo = col2[0].split("-")
                hsgene = hsInfo[0]
                mmInfo = col2[1].split("-")
                mmgene = mmInfo[0]
                if database!="ensembl":
                    hsTr2 = hsInfo[-1]
                    mmTr2 = mmInfo[-1]
                if database == "ensembl":
                    hsTr2 = hsInfo[1]
                    mmTr2 = mmInfo[1]
                if "." in hsTr2:
                    hsTr,lixo = hsTr2.split(".")
                    mmTr,lixo2 = mmTr2.split(".")
                if "." not in hsTr2:
                    hsTr = hsTr2
                    mmTr = mmTr2
            needleInfo = "\t".join(col2[2:])
            pair=hsTr + "\t"+mmTr
            pair2 = hsgene+"\t"+mmgene
            appr="needle"
            if pair2 in alignD.keys():
                alignD[pair2]["needleTr"].append(pair)
                alignD[pair2]["needleInfo"].append(needleInfo)
                alignD[pair2]["rbhTr"].append(rbhNA2)
                alignD[pair2]["rbhInfo"].append(rbhNA)
                alignD[pair2]["method"].append(appr)

    NeedleFilterOutFileName.close()

    # 2. APP1 - RBH
    App1OutFilePath = "analysis/app1_rbh/"
    
    # 2.1 app1 Total
    App1TotalFile = open(database+"/diamond_needle/RBH/hsap2/mmus2.rbh").read().split("\n")[:-1]
    App1TotalOutFile = open(App1OutFilePath+database+"_app1_RBH_total.txt","w")
    for line3 in App1TotalFile:
        line3 = line3.replace("#","").replace("Hsap_","").replace("Mmus_","")
        App1TotalOutFile.write(line3)
        App1TotalOutFile.write("\n")
    App1TotalOutFile.close()

    # 2.2 app1 filtered tcover/qcover ranging from 60 to 90
    App1File = open(database+"/diamond_needle/analysis_RBH/"+database+"_RBH_ortho_60-90_tab.txt").read().split("\n")[:-1]
    App1OutFile = open(App1OutFilePath+database+"_app1_RBH_filtered.txt","w")
    rbhHeader = App1File[0]
    rbhHeaderItens = rbhHeader.split("\t")

    for line4 in App1File:
        App1OutFile.write(line4)
        App1OutFile.write("\n")

        if App1File.index(line4) != 0:
            col4 = line4.split("\t")
            hsgene2 = col4[0].strip()
            mmgene2 = col4[1].strip()
            hsTr2 = col4[2].strip()
            mmTr2 = col4[3].strip()
            appr2 = "rbh"
            rbhInfo = "\t".join(col4[4:-1])
            pair_1 = hsTr2+"\t"+mmTr2
            pair_2 = hsgene2 +"\t"+ mmgene2

            if pair_2 in alignD.keys():
                alignD[pair_2]["rbhTr"].append(pair_1)
                alignD[pair_2]["rbhInfo"].append(rbhInfo)
                alignD[pair_2]["needleTr"].append(needleNA2)
                alignD[pair_2]["needleInfo"].append(needleNA)
                alignD[pair_2]["method"].append(appr2)

    App1OutFile.close()

    headerFinal = "\t".join(["HsGeneID","MmGeneID","\t".join(needleHeaderItens),"\t".join(rbhHeaderItens[2:]),"method"])
    OutFileTab.write(headerFinal)
    OutFileTab.write("\n")
   
    # App1 Final Tab
    for k,v in alignD.items():
        if v["rbhTr"] or v["needleTr"] != [] :
            for i in range(len(v["method"])):
                aa = "\t".join([k,v["needleTr"][i],v["needleInfo"][i],v["rbhTr"][i],v["rbhInfo"][i],v["method"][i]])
                OutFileTab.write(aa)
                OutFileTab.write("\n")
    OutFileTab.close()


    # Part II Comparison
    # 1. App1 x App2 

    #OutFileCompare = open("analysis/finalTab/"+database+"_db_compare_App1_App2.txt" ,"w")
    
    # 1.1 Load app2 total
    App2OutFilePath = "analysis/app2_rbh/"
    App2TotalFile = open(database+"/app2_diamond/RBH/hsap/mmus.rbh").read().split("\n")[:-1]
    App2TotalOutFile = open(App2OutFilePath+database+"_app2_RBH_total.txt","w")

    for line5 in App2TotalFile:
        line5 = line5.replace("#","")
        App2TotalOutFile.write(line5)
        App2TotalOutFile.write("\n")
    App2TotalOutFile.close()

    App2File = open(database+"/app2_diamond/analysis_RBH/"+database+"_RBH_ortho_60-90_tab.txt").read().split("\n")[:-1]
    App2FileOutFile = open(App2OutFilePath+database+"_app2_RBH_filtered.txt","w")
    for line6 in App2File:
        line6 = line6.replace("#","")
        App2FileOutFile.write(line6)
        App2FileOutFile.write("\n")
        if App2File.index(line6) != 0:
            col6 = line6.split("\t")
            hsG3 = col6[0]
            mmG3 = col6[1]
            hsTr3 = col6[2]
            mmTr3 = col6[3]
            info3 = col6[4:]
            pair_a = hsTr3+"\t"+mmTr3
            pair_b = hsG3 +"\t"+ mmG3
            if pair_b in alignD.keys():
               alignD[pair_b]["app2Tr"].append(pair_a)
               alignD[pair_b]["app2Info"].append(info3)

    App2FileOutFile.close()

    # lista de números
    needle_list = []
    rbh_list = []
    app1_list = []
    app2_list = []
    app1_2_list = []
    app2_dif_app1_list = []
    app1_dif_app2_list = []
    app1_symdiff_app2_list = []
    needle_rbh2_list = []
    rbh_rbh2_list = []
    needle_dif_rbh2_list = []
    rbh2_dif_needle_list = []
    rbh_dif_rbh2_list = []
    rbh2_dif_rbh_list = []
    rbh2_needle_inter_rbh1_list = []
    n_dif_rbh2_dif_rbh1_rbh2_app1_list = []
    n_dif_rbh2_dif_rbh1_rbh2_app2_list = []
    # lista de transcritos
    rbh2_needle_inter_rbh1_list_2 =[]
    needle_list_2 = []
    rbh_list_2 = []
    app1_list_2 = []
    app2_list_2 = []
    app1_2_list_2 = []
    app2_dif_app1_list_2 = []
    app1_dif_app2_list_2 = []
    app1_symdiff_app2_list_2 = []
    needle_rbh2_list_2 = []
    rbh_rbh2_list_2 = []
    needle_dif_rbh2_list_2 = []
    rbh2_dif_needle_list_2 = []
    rbh_dif_rbh2_list_2 = []
    rbh2_dif_rbh_list_2 = []
    n_dif_rbh2_dif_rbh1_rbh2_app1_list_2 = []
    n_dif_rbh2_dif_rbh1_rbh2_app2_list_2 = []
    # lista id modificados
    n_dif_rbh2_dif_rbh1_rbh2_app1_list_3 = []
    n_dif_rbh2_dif_rbh1_rbh2_app2_list_3 = []
    rbh2_needle_inter_rbh1_list_3 = []
    needle_list_3 = []
    rbh_list_3 = []
    app1_list_3 = []
    app2_list_3 = []
    app1_2_list_3 = []
    app2_dif_app1_list_3 = []
    app1_dif_app2_list_3 = []
    app1_symdiff_app2_list_3 = []
    needle_rbh2_list_3 = []
    rbh_rbh2_list_3 = []
    needle_dif_rbh2_list_3 = []
    rbh2_dif_needle_list_3 = []
    rbh_dif_rbh2_list_3 = []
    rbh2_dif_rbh_list_3 = []
    for x,y in alignD.items():
        needle2 = y["needleTr"] #.remove("NA\tNA")
        rbh2 = y["rbhTr"]
        app2_2 = y["app2Tr"]
        if "NA\tNA" in needle2:
            while "NA\tNA" in needle2:
                needle2.remove("NA\tNA")
        if "NA\tNA" in rbh2:
            while "NA\tNA" in rbh2:
                rbh2.remove("NA\tNA")
        if "NA\tNA" in app2_2:
            while "NA\tNA" in app2_2:
                app2_2.remove("NA\tNA")
    
        needle = set(filter(None,needle2))
        rbh = set(filter(None,rbh2))
        app2 = set(filter(None,app2_2))
        # app1 = union needle e rbh
        app1 = set(needle.union(rbh))
        app1_intersection_app2 = list(app1.intersection(app2))
        needle_rbh2 = list(needle.intersection(app2))
        rbh_rbh2 = list(rbh.intersection(app2))
        #### abordagem que comparar app1 e app2 
        needle_diff_rbh2 = list(needle.difference(app2))
        rbh2_diff_needle = list(app2.difference(needle))
        n_dif_rbh2_dif_rbh1_rbh2_app1 = list(set(rbh).difference(set(rbh_rbh2)))
        n_dif_rbh2_dif_rbh1_rbh2_app2 = list(set(rbh2_diff_needle).difference(set(rbh_rbh2)))

        rbh2_needle_inter_rbh1 = list(set(rbh2_diff_needle).intersection(set(rbh)))
        rbh_diff_rbh2 = list(rbh.difference(app2))
        rbh2_diff_rbh = list(app2.difference(rbh))
        app1_symdiff_app2 = list(app1.symmetric_difference(app2))
        #ver qualidade dos pares identificados apenas em app2
        app2_diff_app1 = list(app2.difference(app1))
        app1_diff_app2 = list(app1.difference(app2))
        for i in app1:
            app1_list_2.append(i)
            ii = i.replace("\t",",")
            app1_list_3.append(ii)
            FileOutFinalListApp1.write(ii)
            FileOutFinalListApp1.write("\n")
        for i2 in app2:
            app2_list_2.append(i2)
            ii2 = i2.replace("\t",",")
            app2_list_3.append(ii2)
            FileOutFinalListApp2.write(ii2)
            FileOutFinalListApp2.write("\n")
        for i3 in needle:
            needle_list_2.append(i3)
            ii3 = i3.replace("\t",",")
            needle_list_3.append(ii3)
            FileOutFinalListNeedle.write(ii3)
            FileOutFinalListNeedle.write("\n")
        for i4 in rbh:
            rbh_list_2.append(i4)
            ii4 = i4.replace("\t",",")
            rbh_list_3.append(ii4)
            FileOutFinalListRBH.write(ii4)
            FileOutFinalListRBH.write("\n")
        for i5 in app1_intersection_app2:
            app1_2_list_2.append(i5)
            ii5 = i5.replace("\t",",")
            app1_2_list_3.append(ii5)
        for i6 in app1_symdiff_app2:
            app1_symdiff_app2_list_2.append(i6)
            ii6 = i6.replace("\t",",")
            app1_symdiff_app2_list_3.append(ii6)
        for i7 in app2_diff_app1:
            app2_dif_app1_list_2.append(i7)
            ii7 = i7.replace("\t",",")
            app2_dif_app1_list_3.append(ii7)
        for i8 in app1_diff_app2:
            app1_dif_app2_list_2.append(i8)
            ii8 = i8.replace("\t",",")
            app1_dif_app2_list_3.append(ii8)
        for i9 in needle_rbh2:
            needle_rbh2_list_2.append(i9)
            ii9 = i9.replace("\t",",")
            needle_rbh2_list_3.append(ii9)
        for i10 in rbh_rbh2:
            rbh_rbh2_list_2.append(i10)
            ii10 = i10.replace("\t",",")
            rbh_rbh2_list_3.append(ii10)
        for i11 in needle_diff_rbh2:
            needle_dif_rbh2_list_2.append(i11)
            ii11 = i11.replace("\t",",")
            needle_dif_rbh2_list_3.append(ii11)
            iii11 = ii11+"\n"
            FileOutFinalListApp2_N.write(iii11)
        for i12 in rbh2_diff_needle:
            rbh2_dif_needle_list_2.append(i12)
            ii12 = i12.replace("\t",",")
            rbh2_dif_needle_list_3.append(ii12)
        for i13 in rbh_diff_rbh2:
            rbh_dif_rbh2_list_2.append(i13)
            ii13 = i13.replace("\t",",")
            rbh_dif_rbh2_list_3.append(ii13)
        for i14 in rbh2_diff_rbh:
            rbh2_dif_rbh_list_2.append(i14)
            ii14 = i14.replace("\t",",")
            rbh2_dif_rbh_list_3.append(ii14)

        for i15 in rbh2_needle_inter_rbh1:
            rbh2_needle_inter_rbh1_list_2.append(i15)
            ii15 = i15.replace("\t",",")
            #iii15 = ii15+"\n"
            rbh2_needle_inter_rbh1_list_3.append(ii15)
            #FileOutFinalListApp2_N.write(iii15)
        for i16 in n_dif_rbh2_dif_rbh1_rbh2_app1:
            n_dif_rbh2_dif_rbh1_rbh2_app1_list_2.append(i16)
            ii16 = i16.replace("\t",",")
            iii16 = ii16+"\n"
            n_dif_rbh2_dif_rbh1_rbh2_app1_list_3.append(ii16)
            FileOutRBH1.write(iii16)

        for i17 in n_dif_rbh2_dif_rbh1_rbh2_app2:
            n_dif_rbh2_dif_rbh1_rbh2_app2_list_2.append(i17)
            ii17 = i17.replace("\t",",")
            iii17 = ii17+"\n"
            n_dif_rbh2_dif_rbh1_rbh2_app2_list_3.append(ii17)
            FileOutRBH2.write(iii17)

        needle_rbh2_list.append(len(list(needle_rbh2)))
        rbh_rbh2_list.append(len(list(rbh_rbh2)))
        needle_list.append(len(list(needle)))
        rbh_list.append(len(list(rbh)))
        app1_2_list.append(len(app1_intersection_app2))
        app2_list.append(len(list(app2)))
        app1_list.append(len(list(app1)))
        app2_dif_app1_list.append(len(app2_diff_app1))
        app1_dif_app2_list.append(len(app1_diff_app2))
        app1_symdiff_app2_list.append(len(app1_symdiff_app2))

    hsApp1List = []
    mmApp1List = [] 
    for app1Pair in app1_list_3:
        hs,mm = app1Pair.split(",")
        hsApp1List.append(hs)
        mmApp1List.append(mm)
    
    print(">>>>>", database)
    for rbhPair in app2_dif_app1_list_3:
        RBHOnlyFile.write(rbhPair)
        RBHOnlyFile.write("\n")
        hs1,mm1 = rbhPair.split(",")
        if hs1 not in hsApp1List and mm1 not in mmApp1List:
            print(rbhPair)
    

    dbDict[database] = {"rbh2_needle_inter_rbh1_list_2":rbh2_needle_inter_rbh1_list_2,
                        "needle_list_2": needle_list_2,
                        "rbh_list_2": rbh_list_2,
                        "app1_list_2": app1_list_2,
                        "app2_list_2": app2_list_2,
                        "app1_2_list_2": app1_2_list_2,
                        "app2_dif_app1_list_2": app2_dif_app1_list_2,
                        "app1_dif_app2_list_2": app1_dif_app2_list_2,
                        "app1_symdiff_app2_list_2": app1_symdiff_app2_list_2,
                        "needle_rbh2_list_2":needle_rbh2_list_2,
                        "rbh_rbh2_list_2":rbh_rbh2_list_2,
                        "needle_dif_rbh2_list_2":needle_dif_rbh2_list_2,
                        "rbh2_dif_needle_list_2":rbh2_dif_needle_list_2,
                        "rbh_dif_rbh2_list_2":rbh_dif_rbh2_list_2,
                        "rbh2_dif_rbh_list_2 ":rbh2_dif_rbh_list_2}
    
    RBHOnlyFile.close()
    a1 = "\t".join([">>>",database,"\n"])
    a2 = "\t".join(["app1:",str(len(app1_list_3)),"\n"])
    a22 = "".join([";".join(app1_list_3),"\n","#############","\n"])
    a3 = "\t".join(["app2:",str(len(app2_list_3)),"\n"])
    a33 ="".join([";".join(app2_list_3),"\n","#############","\n"])
    a4 = "\t".join(["app1_intersection_app2:",str(len(app1_2_list_3)),"\n"])
    a44 ="".join([";".join(app1_2_list_3),"\n","#############","\n"])
    a5 = "\t".join(["app1_diff_app2:",str(len(app1_dif_app2_list_3)),"\n"])
    a55 ="".join([";".join(app1_dif_app2_list_3),"\n","#############","\n"])
    a6 = "\t".join(["app2_diff_app1:",str(len(app2_dif_app1_list_3)),"\n"])
    a66 = "".join([";".join(app2_dif_app1_list_3),"\n","#############","\n"])
    a7 = "\t".join(["app1_symdiff_app2:",str(len(app1_symdiff_app2_list_3)),"\n"])
    a77 = "".join([";".join(app1_symdiff_app2_list_3),"\n","#############","\n"])
    a8 = "\t".join(["needle_rbh2:",str(len(needle_rbh2_list_3)),"\n"])
    a88 = "".join([";".join(needle_rbh2_list_3),"\n","#############","\n"])
    a9 = "\t".join(["rbh_rbh2:",str(len(rbh_rbh2_list_3)),"\n"])
    a99 ="".join([";".join(rbh_rbh2_list_3),"\n","#############","\n"])
    a10 ="\t".join(["needle_diff_rbh2:",str(len(needle_dif_rbh2_list_3)),"\n"])
    a110 ="".join([";".join(needle_dif_rbh2_list_3),"\n","#############","\n"])
    a11 = "\t".join(["rbh2_diff_needle:",str(len(rbh2_dif_needle_list_3)),"\n"])
    a111 ="".join([";".join(rbh2_dif_needle_list_3),"\n","#############","\n"])
    a12 ="\t".join(["rbh_diff_rbh2:",str(len(rbh_dif_rbh2_list_3)),"\n"])
    a112 = "".join([";".join(rbh_dif_rbh2_list_3),"\n","#############","\n"])
    a13 = "\t".join(["rbh2_diff_rbh:",str(len(rbh2_dif_rbh_list_3)),"\n"])
    a113 = "".join([";".join(rbh2_dif_rbh_list_3),"\n","#############","\n"])
    a14 = "\t".join(["rbh2_needle_inter_rbh1:",str(len(rbh2_needle_inter_rbh1_list_3)),"\n"])
    a114 = "".join([";".join(rbh2_needle_inter_rbh1_list_3),"\n","#############","\n"])
    a20 = "-------------------------------------------------------\n"
    OutFileCompare3.write(a1)
    OutFileCompare3.write(a2)
    OutFileCompare3.write(a22)
    OutFileCompare3.write(a3)
    OutFileCompare3.write(a33)
    OutFileCompare3.write(a4)
    OutFileCompare3.write(a44)
    OutFileCompare3.write(a5)
    OutFileCompare3.write(a55)
    OutFileCompare3.write(a6)
    OutFileCompare3.write(a66)
    OutFileCompare3.write(a7)
    OutFileCompare3.write(a77)
    OutFileCompare3.write(a8)
    OutFileCompare3.write(a88)
    OutFileCompare3.write(a9)
    OutFileCompare3.write(a99)
    OutFileCompare3.write(a10)
    OutFileCompare3.write(a110)
    OutFileCompare3.write(a11)
    OutFileCompare3.write(a111)
    OutFileCompare3.write(a12)
    OutFileCompare3.write(a112)
    OutFileCompare3.write(a13)
    OutFileCompare3.write(a113)
    OutFileCompare3.write(a14)
    OutFileCompare3.write(a114)
    OutFileCompare3.write(a20)

    x3 = "\t".join([database,str(len(needle_list_3)),str(len(rbh_list_3)),str(len(app1_list_3)),str(len(app2_list_3)),str(len(app1_2_list_3)), \
                                        str(len(app1_dif_app2_list_3)),str(len(app2_dif_app1_list_3)), str(len(app1_symdiff_app2_list_3)), \
                                        str(len(needle_rbh2_list_3)),str(len(rbh_rbh2_list_3)),str(len(needle_dif_rbh2_list_3)), \
                                        str(len(rbh2_dif_needle_list_3)),str(len(rbh_dif_rbh2_list_3)),str(len(rbh2_dif_rbh_list_3)),"\n"])
    
    OutFileCompare4.write(x3)
    FileOutRBH1.close()
    FileOutRBH2.close()
    FileOutFinalListApp1.close()
    FileOutFinalListApp2.close()
    FileOutFinalListNeedle.close()
    FileOutFinalListRBH.close()
    FileOutFinalListApp2_N.close()
OutFileCompare3.close()
OutFileCompare4.close()

DictFile = open("analysis/finalTab/dbDict_App1-App2_Comparison.txt","w")
pprint.pprint(dbDict,DictFile)
DictFile.close()

## comparar com dados datasets
FileCountValid = open("analysis/validation/Result_intersection_app_datasets.txt","w")
c11 = "\t".join(["database","dataset","app1","intersection_app1_dset","app2","intersection_app2_dset","\n"])
FileCountValid.write(c11)
for db,info in dbDict.items():
    if "ccds" in db:
        ds2_app1_2 = []
        ds2_app2_2 = []
        FileOutValid = open("analysis/validation/dataset2_ccds_app1.txt","w")
        FileOutValid2 = open("analysis/validation/dataset2_ccds_app2.txt","w")
        app1_list_2 = info["app1_list_2"]
        app2_list_2 = info["app2_list_2"]
        ds2_app1 = list(dataset2.intersection(set(app1_list_2)))
        ds2_app2 = list(dataset2.intersection(set(app2_list_2)))
        for pair1 in ds2_app1:
            pair11 = pair1.replace("\t",",")
            ds2_app1_2.append(pair11)
        for pair2 in ds2_app2:
            pair22 = pair2.replace("\t",",") 
            ds2_app2_2.append(pair22)
        b1 = "\n".join(ds2_app1_2)
        b2 = "\n".join(ds2_app2_2)
        h1 = "##### dataset2.intersection.app1 ####\n"
        h2 = "##### dataset2.intersection.app2 ####\n"
        FileOutValid.write(h1)
        FileOutValid.write(b1)
        #FileOutValid.write("\n----------------------------------------------\n")
        FileOutValid2.write(h2)
        FileOutValid2.write(b2)
        #FileOutValid.write("\n----------------------------------------------\n")
            
        c1 = str(len(app1_list_2))
        c2 = str(len(app2_list_2))
        c3 = str(len(ds2_app1))
        c4 = str(len(ds2_app2))
        c5 = str(len(dataset2))
        #c11 = "\t".join(["database","dataset2","app1","intersection_app1_ds2","app2","intersection_app2_ds2","\n"])
        c12 = "\t".join([db,c5,c1,c3,c2,c4,"\n"])

        #FileCountValid.write(c11)
        FileCountValid.write(c12)
        FileOutValid.close()   
        FileOutValid2.close()
    if "refseq" in db:
        ds3_app1_2 = []
        ds3_app2_2 = []
        FileOutValid = open("analysis/validation/dataset3_refseq_app1.txt","w")
        FileOutValid2 = open("analysis/validation/dataset3_refseq_app2.txt","w")
        app1_list_2 = info["app1_list_2"]
        app2_list_2 = info["app2_list_2"]
        ds3_app1 = list(dataset3.intersection(app1_list_2))
        ds3_app2 = list(dataset3.intersection(app2_list_2))
        for pair1 in ds3_app1:
            pair11 = pair1.replace("\t",",")
            ds3_app1_2.append(pair11)
        for pair2 in ds3_app2:
            pair22 = pair2.replace("\t",",")
            ds3_app2_2.append(pair22)
        b3 = "\n".join(ds3_app1_2)
        b4 = "\n".join(ds3_app2_2)
        h1 = "##### dataset3.intersection.app1 ####\n"
        h2 = "##### dataset3.intersection.app2 ####\n"
        FileOutValid.write(h1)
        FileOutValid.write(b3)
        #FileOutValid.write("\n----------------------------------------------\n")
        FileOutValid2.write(h2)
        FileOutValid2.write(b4)
        #FileOutValid.write("\n----------------------------------------------\n")

        d1 = str(len(app1_list_2))
        d2 = str(len(app2_list_2))
        d3 = str(len(ds3_app1))
        d4 = str(len(ds3_app2))
        d5 = str(len(dataset3))
        #d11 = "\t".join(["database","dataset3","app1","intersection_app1_ds2","app2","intersection_app2_ds2","\n"])
        d12 = "\t".join([db,d5,d1,d3,d2,d4,"\n"])
        #FileCountValid.write(c11)
        FileCountValid.write(d12)
        FileOutValid.close()
        FileOutValid2.close()
    if "ensembl" in db:
        ds4_app1_2 = []
        ds4_app2_2 = []
        FileOutValid = open("analysis/validation/dataset4_ensembl_app1.txt","w")
        FileOutValid2 = open("analysis/validation/dataset4_ensembl_app2.txt","w")
        app1_list_2 = info["app1_list_2"]
        app2_list_2 = info["app2_list_2"]
        ds4_app1 = list(dataset4.intersection(app1_list_2))
        ds4_app2 = list(dataset4.intersection(app2_list_2))
        for pair1 in ds4_app1:
            pair11 = pair1.replace("\t",",")
            ds4_app1_2.append(pair11)
        for pair2 in ds4_app2:
            pair22 = pair2.replace("\t",",")
            ds4_app2_2.append(pair22)
        b5 = "\n".join(ds4_app1)
        b6 = "\n".join(ds4_app2)
        h1 = "##### dataset4.intersection.app1 ####\n"
        h2 = "##### dataset4.intersection.app2 ####\n"
        FileOutValid.write(h1)
        FileOutValid.write(b5)
        #FileOutValid.write("\n----------------------------------------------\n")
        FileOutValid2.write(h2)
        FileOutValid2.write(b6)
        #FileOutValid.write("\n----------------------------------------------\n")

        d1 = str(len(app1_list_2))
        d2 = str(len(app2_list_2))
        d3 = str(len(ds4_app1))
        d4 = str(len(ds4_app2))
        d5 = str(len(dataset4))
        #d11 = "\t".join(["database","dataset4","app1","intersection_app1_ds2","app2","intersection_app2_ds2","\n"])
        d13 = "\t".join([db,d5,d1,d3,d2,d4,"\n"])
        #FileCountValid.write(c11)
        FileCountValid.write(d13)
        FileOutValid.close()
        FileOutValid2.close()
FileCountValid.close()


