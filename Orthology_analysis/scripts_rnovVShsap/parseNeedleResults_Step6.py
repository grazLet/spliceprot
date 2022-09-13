import sys
import os
import itertools
import subprocess
import pickle
import pprint

# dos pares rejeitados pelo needle
# e aceitos pelo rbh (app1)
# quais os valores de similarity, identity, gap etc pelo needle?

database_list = ["ccds","ensembl","spliceprot","refseq"]
for database in database_list:
    # 1. Needle Final List
    #NeedleFinalListFile = open("analysis/RBH-Only_analise/"+database+"_FinalList-Needle.txt").read().split("\n")[:-1]
    FileOut = open("analysis/RBH-Only_analise/"+database+"_NeedleInfo_RBH1-Only-pairs.txt","w")
    FileOut2 = open("analysis/RBH-Only_analise/"+database+"_NeedleInfo_RBH2-Only-pairs.txt","w")
    FileOut3 = open("analysis/RBH-Only_analise/"+database+"_NeedleInfo_RBH1-Total-pairs.txt","w")
    FileOut4 = open("analysis/RBH-Only_analise/"+database+"_NeedleInfo_App1-Total-pairs.txt","w")
    FileOut5 = open("analysis/RBH-Only_analise/"+database+"_NeedleInfo_App2-Total-pairs.txt","w")
    FileOut6 = open("analysis/RBH-Only_analise/"+database+"_NeedleInfo_App1-Only-pairs.txt","w")
    FileOut7 = open("analysis/RBH-Only_analise/"+database+"_NeedleInfo_App2-Needle_RBH1.txt","w")
    FileOut8 = open("analysis/RBH-Only_analise/"+database+"_NeedleInfo_App2-Needle_RBH2.txt","w")
    # 2. Needle DB
    NeedleDBFile = open(database+"/needle/Needle_db/Needle_complete_db_table.txt").read().split("\n")[:-1]
    header = NeedleDBFile.pop(0)+"\n"
    FileOut.write(header)
    FileOut2.write(header)
    FileOut3.write(header)
    FileOut4.write(header)
    FileOut5.write(header)
    FileOut6.write(header)
    FileOut7.write(header)
    FileOut8.write(header)
    needleDict = {}
    hsDict = {}
    mmDict = {}
    for line in NeedleDBFile:
        col = line.split("\t")
        if database == "spliceprot":
            pair  = col[0]+","+col[1]
            info = "\t".join(col[2:])
            hs = col[0]
            mm = col[1]
        if database != "spliceprot":
            id_info1 = col[0].split("-")
            id_info2 = col[1].split("-")
            info = "\t".join(col[2:])
            if database == "refseq":
                pair = id_info1[2]+","+id_info2[2]
                hs = id_info1[2]
                mm = id_info2[2]
            if database == "ccds" or database == "ensembl":
                pair = id_info1[1]+","+id_info2[1]
                hs = id_info1[1]
                mm = id_info2[1]
        needleDict[pair] = info

    # 3. RBH (diamond)
    RBHFinalListFile = open("analysis/finalTab/"+database+"_FinalList-RBH.txt").read().split("\n")[:-1]
    App2FinalListFile = open("analysis/finalTab/"+database+"_FinalList-App2.txt").read().split("\n")[:-1]
    App1FinalListFile = open("analysis/finalTab/"+database+"_FinalList-App1.txt").read().split("\n")[:-1]
#    App2_N_FinalListFile = open("analysis/finalTab/"+database+"_FinalList_App2-Needle.txt").read().split("\n")[:-1]
    RBH1 = open("analysis/RBH-Only_analise/"+database+"_FinalList-RBH1-rbh2xrb1.txt").read().split("\n")[:-1]
    RBH2 = open("analysis/RBH-Only_analise/"+database+"_FinalList-RBH2-rbh2xrb1.txt").read().split("\n")[:-1]

    #App2 - needle intersection com rbh1
    for rbh1 in RBH1: 
        if rbh1 in needleDict.keys():
            rbh1a = rbh1.replace(",","\t")
            infoa = rbh1a+ "\t"+needleDict[rbh1]+"\n"
            FileOut7.write(infoa)
    FileOut7.close()

    for rbh2 in RBH2:
        if rbh2 in needleDict.keys():
            rbh2b = rbh2.replace(",","\t")
            infob = rbh2b+ "\t"+needleDict[rbh2]+"\n"
            FileOut8.write(infob)
    FileOut8.close()
    # RBH1-Only
    app1_rbh_exclusive = list(set(RBHFinalListFile).difference(set(App2FinalListFile)))
    for pair2 in app1_rbh_exclusive:
        if pair2 in needleDict.keys():
            p = pair2.replace(",","\t")
            info2 = p + "\t"+needleDict[pair2]+"\n"
            FileOut.write(info2)
    FileOut.close()
    # RBH2-Only
    app2_exclusive = list(set(App2FinalListFile).difference(set(App1FinalListFile)))
    for pair3 in app2_exclusive:
        if pair3 in needleDict.keys():
            p3 = pair3.replace(",","\t")
            info3 = p3 + "\t"+needleDict[pair3]+"\n"
            FileOut2.write(info3)
    FileOut2.close()
    
    # RBH1 - Total
    for pair6 in RBHFinalListFile:
        if pair6 in needleDict.keys():
            p6 = pair6.replace(",","\t")
            info6 = p6 + "\t"+needleDict[pair6]+"\n"
            FileOut3.write(info6)
    FileOut3.close()

    # total app1
    for pair4 in App1FinalListFile:
        if pair4 in needleDict.keys():
            p4 = pair4.replace(",","\t")
            info4 = p4 + "\t"+needleDict[pair4]+"\n"
            FileOut4.write(info4)
    FileOut4.close()

    # total app2
    for pair5 in App2FinalListFile:
        if pair5 in needleDict.keys():
            p5 = pair5.replace(",","\t")
            info5 = p5+ "\t"+needleDict[pair5]+"\n"
            FileOut5.write(info5)
    FileOut5.close()

    ## app1 exclusive
    ## app1 - app2
    app1_exclusive = list(set(App1FinalListFile).difference(set(App2FinalListFile)))
    for pairA in app1_exclusive:
        if pairA in needleDict.keys():
            pA = pairA.replace(",","\t")
            infoA = pA + "\t"+needleDict[pairA]+"\n"
            FileOut6.write(infoA)
    FileOut6.close()

