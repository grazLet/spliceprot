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
    
    FileOut = open("analysis/RBH-Only_analise/"+database+"_RBHInfo_RBH1-Only-pairs.txt","w")
    FileOut2 = open("analysis/RBH-Only_analise/"+database+"_RBHInfo_RBH2-Only-pairs.txt","w")
    FileOut6 = open("analysis/RBH-Only_analise/"+database+"_RBHInfo_App1-Only-pairs.txt","w")
    FileOut7 = open("analysis/RBH-Only_analise/"+database+"_RBHInfo_App2-Needle_RBH1.txt","w")
    FileOut8 = open("analysis/RBH-Only_analise/"+database+"_RBHInfo_App2-Needle_RBH2.txt","w")
    # RBH (diamond)
    RBHFinalListFile = open("analysis/finalTab/"+database+"_FinalList-RBH.txt").read().split("\n")[:-1]
    App2FinalListFile = open("analysis/finalTab/"+database+"_FinalList-App2.txt").read().split("\n")[:-1]
    App1FinalListFile = open("analysis/finalTab/"+database+"_FinalList-App1.txt").read().split("\n")[:-1]
    RBH1 = open("analysis/RBH-Only_analise/"+database+"_FinalList-RBH1-rbh2xrb1.txt").read().split("\n")[:-1]
    RBH2 = open("analysis/RBH-Only_analise/"+database+"_FinalList-RBH2-rbh2xrb1.txt").read().split("\n")[:-1]

    RBH1File = open(database+"/diamond_needle/analysis_RBH/"+database+"_RBH_ortho_60-90_tab.txt").read().split("\n")[:-1]
    headerItens = RBH1File.pop(0).split("\t")
    headerItens.pop(0)
    headerItens.pop(0)
    header = "\t".join(headerItens)+"\n"
    RBH2File = open(database+"/app2_diamond/analysis_RBH/"+database+"_RBH_ortho_60-90_tab.txt").read().split("\n")[1:-1]
    rbh1Dict = {}
    rbh2Dict = {}
    
    for line in RBH1File:
        col = line.split("\t")
        pair  = col[2]+","+col[3]
        info = "\t".join(col[4:])
        rbh1Dict[pair] = info

    for line2 in RBH2File:
        col2 = line2.split("\t")
        pair8 = col2[2]+","+col2[3]
        info8 = "\t".join(col2[4:])
        rbh2Dict[pair8] = info8
    
    # RBH1-Only
    FileOut.write(header)
    FileOut2.write(header)
    FileOut6.write(header)
    FileOut7.write(header)
    FileOut8.write(header)
    for rbh1 in rbh1Dict.keys():
        rbh1a = rbh1.replace(",","\t")
        infoa = rbh1a+ "\t"+rbh1Dict[rbh1]+"\n"
        FileOut7.write(infoa)
    FileOut7.close()
    for rbh2 in rbh2Dict.keys():
        rbh2b = rbh2.replace(",","\t")
        infob = rbh2b+ "\t"+rbh2Dict[rbh2]+"\n"
        FileOut8.write(infob)
    FileOut8.close()

    app1_rbh_exclusive = list(set(RBHFinalListFile).difference(set(App2FinalListFile)))

    for pair2 in app1_rbh_exclusive:
        if pair2 in rbh1Dict.keys():
            p = pair2.replace(",","\t")
            info2 = p + "\t"+rbh1Dict[pair2]+"\n"
            FileOut.write(info2)
    FileOut.close()
    # RBH2-Only
    app2_exclusive = list(set(App2FinalListFile).difference(set(App1FinalListFile)))
    for pair3 in app2_exclusive:
        if pair3 in rbh2Dict.keys():
            p3 = pair3.replace(",","\t")
            info3 = p3 + "\t"+rbh2Dict[pair3]+"\n"
            FileOut2.write(info3)
    FileOut2.close()

    app1_exclusive = list(set(App1FinalListFile).difference(set(App2FinalListFile)))
    for pairA in app1_exclusive:
        if pairA in rbh1Dict.keys():
            pA = pairA.replace(",","\t")
            infoA = pA + "\t"+rbh1Dict[pairA]+"\n"
            FileOut6.write(infoA)
    FileOut6.close()
    

