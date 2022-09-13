import sys, os
import itertools
import subprocess


# estrutura de diretorio
# diamond e needle
# diamond/ contem resultados de toda a base de dados 
# needle/ e diamond_needle/

# usar arquivos pos-analise do RBH_analysis.py
# q filtra os pares de genes ortologos


def compareDataset2(dataset2FileName):
    dataset2File = open(dataset2FileName).read().split("\n")[1:-1]
    blanqList = []
    for line in dataset2File:
        hsCCDS,mmCCDS = line.split()
        hsCCDS = hsCCDS.strip()
        mmCCDS = mmCCDS.strip()
        info2 = hsCCDS+"-"+mmCCDS
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
            pair = hs+"-"+mm
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
        pair = hsT+"-"+mmT
        dataset3List.append(pair)
    return(dataset3List)

####################
# parsear saida Tabela final 
#app1

#def parseFinalTable(FinalTableFileName):
    #FinalTable = open(FinalTableFileName).read().split("\n")[:-1]

folderFile = sys.argv[1]
folder,filename = folderFile.split('/')

out = subprocess.Popen(['ls', folder],
      stdout = subprocess.PIPE,
      stderr = subprocess.STDOUT, encoding='utf-8')
stdout2,stderr2 = out.communicate()
fileList = stdout2.split()

folderFileapp2 = sys.argv[2]
folder2,filename2 = folderFileapp2.split('/')
out2 = subprocess.Popen(['ls', folder2],
            stdout = subprocess.PIPE,
            stderr = subprocess.STDOUT, encoding='utf-8')
stdout3,stderr3 = out2.communicate()
fileListapp2 = stdout3.split()

TableOrtho = open("/home/flavia.freitas/disk1/flavia.freitas/data/ortho_ensembl_101/Table_ID_Hsap_Mmus_Ortho.tsv").read().split("\n")[1:-1]

OrthoFinal={}
app2D = {}
for line in TableOrtho:
    col = line.split()
    #  considerar apenas ortologo 1:1
    hsapID = col[0]
    mmusID = col[1]
    #    pair = hsapID+"-"+mmusID
    if "," not in mmusID:
        pair = hsapID+"-"+mmusID
        OrthoFinal[pair] = {"ccds_prot":[],"ccds_method":[],
                            "refseq_prot":[],"refseq_method":[],
                            "spliceprot_prot":[], "spliceprot_method":[], 
                            "ensembl_prot":[],"ensembl_method":[]}
        app2D[pair] = {"ccds":[], "refseq":[],"spliceprot":[],"ensembl":[]}

for FinalFileName in fileList:
    FinalTab = open(folder+"/"+FinalFileName).read().split("\n")[:-1]
    Header1 = FinalTab.pop(0)
    db,lixo = FinalFileName.split("_",1)

    for line in FinalTab:
        col = line.split("\t")
        Hsgene = col[0].strip()
        Mmgene = col[1].strip()
        ndHsTr = col[2].strip()
        ndMmTr = col[3].strip()
        rbHsTr = col[11].strip()
        rbMmTr = col[12].strip()
        method = col[-1].strip()
        genePair = Hsgene+"-"+Mmgene
        
        if ndHsTr == "NA":
            pair = rbHsTr+"-"+rbMmTr
        if rbHsTr == "NA":
            pair = ndHsTr+"-"+ndMmTr

        if "ccds" in db:
            if genePair in OrthoFinal.keys():
                OrthoFinal[genePair]["ccds_prot"].append(pair)
                OrthoFinal[genePair]["ccds_method"].append(method)

        if "spliceprot" in db:
            if genePair in OrthoFinal.keys():
                OrthoFinal[genePair]["spliceprot_prot"].append(pair)
                OrthoFinal[genePair]["spliceprot_method"].append(method)

        if "refseq" in db:
            if genePair in OrthoFinal.keys():
                OrthoFinal[genePair]["refseq_prot"].append(pair)
                OrthoFinal[genePair]["refseq_method"].append(method)

        if "ensembl" in db:
            if genePair in OrthoFinal.keys():
                OrthoFinal[genePair]["ensembl_prot"].append(pair)
                OrthoFinal[genePair]["ensembl_method"].append(method)
#app2D

for appFileName in fileListapp2:
    app2Tab = open(folder+"/"+FinalFileName).read().split("\n")[:-1]
    Header1 = app2Tab.pop(0)
    db,lixo = appFileName.split("_",1)
    col = line.split("\t")
    Hsgene = col[0].strip()
    Mmgene = col[1].strip()
    HsTr = col[11].strip()
    MmTr = col[12].strip()
    genePair = Hsgene+"-"+Mmgene
    pair = HsTr+"-"+MmTr

    if "ccds" in db:
        if genePair in app2D.keys():
            app2D[genePair]["ccds"].append(pair)

    if "spliceprot" in db:
        if genePair in app2D.keys():
            app2D[genePair]["spliceprot"].append(pair)

    if "refseq" in db:
        if genePair in app2D.keys():
            app2D[genePair]["refseq"].append(pair)

    if "ensembl" in db:
        if genePair in app2D.keys():
            app2D[genePair]["ensembl"].append(pair)

fileDs2 = "/home/flavia.freitas/disk1/flavia.freitas/results/validation/dataset2_blanquart_ccds_ortho_pairs.txt"
fileDs4 = "/home/flavia.freitas/disk1/flavia.freitas/results/validation/dataset4_ensembl_ortho_pairs.txt"
fileDs3 = "/home/flavia.freitas/disk1/flavia.freitas/results/validation/dataset3_refseq_zambelli_ortho_pairs.txt"

ccds_list = []
refseq_list = []
ensembl_list = []
spliceprot_list = []
for gene,info in OrthoFinal.items():
    for entry in info["ccds_prot"]:
        if entry != []:
            ccds_list.append(entry) #info["ccds_prot"])
    for refseq in info["refseq_prot"]:        
        if refseq != []:
            refseq_list.append(refseq)
    for ensembl in info["ensembl_prot"]:
        if ensembl != []:
            ensembl_list.append(ensembl)
    for spliceprot in info["spliceprot_prot"]:
        if spliceprot != []:
            spliceprot_list.append(spliceprot)

app2ccds_list = []
app2refseq_list = []
app2ensembl_list = []
app2spliceprot_list = []
for gene2,info2 in app2D.items():
    app2ccds_list.append(info2["ccds"])
    app2refseq_list.append(info2["refseq"])
    app2ensembl_list.append(info2["ensembl"])
    app2spliceprot_list.append(info2["spliceprot"])

dataset2 = compareDataset2(fileDs2)
dataset3 = compareDataset3(fileDs3)
dataset4 = compareDataset4(fileDs4)
#print(refseq_list)
#print(">",dataset3)

ds2Intersection = sorted(set(ccds_list).intersection(dataset2))
ds3Intersection = sorted(set(refseq_list).intersection(dataset3))
ds4Intersection = sorted(set(ensembl_list).intersection(dataset4))
print(app2ccds_list)
print(">>","db_list","dataset","inter")
print("d2",len(ccds_list),len(dataset2),len(ds2Intersection))
print("d3",len(refseq_list),len(dataset3),len(ds3Intersection))
print("d4",len(ensembl_list),len(dataset4),len(ds4Intersection))



