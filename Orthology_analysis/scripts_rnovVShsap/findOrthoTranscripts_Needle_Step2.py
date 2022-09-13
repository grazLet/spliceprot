import sys
import os
from Bio.Emboss.Applications import NeedleCommandline
from Bio import AlignIO
import subprocess
import pickle
import pprint
import shutil
import itertools

##### ETAPA 3 #####
# apontar diretorio onde os arquivos FASTA estao armazenados
# arquivos fasta q contem as proteoformas variantes de cada gene
# para as especies de interesse 
# objetivo da analise: para cada gene, comparar as proteoformas de cada especie
# com o objetivo de encontrar proteoformas ortologas com base na similaridade de sequencia
# utiliza o needle para fazer alinhamentp
# gera diretorios para amazenamento dos resultados e arquivos temp
# sendo o principal diretorio "align_results/"
# com os resultados: OrthoAnalysis_alignment.txt e OrthoAnalysis_table.txt

hs_list=[]
rn_list=[]

hs_info=[]
rn_info=[]

# definir funcao q abre multifasta de cada gene
# com varios transcritos e cria dicionario

def multiFasta_printDict_function(multiFasta,dbName):
    hsDict={}
    rnDict={}
    hsFastaLista = []
    rnFastaLista = []
    
    for fasta in multiFasta:
        header,seq = fasta.split("\n",1)
        seq = seq.replace("\n","")
        header_itens = header.split("|")
        sp = header_itens.pop(0)
        
        if "spliceprot" in dbName:
            transEnsL = header_itens.pop(-1)
            transcriptID = "-".join(header_itens)

        if "ccds" in dbName:
            transcriptID = "-".join(header_itens)#header_itens[1]+"_"+header_itens[0]

        if "refseq" in dbName:
            transcriptID = "-".join(header_itens)#header_itens[0]+"_"+header_itens[1]

        if "ensembl" in dbName:    
            transcriptID = "-".join(header_itens)#header_itens[0]+"_"+header_itens[2]

        fasta2 = ">"+transcriptID+"\n"+seq+"\n"
        
        if "Hsap" in sp:
            hsFastaLista.append(fasta2)
        if "Rnov" in sp:
            rnFastaLista.append(fasta2)

    return(hsFastaLista,rnFastaLista)




# inicialmente selecionar apenas seqs de hsap e rnov
# rodar o needle para cada combinacao
# depois ampliar para rnov e automatizar para todos os genes

# passar diretorio onde estao os dados q devem ser analisados
# dados = arquivos fasta com os transcritos dos genes ortologos
# das duas ou tres especies analisadas
# rodar ls comando
folder = sys.argv[1]
folder2 = folder.replace("HsapRnov_OrthoProt/","")
folder_itens = folder.split("/")
database = folder_itens[0]
out = subprocess.Popen(['ls','-Sr', folder], 
           stdout = subprocess.PIPE, 
           stderr = subprocess.STDOUT, 
           encoding='utf-8')
stdout2,stderr2 = out.communicate()

fileList = stdout2.split()

#os.mkdir(folder2+"needle_files")
#os.mkdir(folder2+"temp_fasta")
#os.mkdir(folder2+"align_results")
#os.mkdir(folder2+"temp_results")

# guardar arquivos com seqs de proteoformas
# nao-ortologas
#os.mkdir("non_orthologs_fasta")

fileout4 = open(folder2+"align_results/OrthoAnalysis_table.txt","w")

alignheader = "hsTranscriptID"+"\t"+"rnTranscriptID"+"\t"+"similarity"+"\t"+"similarityRatio"+"\t"+"identity"+"\t"+"identityRatio"+"\t"+"gap"+"\t"+"gapRatio"+"\t"+"score"+"\n"

fileout4.write(alignheader)
#os.mkdir(folder2+"Needle_db")
fileNeedleOut = open(folder2+"Needle_db/Needle_complete_db_table.txt","w")
fileNeedleOut.write(alignheader)
fileAlignNeedleOut = open(folder2+"Needle_db/Needle_complete_db_align.txt","w")
fileout5 = open(folder2+"align_results/OrthoAnalysis_alignment.txt","w")

#headerTable = []
# cada arquivo um gene
# com todas as proteoformas

ProtOrt = {}
count  = 0
count2 = 0
count3 = 0
for files2 in fileList:
        #count3 += 1
        #msg3 = "File Number " + str(count3)
        #print(msg3)
        files = folder+files2
        multiFasta2 = open(files).read().split(">")[1:]
        multiFasta3 = multiFasta_printDict_function(multiFasta2,database)
        hsL = multiFasta3[0]
        rnL = multiFasta3[1]

    # criar lista para guardar ID de Rnov 
    # observar se um mesmo Rnov Transcript ID
    # sera associado a mais de um Hsap Transcript ID 
        RnovAlignD = {}
    # abrir arquivo para guardar todos os resultados do needle
    # independentemente se passa ou nao pelos citerios 0.8 e 0.01
    # arquivo para consulta
        F_TempResultName = folder2+"temp_results/"+files2.replace("ortho_proteins.fasta","")+"align_results.txt" 
        F_TempResult = open(F_TempResultName,"w")
    # key hs1 -- transcript ID hs
    # value hs2["fasta"] e hs2["geneID"]
        hsAlignDict = {}
        AlignSeqDict = {}
        for hsFasta2 in hsL:
                hs_header2,hs_seq2 = hsFasta2.split("\n",1)
                hs1_ID_1 = hs_header2.replace(">","").split("-")
                hs1_ID = "hs_"+hs1_ID_1[-1]
                # criar arquivo1 q sera lido no alinhamento
                if len(hs1_ID) >= 50:
                    #hs1_original = hs1
                        x = 50 - len(hs1_ID)
                        hs1_new = hs1_ID[:x]
                        hs1_ID = hs1_new
                if len(hs1_ID) <  50 :
                        hs1_ID = hs1_ID
                hsFilename2 = folder2+"temp_fasta/"+hs1_ID+"_temp.fa"
                hsfile2 = open(hsFilename2,"w")
                hsfile2.write(hsFasta2)
                hsfile2.close()

        for rnFasta2 in rnL:
            rn_header2,rn_seq2 = rnFasta2.split("\n",1)
            rn1_ID_1 = rn_header2.replace(">","").split("-")
            rn1_ID = "rn_"+rn1_ID_1[-1]

            if len(rn1_ID) >= 50:
                xx = 50 - len(rn1_ID)
                rn1_new = rn1_ID[:xx]
                rn1_ID = rn1_new
            if len(rn1_ID) < 50:
                rn1_ID = rn1_ID

            # criar arquivo1 q sera lido no alinhamento
            rnFilename2 = folder2+"temp_fasta/"+rn1_ID+"_temp.fa"
            rnfile2 = open(rnFilename2,"w")
            rnfile2.write(rnFasta2)
            rnfile2.close()


        for hsFASTA,rnFASTA in itertools.product(hsL,rnL):
                hs_header2,hs_seq2 = hsFASTA.split("\n",1)
                hs1_1 = hs_header2.replace(">","").split("-")
                hs1 = "hs_"+hs1_1[-1]
                rn_header2,rn_seq2 = rnFASTA.split("\n",1)
                rn1_1 = rn_header2.replace(">","").split("-")
                rn1 = "rn_"+rn1_1[-1]
                #rnFilename = folder2+"temp_fasta/"+rn1+"_temp.fa"
            
        # criar nome do arquivo q vai guardar resultado do alinhament
                if len(hs1) >= 50:
                    x = 50 - len(hs1)
                    hs2_new = hs1[:x]
                    hsFilename = folder2+"temp_fasta/"+hs2_new+"_temp.fa"
                        #outFname = folder2+"needle_files/"+hs2_new+"_"+rn1+"_needle.txt"
                if len(hs1) < 50:
                    hs2_new = hs1
                    hsFilename = folder2+"temp_fasta/"+hs2_new+"_temp.fa"
                        #outFname = folder2+"needle_files/"+hs2_new+"_"+rn1+"_needle.txt"
                if len(rn1) >= 50:
                    xx= 50 - len(rn1)
                    rn2_new = rn1[:xx]
                    rnFilename = folder2+"temp_fasta/"+rn2_new+"_temp.fa"
                if len(rn1) < 50:
                    rn2_new = rn1
                    rnFilename = folder2+"temp_fasta/"+rn2_new+"_temp.fa"
                outFname = folder2+"needle_files/"+hs2_new+"_"+rn2_new+"_needle.txt"


        # rodar o alinhamento needle
                needle_cline = NeedleCommandline(asequence=hsFilename, bsequence=rnFilename, gapopen=10, gapextend=0.5, outfile=outFname)
                stdout,stderr = needle_cline()
        
        # parse o resultado do alinhamento
        # guardado no texto criado na variavel outFname
                itens2 = open(outFname).read()
                fileAlignNeedleOut.write(itens2)
                fileAlignNeedleOut.write("\n")
                itens = itens2.split("#=======================================")
                align_info = itens[1].split("\n")
                seqID_1 = align_info[3].replace("# 1:","").strip()
                seqID_2 = align_info[4].replace("# 2:","").strip()
                lenght = align_info[9].replace("# Length:","").strip()
                identity_info = align_info[10].replace("# Identity:","").split()
                identity = identity_info[0].strip()
                ident1,ident2 = identity.split("/")
                ident_ratio = format((int(ident1)/int(ident2)),".3f")
                identity_perc = identity_info[1].replace("(","").replace(")","").strip()
                simil_info = align_info[11].replace("# Similarity:","").split()
                similarity = simil_info[0].strip()
                sim1,sim2 = similarity.split("/")
                sim_ratio = format((int(sim1)/int(sim2)),".3f")
                similarity_perc = simil_info[1].replace("(","").replace(")","").strip()
                gaps_info = align_info[12].replace("# Gaps:","").split()
                gaps = gaps_info[0].strip()
                gap1,gap2 = gaps.split("/")
                gap_ratio = format((int(gap1)/int(gap2)),".3f")
                gaps_perc = gaps_info[1].replace("(","").replace(")","").strip()
                score = float(align_info[13].replace("# Score:","").strip())

        #isolar alinhamento para conferir depois
                alignment_itens = itens[2].split("\n")[2:-5]
                alignment_seq = "\n".join(alignment_itens)
        # usar Hs transcriptID como key
        # criar lista com Rn transcript ID, similarity, identity, score, gaps
        # preencher dict hsAlignDict com os resultados do needle
                os.remove(outFname)
                hs1_original = hs_header2.replace(">","")
                rn1_original = rn_header2.replace(">","")
                if hs1_original in hsAlignDict.keys():
                        rn1List = []
                        rn1List = hsAlignDict[hs1_original]["rn1"]
                        rn1List.append(rn1_original)
                        hsAlignDict[hs1_original]["rn1"] = rn1List

                        identRatioList = []
                        identRatioList = hsAlignDict[hs1_original]["ident_ratio"]
                        identRatioList.append(ident_ratio)
                        hsAlignDict[hs1_original]["ident_ratio"] = identRatioList

                        simRatioList = []
                        simRatioList = hsAlignDict[hs1_original]["sim_ratio"]
                        simRatioList.append(sim_ratio)
                        hsAlignDict[hs1_original]["sim_ratio"] = simRatioList

                        scoreList = []
                        scoreList = hsAlignDict[hs1_original]["score"]
                        scoreList.append(score)
                        hsAlignDict[hs1_original]["score"] = scoreList

                        gapRatioList = []
                        gapRatioList = hsAlignDict[hs1_original]["gap_ratio"]
                        gapRatioList.append(gap_ratio)
                        hsAlignDict[hs1_original]["gap_ratio"] = gapRatioList
            
                        alignSeqList = []
                        alignSeqList = hsAlignDict[hs1_original]["align_seq"]
                        alignSeqList.append(alignment_seq)
                        hsAlignDict[hs1_original]["align_seq"] = alignSeqList

                        similarityList = []
                        similarityList = hsAlignDict[hs1_original]["similarity"]
                        similarityList.append(similarity)
                        hsAlignDict[hs1_original]["similarity"] = similarityList

                        identityList = []
                        identityList = hsAlignDict[hs1_original]["identity"]
                        identityList.append(identity)
                        hsAlignDict[hs1_original]["identity"] = identityList

                        gapsList = []
                        gapsList = hsAlignDict[hs1_original]["gap"]
                        gapsList.append(gaps)
                        hsAlignDict[hs1_original]["gap"] = gapsList

                if hs1_original not in hsAlignDict.keys():
                        hsAlignDict[hs1_original] = {"rn1":[rn1_original],
                    "ident_ratio":[ident_ratio],
                    "sim_ratio":[sim_ratio],
                    "score":[score],
                    "gap_ratio":[gap_ratio],
                    "align_seq":[alignment_seq],
                    "similarity":[similarity],
                    "identity":[identity],
                    "gap":[gaps]}

        # para cada transcrito hs
        # filtrar melhor proteoforma de Rn
         
        for hsT1,AlignInfo in hsAlignDict.items():
            
        # hsT1 eh o transcriptID Hs
        # recuperar valores de similarity na lista SimValueList
        # identificar index do valor mais alto na lista
        # usar o index para recuperar ID dos transcritos Hs e Rn
        # e alinhamento
        # repetir procedimento para score e ident
                for aa in range(len(AlignInfo["rn1"])):
                        
                        alignIDinfoLista2 = [hsT1,AlignInfo["rn1"][aa],AlignInfo["similarity"][aa], \
                            AlignInfo["sim_ratio"][aa],AlignInfo["identity"][aa],AlignInfo["ident_ratio"][aa],\
                            AlignInfo["gap"][aa],AlignInfo["gap_ratio"][aa],AlignInfo["score"][aa]]
                        alignIDinfo2 = "\t".join(map(str,alignIDinfoLista2))

                        F_TempResult.write(alignIDinfo2)
                        F_TempResult.write("\n")
                        fileNeedleOut.write(alignIDinfo2)
                        fileNeedleOut.write("\n")
                #count += 1
                #msg = "Hsap sequence number " + str(count)
                #print(msg)
                SimValueList = AlignInfo["sim_ratio"]
                ScoreValueList = AlignInfo["score"]
                GapValueList = AlignInfo["gap_ratio"]
                indexlist = []
                #GapValueList2 = []
                #ScoreValueList2 = []
                # ETAPA filtragem
                for simIndex,simvalue in enumerate(SimValueList):
                        
                        if float(simvalue) >= 0.8 and \
                                float(AlignInfo["gap_ratio"][simIndex]) <= 0.01:
                                    if SimValueList.count(max(SimValueList)) > 1 :               
                                            if float(simvalue) == float(max(SimValueList)):
                                                    indexlist.append(simIndex)
                
                                    if SimValueList.count(max(SimValueList)) == 1 :
                                            if float(simvalue) == float(max(SimValueList)):
                                                    indexlist = [simIndex]
                ScoreSubconj = []
                GapSubconj = []
                finalIndex = ""
                if len(indexlist) == 1:
                        finalIndex = indexlist[0]
                       # print(">>>",finalIndex)
                if len(indexlist) > 1:
                        for indexSim in indexlist:
                                #seleciona valores de score com mesmo index dos maiores valores de similaridade
                                ScoreSubconj.append(ScoreValueList[indexSim])
                                GapSubconj.append(GapValueList[indexSim])
                                # identifica o maior valor de score dentre os pares com mesmo valor de similaridade
                                # identifica o index do maior valor de score
                                # para usa-lo para recuperar info e preencher o dicionario
                        maxScoreSubconj = max(ScoreSubconj)
                        #indexMaxScoreTotal = ScoreValueList.index(maxScoreSubconj)
                        minGapSubconj = min(GapSubconj)
                        #indexMinGapTotal = GapValueList.index(minGapSubconj)
                        if GapSubconj.count(minGapSubconj) == 1:
                            finalIndex = GapValueList.index(minGapSubconj)
                        if GapSubconj.count(minGapSubconj) > 1:
                            finalIndex = ScoreValueList.index(maxScoreSubconj)
                
                    # usar um segundo dicitonario para evitar 
                    # q uma mesma proteoforma de Rnov seja eleita 
                    # como ortologa de mais de uma proteoforma de Hsap                 

                if type(finalIndex) != str:
                        if AlignInfo["rn1"][finalIndex] in RnovAlignD.keys():
                                rn1List2 = []
                                rn1List2 = RnovAlignD[AlignInfo["rn1"][finalIndex]]["hsT1"]
                                rn1List2.append(hsT1)
                                RnovAlignD[AlignInfo["rn1"][finalIndex]]["hsT1"] = rn1List2
            
                                identRatioList2 = []
                                identRatioList2 = RnovAlignD[AlignInfo["rn1"][finalIndex]]["ident_ratio"]
                                identRatioList2.append(float(AlignInfo["ident_ratio"][finalIndex]))                            
                                RnovAlignD[AlignInfo["rn1"][finalIndex]]["ident_ratio"] = identRatioList2
    
                                simRatioList2 = []
                                simRatioList2 = RnovAlignD[AlignInfo["rn1"][finalIndex]]["sim_ratio"]
                                simRatioList2.append(float(AlignInfo["sim_ratio"][finalIndex]))
                                RnovAlignD[AlignInfo["rn1"][finalIndex]]["sim_ratio"] = simRatioList2
    
                                scoreList2 = []
                                scoreList2 = RnovAlignD[AlignInfo["rn1"][finalIndex]]["score"]
                                scoreList2.append(AlignInfo["score"][finalIndex])
                                RnovAlignD[AlignInfo["rn1"][finalIndex]]["score"] = scoreList2

                                gapRatioList2 = []
                                gapRatioList2 = RnovAlignD[AlignInfo["rn1"][finalIndex]]["gap_ratio"]
                                gapRatioList2.append(float(AlignInfo["gap_ratio"][finalIndex]))
                                RnovAlignD[AlignInfo["rn1"][finalIndex]]["gap_ratio"] = gapRatioList2
                                
                                alignSeqList2 = []
                                alignSeqList2 = RnovAlignD[AlignInfo["rn1"][finalIndex]]["align_seq"]
                                alignSeqList2.append(AlignInfo["align_seq"][finalIndex])
                                RnovAlignD[AlignInfo["rn1"][finalIndex]]["align_seq"] = alignSeqList2
                
                                similarityList2 = []
                                similarityList2 = RnovAlignD[AlignInfo["rn1"][finalIndex]]["similarity"]
                                similarityList2.append(AlignInfo["similarity"][finalIndex])
                                RnovAlignD[AlignInfo["rn1"][finalIndex]]["similarity"] = similarityList2

                                identityList2 = []
                                identityList2 = RnovAlignD[AlignInfo["rn1"][finalIndex]]["identity"]
                                identityList2.append(AlignInfo["identity"][finalIndex])
                                RnovAlignD[AlignInfo["rn1"][finalIndex]]["identity"] = identityList2

                                gapsList2 = []
                                gapsList2 = RnovAlignD[AlignInfo["rn1"][finalIndex]]["gap"]
                                gapsList2.append(AlignInfo["gap"][finalIndex])
                                RnovAlignD[AlignInfo["rn1"][finalIndex]]["gap"] = gapsList2
                
                        if AlignInfo["rn1"][finalIndex] not in RnovAlignD.keys():
                                RnovAlignD[AlignInfo["rn1"][finalIndex]] = {"hsT1": [hsT1],
                                                                       "similarity":[AlignInfo["similarity"][finalIndex]],
                                                                       "sim_ratio":[float(AlignInfo["sim_ratio"][finalIndex])],
                                                                       "identity":[AlignInfo["identity"][finalIndex]],
                                                                       "ident_ratio":[float(AlignInfo["ident_ratio"][finalIndex])],
                                                                       "gap":[AlignInfo["gap"][finalIndex]],
                                                                       "gap_ratio":[float(AlignInfo["gap_ratio"][finalIndex])],
                                                                        "align_seq":[AlignInfo["align_seq"][finalIndex]],
                                                                        "score":[AlignInfo["score"][finalIndex]]}

        removePath1 = folder2+"temp_fasta/*"
        removePath2 = " ".join(["rm",removePath1])
        #os.system("rm temp_fasta/*")
        os.system(removePath2)

        for rnT1,RnovAlignInfo in RnovAlignD.items():
                SimValueList3 = RnovAlignInfo["sim_ratio"]
                ScoreValueList3 = RnovAlignInfo["score"]
                GapValueList3 = RnovAlignInfo["gap_ratio"]
                indexListRn = []
                       
                for SimIndex3,SimValue3 in enumerate(SimValueList3):
                        if float(SimValue3) >= 0.8 and \
                            float(RnovAlignInfo["gap_ratio"][SimIndex3]) <= 0.01:
                                if SimValueList3.count(max(SimValueList3)) > 1:
                                        if float(SimValue3) == float(max(SimValueList3)):
                                                indexListRn.append(SimIndex3) 

                                if SimValueList3.count(max(SimValueList3)) == 1:
                                        indexListRn = [SimIndex3]
                # alimentar dict com pares q passararam no filtroo
                ScoreSubconjRn = []
                GapSubconjRn = []
                finalIndexRn = ""
                
                if len(indexListRn) == 1:
                        finalIndexRn = indexListRn[0]
                
                if len(indexListRn) > 1:
                        for indexSimRn in indexListRn:
                                ScoreSubconjRn.append(ScoreValueList3[indexSimRn])
                                GapSubconjRn.append(GapValueList3[indexSimRn])
                        maxScoreSubconjRn = max(ScoreSubconjRn)
                        minGapSubconjRn = min(GapSubconjRn) 
                
                        if GapSubconjRn.count(minGapSubconjRn) == 1:
                                finalIndexRn = GapValueList3.index(minGapSubconjRn)

                        if GapSubconjRn.count(minGapSubconjRn) > 1:
                                finalIndexRn = ScoreValueList3.index(maxScoreSubconjRn)

                if type(finalIndexRn) != str:
                        hsTransId = RnovAlignInfo["hsT1"][finalIndexRn]

                        if hsTransId in ProtOrt.keys():
                                rnT1Lista = []
                                rnT1Lista = ProtOrt[hsTransId]
                                rnT1Lista.append(rnT1)
                                ProtOrt[hsTransId] = rnT1Lista

                        if hsTransId not in ProtOrt.keys():
        
                                ProtOrt[hsTransId] = [rnT1]
                                alignIDinfoLista2 = [RnovAlignInfo["hsT1"][finalIndexRn],rnT1,RnovAlignInfo["similarity"][finalIndexRn], \
                                    RnovAlignInfo["sim_ratio"][finalIndexRn],RnovAlignInfo["identity"][finalIndexRn], \
                                    RnovAlignInfo["ident_ratio"][finalIndexRn],RnovAlignInfo["gap"][finalIndexRn], \
                                    RnovAlignInfo["gap_ratio"][finalIndexRn],RnovAlignInfo["score"][finalIndexRn]]

                                alignIDinfo2 = "\t".join(map(str,alignIDinfoLista2))
                                fileout4.write(alignIDinfo2)
                                fileout4.write("\n")
                                # gerar arquivo com alinhamento
                                fileout5.write(alignheader)
                                fileout5.write(alignIDinfo2)
                                fileout5.write("\n")
                                fileout5.write("\n")
                                fileout5.write(RnovAlignInfo["align_seq"][finalIndexRn])
                                fileout5.write("\n")
                                fileout5.write("\n")
                                fileout5.write("#------------------------------------------------\n")
                                fileout5.write("\n")
                                fileout5.write("\n")
                
        F_TempResult.close()
#print ("##############################################")
#print("remocao diretorios temp")
#print("\n\n")
os.rmdir(folder2+"temp_fasta/")
os.rmdir(folder2+"needle_files/") 

fileNeedleOut.close()
fileAlignNeedleOut.close()
fileout4.close()
fileout5.close()
#print ("##############################################")
#print("finalizacao...")
#print("\n\n")
# salvar dicionario com pares de proteoformas
#print ("##############################################")
#print("criando dicionario... Etapa Final")
#print("\n\n")

#os.mkdir(folder2+"ortho-pairs_app1")

pickleFileName = folder2+"ortho-pairs_app1/HsRn-ProteoOrtho_v1.pickle"
pickleFileout =  open(pickleFileName, "wb")
pickle.dump(ProtOrt,pickleFileout)
pickleFileout.close()

ProtOrtD = open(folder2+"ortho-pairs_app1/HsRn-ProteoOrtho_v1.txt","w")
ProtValues = list(ProtOrt.values())
ProtValuesCount = []
for element in ProtValues:
    if type(element) is not list:
        ProtValuesCount.append(element)
    if type(element) is list:
        for index3 in range(len(element)):
            ProtValuesCount.append(element[index3])

for k,v in ProtOrt.items():
    if type(v) is not list:
        printelement = k+"\t"+v+"\t"+str(ProtValuesCount.count(v))+"\n"
        ProtOrtD.write(printelement)
    if type(v) is list:
        for index4 in range(len(v)):
            printelement = k+"\t"+v[index4]+"\t"+str(ProtValuesCount.count(v[index4]))+"\n"
            ProtOrtD.write(printelement)

ProtOrtD.close()
#pprint.pprint(ProtOrt) 

