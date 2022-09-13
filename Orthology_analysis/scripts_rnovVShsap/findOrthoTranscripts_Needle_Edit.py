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
mm_list=[]

hs_info=[]
mm_info=[]

# definir funcao q abre multifasta de cada gene
# com varios transcritos e cria dicionario

def multiFasta_printDict_function(multiFasta,dbName):
    hsDict={}
    mmDict={}
    hsFastaLista = []
    mmFastaLista = []
    
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
        if "Mmus" in sp:
            mmFastaLista.append(fasta2)

    return(hsFastaLista,mmFastaLista)




# inicialmente selecionar apenas seqs de hsap e mmus
# rodar o needle para cada combinacao
# depois ampliar para rnov e automatizar para todos os genes

# passar diretorio onde estao os dados q devem ser analisados
# dados = arquivos fasta com os transcritos dos genes ortologos
# das duas ou tres especies analisadas
# rodar ls comando
folder = sys.argv[1]
folder2 = folder.replace("HsapMmus_OrthoProt/","")
folder_itens = folder.split("/")
database = folder_itens[0]
out = subprocess.Popen(['ls','-Sr', folder], 
           stdout = subprocess.PIPE, 
           stderr = subprocess.STDOUT, encoding='utf-8')
stdout2,stderr2 = out.communicate()

fileList = stdout2.split()

os.mkdir(folder2+"needle_files")
os.mkdir(folder2+"temp_fasta")
os.mkdir(folder2+"align_results")
os.mkdir(folder2+"temp_results")

# guardar arquivos com seqs de proteoformas
# nao-ortologas
#os.mkdir("non_orthologs_fasta")

fileout4 = open(folder2+"align_results/OrthoAnalysis_table.txt","w")

alignheader = "hsTranscriptID"+"\t"+"mmTranscriptID"+"\t"+"similarity"+"\t"+"similarityRatio"+"\t"+"identity"+"\t"+"identityRatio"+"\t"+"gap"+"\t"+"gapRatio"+"\t"+"score"+"\n"

fileout4.write(alignheader)
os.mkdir(folder2+"Needle_db")
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
        mmL = multiFasta3[1]

    # criar lista para guardar ID de Mmus 
    # observar se um mesmo Mmus Transcript ID
    # sera associado a mais de um Hsap Transcript ID 
        MmusAlignD = {}
    # abrir arquivo para guardar todos os resultados do needle
    # independentemente se passa ou não pelos citerios 0.8 e 0.01
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

        for mmFasta2 in mmL:
            mm_header2,mm_seq2 = mmFasta2.split("\n",1)
            mm1_ID_1 = mm_header2.replace(">","").split("-")
            mm1_ID = "mm_"+mm1_ID_1[-1]

            if len(mm1_ID) >= 50:
                xx = 50 - len(mm1_ID)
                mm1_new = mm1_ID[:xx]
                mm1_ID = mm1_new
            if len(mm1_ID) < 50:
                mm1_ID = mm1_ID

            # criar arquivo1 q sera lido no alinhamento
            mmFilename2 = folder2+"temp_fasta/"+mm1_ID+"_temp.fa"
            mmfile2 = open(mmFilename2,"w")
            mmfile2.write(mmFasta2)
            mmfile2.close()


        for hsFASTA,mmFASTA in itertools.product(hsL,mmL):
                hs_header2,hs_seq2 = hsFASTA.split("\n",1)
                hs1_1 = hs_header2.replace(">","").split("-")
                hs1 = "hs_"+hs1_1[-1]
                mm_header2,mm_seq2 = mmFASTA.split("\n",1)
                mm1_1 = mm_header2.replace(">","").split("-")
                mm1 = "mm_"+mm1_1[-1]
                #mmFilename = folder2+"temp_fasta/"+mm1+"_temp.fa"
            
        # criar nome do arquivo q vai guardar resultado do alinhament
                if len(hs1) >= 50:
                    x = 50 - len(hs1)
                    hs2_new = hs1[:x]
                    hsFilename = folder2+"temp_fasta/"+hs2_new+"_temp.fa"
                        #outFname = folder2+"needle_files/"+hs2_new+"_"+mm1+"_needle.txt"
                if len(hs1) < 50:
                    hs2_new = hs1
                    hsFilename = folder2+"temp_fasta/"+hs2_new+"_temp.fa"
                        #outFname = folder2+"needle_files/"+hs2_new+"_"+mm1+"_needle.txt"
                if len(mm1) >= 50:
                    xx= 50 - len(mm1)
                    mm2_new = mm1[:xx]
                    mmFilename = folder2+"temp_fasta/"+mm2_new+"_temp.fa"
                if len(mm1) < 50:
                    mm2_new = mm1
                    mmFilename = folder2+"temp_fasta/"+mm2_new+"_temp.fa"
                outFname = folder2+"needle_files/"+hs2_new+"_"+mm2_new+"_needle.txt"


        # rodar o alinhamento needle
                needle_cline = NeedleCommandline(asequence=hsFilename, bsequence=mmFilename, gapopen=10, gapextend=0.5, outfile=outFname)
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
        # criar lista com Mm transcript ID, similarity, identity, score, gaps
        # preencher dict hsAlignDict com os resultados do needle
                os.remove(outFname)
                hs1_original = hs_header2.replace(">","")
                mm1_original = mm_header2.replace(">","")
                if hs1_original in hsAlignDict.keys():
                                       
                        hsAlignDict[hs1_original]["mm1"].append(mm1_original)
                        hsAlignDict[hs1_original]["ident_ratio"].append(ident_ratio)
                        hsAlignDict[hs1_original]["sim_ratio"].append(sim_ratio)
                        hsAlignDict[hs1_original]["score"].append(score)
                        hsAlignDict[hs1_original]["gap_ratio"].append(gap_ratio)
                        hsAlignDict[hs1_original]["align_seq"].append(alignment_seq)
                        hsAlignDict[hs1_original]["similarity"].append(similarity)
                        hsAlignDict[hs1_original]["identity"].append(identity)
                        hsAlignDict[hs1_original]["gap"].append(gaps)

                if hs1_original not in hsAlignDict.keys():
                        hsAlignDict[hs1_original] = {"mm1":[mm1_original],
                    "ident_ratio":[ident_ratio],
                    "sim_ratio":[sim_ratio],
                    "score":[score],
                    "gap_ratio":[gap_ratio],
                    "align_seq":[alignment_seq],
                    "similarity":[similarity],
                    "identity":[identity],
                    "gap":[gaps]}

        # para cada transcrito hs
        # filtrar melhor proteoforma de Mm
         
        for hsT1,AlignInfo in hsAlignDict.items():
            
        # hsT1 eh o transcriptID Hs
        # recuperar valores de similarity na lista SimValueList
        # identificar index do valor mais alto na lista
        # usar o index para recuperar ID dos transcritos Hs e Mm
        # e alinhamento
        # repetir procedimento para score e ident
                for aa in range(len(AlignInfo["mm1"])):
                        
                        alignIDinfoLista2 = [hsT1,AlignInfo["mm1"][aa],AlignInfo["similarity"][aa], \
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
                    # q uma mesma proteoforma de Mmus seja eleita 
                    # como ortologa de mais de uma proteoforma de Hsap                 

                if type(finalIndex) != str:
                        if AlignInfo["mm1"][finalIndex] in MmusAlignD.keys():
                                MmusAlignD[AlignInfo["mm1"][finalIndex]]["hsT1"].append(hsT1)
                                MmusAlignD[AlignInfo["mm1"][finalIndex]]["ident_ratio"].append(float(AlignInfo["ident_ratio"][finalIndex]))
                                MmusAlignD[AlignInfo["mm1"][finalIndex]]["sim_ratio"].append(float(AlignInfo["sim_ratio"][finalIndex]))
                                MmusAlignD[AlignInfo["mm1"][finalIndex]]["score"].append(AlignInfo["score"][finalIndex])
                                MmusAlignD[AlignInfo["mm1"][finalIndex]]["gap_ratio"].append(float(AlignInfo["gap_ratio"][finalIndex]))
                                MmusAlignD[AlignInfo["mm1"][finalIndex]]["align_seq"].append(AlignInfo["align_seq"][finalIndex])
                                MmusAlignD[AlignInfo["mm1"][finalIndex]]["similarity"].append(AlignInfo["similarity"][finalIndex])
                                MmusAlignD[AlignInfo["mm1"][finalIndex]]["identity"].append(AlignInfo["identity"][finalIndex])
                                MmusAlignD[AlignInfo["mm1"][finalIndex]]["gap"].append(AlignInfo["gap"][finalIndex])
                
                        if AlignInfo["mm1"][finalIndex] not in MmusAlignD.keys():
                                MmusAlignD[AlignInfo["mm1"][finalIndex]] = {"hsT1": [hsT1],
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

        for mmT1,MmusAlignInfo in MmusAlignD.items():
                SimValueList3 = MmusAlignInfo["sim_ratio"]
                ScoreValueList3 = MmusAlignInfo["score"]
                GapValueList3 = MmusAlignInfo["gap_ratio"]
                indexListMm = []
                       
                for SimIndex3,SimValue3 in enumerate(SimValueList3):
                        if float(SimValue3) >= 0.8 and \
                            float(MmusAlignInfo["gap_ratio"][SimIndex3]) <= 0.01:
                                if SimValueList3.count(max(SimValueList3)) > 1:
                                        if float(SimValue3) == float(max(SimValueList3)):
                                                indexListMm.append(SimIndex3) 

                                if SimValueList3.count(max(SimValueList3)) == 1:
                                        if float(SimValue3) == float(max(SimValueList3)):
                                                indexListMm = [SimIndex3]
                # alimentar dict com pares q passararam no filtroo
                ScoreSubconjMm = []
                GapSubconjMm = []
                finalIndexMm = ""
                
                if len(indexListMm) == 1:
                        finalIndexMm = indexListMm[0]
                
                if len(indexListMm) > 1:
                        for indexSimMm in indexListMm:
                                ScoreSubconjMm.append(ScoreValueList3[indexSimMm])
                                GapSubconjMm.append(GapValueList3[indexSimMm])
                        maxScoreSubconjMm = max(ScoreSubconjMm)
                        minGapSubconjMm = min(GapSubconjMm) 
                
                        if GapSubconjMm.count(minGapSubconjMm) == 1:
                                finalIndexMm = GapValueList3.index(minGapSubconjMm)

                        if GapSubconjMm.count(minGapSubconjMm) > 1:
                                finalIndexMm = ScoreValueList3.index(maxScoreSubconjMm)

                if type(finalIndexMm) != str:
                        hsTransId = MmusAlignInfo["hsT1"][finalIndexMm]

                        if hsTransId in ProtOrt.keys():
                                mmT1Lista = []
                                mmT1Lista = ProtOrt[hsTransId]
                                mmT1Lista.append(mmT1)
                                ProtOrt[hsTransId] = mmT1Lista

                        if hsTransId not in ProtOrt.keys():
        
                                ProtOrt[hsTransId] = [mmT1]
                                alignIDinfoLista2 = [MmusAlignInfo["hsT1"][finalIndexMm],mmT1,MmusAlignInfo["similarity"][finalIndexMm], \
                                    MmusAlignInfo["sim_ratio"][finalIndexMm],MmusAlignInfo["identity"][finalIndexMm], \
                                    MmusAlignInfo["ident_ratio"][finalIndexMm],MmusAlignInfo["gap"][finalIndexMm], \
                                    MmusAlignInfo["gap_ratio"][finalIndexMm],MmusAlignInfo["score"][finalIndexMm]]

                                alignIDinfo2 = "\t".join(map(str,alignIDinfoLista2))
                                fileout4.write(alignIDinfo2)
                                fileout4.write("\n")
                                # gerar arquivo com alinhamento
                                fileout5.write(alignheader)
                                fileout5.write(alignIDinfo2)
                                fileout5.write("\n")
                                fileout5.write("\n")
                                fileout5.write(MmusAlignInfo["align_seq"][finalIndexMm])
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

os.mkdir(folder2+"ortho-pairs_app1")

pickleFileName = folder2+"ortho-pairs_app1/HsMm-ProteoOrtho_v1.pickle"
pickleFileout =  open(pickleFileName, "wb")
pickle.dump(ProtOrt,pickleFileout)
pickleFileout.close()

ProtOrtD = open(folder2+"ortho-pairs_app1/HsMm-ProteoOrtho_v1.txt","w")
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

