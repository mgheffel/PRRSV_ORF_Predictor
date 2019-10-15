import os
import pandas as pd
import re
#declare path to test genomes
path=os.getcwd()+'\\testGenomes\\'
#put all test genome flat file names into an array
fastafiles=[]
for filename in os.listdir(path):
    if filename[-6:]==".fasta":
        fastafiles.append(filename)

#load scoring matrices
print("Loading in Matricies")
loadDataPath=os.getcwd()+'\\loadData\\'
orf1matrix=pd.read_pickle(loadDataPath+"orf1matrix.pckl")
allORFmatrix=pd.read_pickle(loadDataPath+"allORFmatrix.pckl")
with open(loadDataPath+"thresholdScores.txt") as f:
    orf1ScoreThreshold=float(f.readline())
    allORFScoreThreshold=float(f.readline())

def isolateGenomeToRNA(filename):
    with open(filename,'r') as f:
        fileString=f.read()
        name=fileString.split("\n")[0]
        fileString=fileString[len(name):].replace("\n","").lower()
    #ISOLATE GENOME
    genome=fileString[fileString.find("ORIGIN")+1:len(fileString)]
    genome=re.sub(r'[^a-z]+',"",genome)
    #print(genome)
    
    #CHANGE T BASE TO U TO REPRESENT AS RNA
    genome=genome.replace('t','u')
    return genome

def scoreSequence(seq,matrix):
    positions=[0,1,2,3,4,5,6,7,8,9,13,14,15,16,17,18,19,20,21,22,23,24]
    seq=seq.upper()
    score=1
    for i in positions:
        try:
            positionScore=matrix[seq[i]][i]
        except:
            positionScore=0
        if positionScore==0:
            score*=.000001
        else:
            score*=positionScore
    return score

#FUNCTION FOR FINDING THE END OF AN ORF BY LOOKING FOR THE FIRST SOPT CODON IN THE READING FRME
def findNextStop(seq,start,end):
        stopCodons=["uag","uga","uaa"]
        for i in range(start,end,3):
            if seq[i:i+3] in stopCodons:
                return i+3

def annotate(rawGenome):
    #dictionary to store open reading frames
    openReadingFrames={}
    
    #Locate ORF1a start
    for i in range(len(rawGenome)):
        if rawGenome[i:i+3]=='aug':
            try:
                seq=rawGenome[i-10:i+15]
            except:
                start1a=i
                break
            score=scoreSequence(seq,orf1matrix)
            if score>orf1ScoreThreshold:
                start1a=i
    
    #Locate ends of ORF1a/b and add orf1a/b to dictionary
    end1a=findNextStop(rawGenome,start1a,len(rawGenome))
    openReadingFrames["1a"]=str(start1a+1) + "," + str(end1a)
    end1b=findNextStop(rawGenome,end1a-1,len(rawGenome))
    openReadingFrames["1b"]=str(start1a+1) + "," + str(end1b)
    
    #Loacte ORF2a start
    for i in range(end1b,len(rawGenome)):
        if rawGenome[i:i+3]=='aug':
            seq=rawGenome[i-10:i+15]
            score=scoreSequence(seq, allORFmatrix)
            if score>allORFScoreThreshold:
                start2a=i+1
                break
    #Locate ORF2a end
    end2a=findNextStop(rawGenome,start2a-1,len(rawGenome))
    #add ORF2a to dictionary
    openReadingFrames["2a"]=str(start2a)+","+str(end2a)
    
    #Locate ORF2b start
    for i in range(start2a+1,len(rawGenome)):
        if rawGenome[i:i+3]=='aug':
            seq=rawGenome[i-10:i+15]
            score=scoreSequence(seq, allORFmatrix)
            if score>allORFScoreThreshold:
                start2b=i+1
                break
    #Locate ORF2b end
    end2b=findNextStop(rawGenome,start2b-1,len(rawGenome))
    #add ORF2b to dictionary
    openReadingFrames["2b"]=str(start2b)+","+str(end2b)
    
    #Locate ORF3 start
    for i in range(end2a-190,end2a):
        if rawGenome[i:i+3]=='aug':
            seq=rawGenome[i-10:i+15]
            score=scoreSequence(seq,allORFmatrix)
            if score>allORFScoreThreshold:
                start3=i+1
                break
    #Locate ORF3 end
    end3=findNextStop(rawGenome,start3-1,len(rawGenome))
    #Add ORF3 to dictionary
    openReadingFrames["3"]=str(start3)+","+str(end3)
    
    #Locate ORF4 start
    start4=-1
    trace=0
    while start4==-1:
        for i in range(end3-(230+trace),end3-100):
            if rawGenome[i:i+3]=='aug':
                seq=rawGenome[i-10:i+15]
                score=scoreSequence(seq,allORFmatrix)
                if score>allORFScoreThreshold:
                    start4=i+1
                    break
        trace+=10
    #Locate ORF4 end
    end4=findNextStop(rawGenome,start4-1,len(rawGenome))
    #Add ORF4 to dictionary
    openReadingFrames["4"]=str(start4)+","+str(end4)
    
    #Locate ORF5a start
    trace=0
    start5=-1
    while start5==-1:
        for i in range(end4+(3+trace),end4+60):
            if rawGenome[i:i+3]=='aug':
                seq=rawGenome[i-10:i+15]
                score=scoreSequence(seq,allORFmatrix)
                if score>allORFScoreThreshold:
                    start5=i+1
                    break
        if trace==0: trace=-13
        elif trace==-13: trace=50
        else: break
    #Locate ORF5a end
    end5=findNextStop(rawGenome,start5-1,len(rawGenome))
    #Add ORF5a to dictionary
    openReadingFrames["5"]=str(start5)+","+str(end5)
    
    #Locate ORF6 start
    for i in range(end5-30,len(rawGenome)):
        if rawGenome[i:i+3]=='aug':
            seq=rawGenome[i-10:i+15]
            score=scoreSequence(seq,allORFmatrix)
            if score>allORFScoreThreshold:
                start6=i+1
                break
    #Locate ORF6 end
    end6=findNextStop(rawGenome,start6-1,len(rawGenome))
    #Add ORF6 to dictionary
    openReadingFrames["6"]=str(start6)+","+str(end6)
    
    #Locate ORF7 start
    for i in range(end6-30,len(rawGenome)):
        if rawGenome[i:i+3]=='aug':
            seq=rawGenome[i-10:i+15]
            score=scoreSequence(seq,allORFmatrix)
            if score>allORFScoreThreshold:
                start7=i+1
                break
    #Locate ORF7 end
    end7=findNextStop(rawGenome,start7-1,len(rawGenome))
    #Add ORF7 to dictionary
    openReadingFrames["7"]=str(start7)+","+str(end7)
    
    return openReadingFrames

print("Annotating Genomes...")
for i in range(len(fastafiles)):
    print(fastafiles[i])
    annotatedGenome=annotate(isolateGenomeToRNA(path+fastafiles[i]))
    for key in annotatedGenome.keys():
        print(key + ": " + annotatedGenome[key])
    print()

print("Done...")
