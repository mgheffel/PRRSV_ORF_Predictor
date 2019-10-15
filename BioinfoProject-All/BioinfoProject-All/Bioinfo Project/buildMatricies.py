#imports and reading in files + path setup
import os
import re
import pandas as pd
flatfiles=[]
for filename in os.listdir("inputGenomes"):
    flatfiles.append(filename)

path=os.getcwd()+'\\inputGenomes\\'

def parseFile(filename, path):
    fileString=open(path+filename).read()
    #AQUIRE START INDICIES
    #CDS tag marks open reading frame annotation
    startIndex=fileString.find("CDS ")
    #arrays for storing annotated genome ORF starts and ends 
    starts=[]
    ends=[]
    #loop to find all annotated ORFs
    while startIndex>0:
        #only search short segement of file after CDS tag
        temp=fileString[startIndex:startIndex+50]
        #regualar expression to find start index
        startEx=r'[ join(]\d+'
        #regualar expression to find end index
        endEx=r'\.\.\d+'
        #print(temp.replace(" ","")[3])
        try:
            starts.append(int(re.search(startEx,temp).group(0)[1:]))
        except:
            if not temp.replace(" ","")[3]=='<':
                print(filename)
        #move to next CDS tag --is set to -1 if no more CDS tags are found, ending the loop
        startIndex=fileString.find("CDS ",startIndex+1)
    
    #ISOLATE GENOME
    genome=fileString[fileString.find("ORIGIN")+6:len(fileString)]
    genome=re.sub(r'[^a-z]+',"",genome)
    
    #CHANGE T BASE TO U TO REPRESENT AS RNA
    genomeRNA=genome.replace("t","u")
    
    #Get all start codons in genome
    startCodonIndicies=[]
    for i in range(len(genomeRNA)):
        if genomeRNA[i:i+3]=='aug':
            startCodonIndicies.append(i+1)
    
    #only add start codon to false starts if it is not in the true starts array
    falsePositiveStarts=[]
    for index in startCodonIndicies:
        if index in starts:
            pass
        else:
            falsePositiveStarts.append(index)
    trueStarts=[]
    for i in range(len(starts)):
        #seperate orf1 start from all other true starts
        if i==0:
            index=starts[i]
            orf1Start=genomeRNA[index-11:index+14]
        else:
            index=starts[i]
            code=genomeRNA[index-11:index+14]+filename
            trueStarts.append(genomeRNA[index-11:index+14]+filename+ " " + str(index))
    falseStarts=[]
    for index in falsePositiveStarts:
        falseStarts.append(genomeRNA[index-11:index+14]+filename+ " " + str(index))
    return {"trueStarts": trueStarts, "falseStarts": falseStarts, "orf1Start": orf1Start}
	
def buildSpaces():
    trueStartSpace=[]
    falseStartSpace=[]
    orf1StartSpace=[]
    #iterate through all flat files and parse them into start codon spaces
    for i in range(len(flatfiles)):
        space=parseFile(flatfiles[i],path)
        trueStartSpace+=space["trueStarts"]
        falseStartSpace+=space["falseStarts"]
        orf1StartSpace.append(space["orf1Start"])
    return {"trueStartSpace": trueStartSpace, "falseStartSpace": falseStartSpace, "orf1StartSpace": orf1StartSpace}

def buildScoringMatrix(startSpace, numNucs):
    #dictionary of positions and their corrosponding nucleotide values
    startIdentities={"A": [0 for i in range(numNucs)], "U": [0 for i in range(numNucs)], "G": [0 for i in range(numNucs)], "C": [0 for i in range(numNucs)]}
    #for each sequence in orf1 start sequences add nucleoides at each position to matrix
    for i in startSpace:
        seq=i.upper()
        for n in range(numNucs):
            try:
                startIdentities[seq[n]][n]+=1
            except:
                #prints exception if a sequence of invalid length is used
                print("except: " + seq)
    #store as dataframe for formatting
    startIDsDF=pd.DataFrame.from_dict(startIdentities)
    #print(startIDsDF.T)
    
    #turn identities matrix into frequency matrix
    for key in startIdentities:
        for i in range(numNucs):
            startIdentities[key][i]=round(startIdentities[key][i]/len(startSpace),2)
    startFreqsDF=pd.DataFrame.from_dict(startIdentities)
    #print(startFreqsDF.T)
    return(startFreqsDF)

def scoreSequence(seq, matrix):
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

#function to help evaluate threshold scores
def getAllScores(trueSpace, falseSpace, matrix):
    trueHigh=0
    trueLow=999
    trueAVG=0
    falseHigh=0
    falseLow=999
    falseAVG=0
    
    #false start scores
    for i in range(len(falseSpace)):
        s=scoreSequence(falseSpace[i],matrix)
        #if s>.39: print(falseStartSpace[i] + str(s))
        if s>falseHigh: falseHigh=s
        if s<falseLow: falseLow=s
        falseAVG+=s
        #print(s)
    falseAVG/=len(falseSpace)
    #true start scores
    for i in range(len(trueSpace)):
        s=scoreSequence(trueSpace[i],matrix)
        #if s<falseHigh: print(trueSpace[i] + str(s))
        if s>trueHigh: trueHigh=s
        if s<trueLow: trueLow=s
        trueAVG+=s
        #print(s)
    trueAVG/=len(trueSpace)
    return {"True High": trueHigh, "True Low": trueLow, "True AVG": trueAVG, "False High": falseHigh, "False Low": falseLow, "False AVG": falseAVG}

#build scoring matricies and save to pickle for fast load & use in annotation script
print("Building Sequence Spaces...")
spaces=buildSpaces()
print("Building Matricies...")
orf1matrix=buildScoringMatrix(spaces["orf1StartSpace"], 25)
orf1matrix.to_pickle("loadData\\orf1matrix.pckl")
print(orf1matrix.T)
allORFmatrix=buildScoringMatrix(spaces["trueStartSpace"], 25)
allORFmatrix.to_pickle("loadData\\allORFmatrix.pckl")
print(allORFmatrix.T)

#determine score threshold values
print("Determining Score Thresholds...")
orf1Scores=getAllScores(spaces["orf1StartSpace"], spaces["falseStartSpace"], orf1matrix)
orf1Threshold=(orf1Scores["True Low"]+orf1Scores["False High"])/2
print(orf1Scores)
allORFScores=getAllScores(spaces["trueStartSpace"], spaces["falseStartSpace"], allORFmatrix)
print(allORFScores)
allORFThreshold=allORFScores["True Low"]
with open("loadData\\thresholdScores.txt",'w') as f:
    f.write(str(orf1Threshold)+"\n")
    f.write(str(allORFThreshold)+"\n")
print("Done...")