import re

def buildScoreMatrix(sequences):
	#build empty matrix
	nucCount=[]
	for i in range(14): nucCount.append([0,0,0,0])
	#count nucleotides in positions from all sequences
	for seq in sequences:
		#take start codon out of sequence for scoring
		seq=seq[0:10]+seq[13:]
		for n in range(len(seq)):
			if seq[n]=='a': nucCount[n][0]+=1
			elif seq[n]=='t': nucCount[n][1]+=1
			elif seq[n]=='g': nucCount[n][2]+=1
			else: nucCount[n][3]+=1
	probabilityMatrix=[]
	for pos in nucCount:
		posProb=[]
		for i in range(len(pos)):
			posProb.append(pos[i]/len(sequences))
		probabilityMatrix.append(posProb)
	return probabilityMatrix
		
def printScoreMatrix(matrix):
	for i in range(len(matrix)):
		if i==10:
			print('--ATG--')
		print('n' + str(i+1) + ': A: ' + str(matrix[i][0])[0:5] + '% T: ' + str(matrix[i][1])[0:5] + '% G: ' + str(matrix[i][2])[0:5] + '% C: ' + str(matrix[i][3])[0:5] + '%')
		
def scoreSequence(seq, scoreMatrix):
	score=0
	#removed atg from sequence
	seq=seq[0:10]+seq[13:]
	for i in range(len(seq)):
		nuc=seq[i]
		if nuc=='a':
			s=scoreMatrix[i][0]
		elif nuc=='t':
			s=scoreMatrix[i][1]
		elif nuc=='g':
			s=scoreMatrix[i][2]
		else:
			s=scoreMatrix[i][3]
		if s==0: s=.0001
		#score*=s
		score+=s
	return score

def getATG(genome):
	positions=[]
	for i in range(len(genome)):
		if genome[i:i+3]=='atg':
			positions.append(i+1)
	return positions


with open('Yeast.txt','r') as f:
	rawData=f.read()

#get all protien loci
pattern=r'CDS\s+\d+\.\.\d+'
match=re.findall(pattern,rawData)
protienPositions=[]
for CDS in match:
	getNums=re.findall(r'[0-9]+', CDS)
	protienPositions.append(getNums[0] + '-' + getNums[1])

#Extract genome from raw data
startIndex=rawData.index('ORIGIN')+6
genome=''
ignore=['0','1','2','3','4','5','6','7','8','9','\t','\n',' ','/']
for c in rawData[startIndex:]:
	if c not in ignore:
		genome+=c
#Get Initiation sequences
initiationSequences=[]
for pos in protienPositions:
	nums=pos.split('-')
	startIndex=int(nums[0])
	initiationSequence=genome[startIndex-11:startIndex+6]
	initiationSequences.append(initiationSequence)

scoreMatrix=buildScoreMatrix(initiationSequences)
printScoreMatrix(scoreMatrix)

avg=0
low=10000
for seq in initiationSequences:
	score=scoreSequence(seq, scoreMatrix)
	avg+=score
	if score<low: low=score
avg/=len(initiationSequences)
print('True start score: ' + str(avg))
print("True low score " + str(low))

i=0
falseATG=[]
protienStarts=[]
for start in initiationSequences:
	protienStarts.append(start.split('-')[0])
for atg in getATG(genome):
	if atg not in protienStarts:
		falseATG.append(atg)

nonInitiationSequences=[]
for atg in falseATG:
	nonInitiationSequences.append(genome[atg-11:atg+6])
avg=0
high=0
for seq in initiationSequences:
	score=scoreSequence(seq, scoreMatrix)
	avg+=score
	if score>high: high=score
avg/=len(nonInitiationSequences)
print('False start avg score: ' + str(avg))
print("False high score " + str(high))

