import pickle
f=open('prrsvDataSet.pckl','rb')
allGenomeData=pickle.load(f)
f.close

def getSequenceFromIndex(index, genome):
	return genome[index-11:index-1]+genome[index+2:index+12]

def getAllTrueStarts(filter=''):
	startSequences=[]
	for key in allGenomeData.keys():
		if key[0:len(filter)]==filter:
			genome=allGenomeData[key]['genome']
			for index in allGenomeData[key]['true_start_positions']:
				seq=getSequenceFromIndex(index, genome)
				if seq not in startSequences and len(seq)==20:
					startSequences.append(seq)
	return startSequences

def getAllFalseStarts(filter=''):
	falseSequences=[]
	for key in allGenomeData.keys():
		if key[0:len(filter)]==filter:
			genome=allGenomeData[key]['genome']
			for index in allGenomeData[key]['false_start_positions']:
				seq=getSequenceFromIndex(index, genome)
				if seq not in falseSequences and len(seq)==20:
					falseSequences.append(seq)
	return falseSequences

def printScoreMatrix(matrix):
	for i in range(len(matrix)):
		if i==10:
			print('--ATG--')
		print('p' + str(i+1) + ': A: ' + str(matrix[i][0])[0:5] + '% T: ' + str(matrix[i][1])[0:5] + '% G: ' + str(matrix[i][2])[0:5] + '% C: ' + str(matrix[i][3])[0:5] + '%')
		
def getConsensusSequence(scoreMatrix):
	seq=''
	count=0
	for row in scoreMatrix:
		if count==10: seq+='ATG'
		max=0
		for i in range(1,4):
			if row[i]>row[max]:
				max=i
		if max==0: seq+='a'
		elif max==1: seq+='t'
		elif max==2: seq+='g'
		else: seq+='c'
		count+=1
	return seq
		
def buildScoreMatrix(sequences):
	#build empty count matrix
	nucCount=[]
	for i in range(20): nucCount.append([0,0,0,0])
	#count nucleotides in positions from all sequences
	for seq in sequences:
		for n in range(20):
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

def scoreCodon(sequence, scoreMatrix):
	score=0
	for i in range(20):
		nuc=sequence[i]
		if nuc=='a':
			s=scoreMatrix[i][0]
		elif nuc=='t':
			s=scoreMatrix[i][1]
		elif nuc=='g':
			s=scoreMatrix[i][2]
		else:
			s=scoreMatrix[i][3]
		#s-=.25
		score+=s
	return score

def printScoreData(ls):
	scores=[]
	for seq in ls:
		scores.append(scoreCodon(seq, scoreMatrix))
	max=-10
	min=10
	avg=0
	for i in range(len(scores)):
		score=scores[i]
		if score>max: max=score
		if score<min: min=score
		avg+=score
	avg/=len(scores)
	print('Max score: ' + str(max))
	print('Min score: ' + str(min))
	print('Avg score: ' + str(avg)[0:5])


trueStartSequences=getAllTrueStarts('G')
falseStartSequences=getAllFalseStarts('G')
scoreMatrix=buildScoreMatrix(trueStartSequences)
printScoreMatrix(scoreMatrix)
print('Consensus Sequence: ' + getConsensusSequence(scoreMatrix))

print('True')
printScoreData(trueStartSequences)
print('False')
printScoreData(falseStartSequences)