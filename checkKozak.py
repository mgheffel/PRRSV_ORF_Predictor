import math
path='C:\\Users\\Matt\\Desktop\\Marthaler\\genomes\\clean\\AB288356Clean.txt'
with open(path,'r') as f:
	genome=f.read()
genomeRNA=""
for char in genome:
	if char=='g' or char=='a' or char=='c':
		genomeRNA+=char
	elif char=='t':
		genomeRNA+='u'
genome=genomeRNA


index=[0,2,4,6,8,10,12,14,16,28,30,32]
def buildScoreMatrix(filename):
	print(filename)
	with open(filename,'r')as f:
		dataset=f.readlines()
	l=len(dataset)
	loci=[]
	for i in range(12):
		loci.append([0,0,0,0])
	index=[0,2,4,6,8,10,12,14,16,28,30,32]
	
	for line in dataset:
		for i in range(len(index)):
			nuc=line[index[i]]
			if nuc=='a':
				loci[i][0]+=1
			elif nuc=='u':
				loci[i][1]+=1
			elif nuc=='g':
				loci[i][2]+=1
			else:
				loci[i][3]+=1
	
	matrix=[]
	for locus in loci:
		temp=[]
		for item in locus:
			temp.append(str(float(str(item/l)[0:5])))
		matrix.append(temp)
	count=0
	for row in matrix:
		rowString='n'+str(count) + ': '
		print('n'+str(count) + ': A: ' + str(row[0]) + ' --U: ' + str(row[1]) + ' --G: ' + str(row[2]) + ' --C: ' + str(row[3]))
		count+=1
	return matrix

def scoreCodon(sequence):
	score=1
	count=0
	nucWeight=[1,0,0,2,1,2,3,1,2,2,1,1]
	for nuc in sequence:
		if nuc=='a':
			s1=float(matrix[count][0])
			s2=float(matrix2[count][0])
		elif nuc=='u':
			s1=float(matrix[count][1])
			s2=float(matrix2[count][1])
		elif nuc=='g':
			s1=float(matrix[count][2])
			s2=float(matrix2[count][3])
		else:
			s1=float(matrix[count][3])
			s2=float(matrix2[count][3])
		if s1==0:
			score-=1
		else:
			score+=s1
		#s1=(s1-.25)/.25
		#s2=(s1-.25)/.25
		#score*=s1
		count+=1
	return score

matrix=buildScoreMatrix('trueStartDataSet.txt')
matrix2=buildScoreMatrix('falseStartDataSet.txt')

def cleanLine(line):
	nucs=[]
	for i in index:
		nucs.append(line[i])
	return nucs
dataset=[]
def setInfo(filename):
	print(filename)
	with open(filename,'r') as f:
		dataset=f.readlines()
		max=-100
		min=100000
		avg=0
		for line in dataset:
			lineScore=scoreCodon(cleanLine(line))
			if lineScore>8.392:
				print(line)
			avg+=lineScore
			if lineScore<min:
				min=lineScore
			if lineScore>max:
				max=lineScore
		print('Codon max score: ' + str(max))
		print('Codon min score: ' + str(min))
		print('Codon avg score: ' + str(avg/len(dataset)))
setInfo('trueStartDataSet.txt')
trueStarts=dataset
setInfo('falseStartDataSet.txt')
falseStarts=dataset

