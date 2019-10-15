with open('falseStartDataSet.txt','r')as f:
	dataset=f.readlines()
temp=[]
for line in dataset:
	if len(line)<10:
		print('!!!!!!')
	else:
		temp.append(line)
dataset=temp
print(len(dataset))
l=len(dataset)
#arrays nums of [a,u,g,c]
u0=[0,0,0,0]
u1=[0,0,0,0]
u2=[0,0,0,0]
u3=[0,0,0,0]
u4=[0,0,0,0]
u5=[0,0,0,0]
u6=[0,0,0,0]
u7=[0,0,0,0]
u8=[0,0,0,0]
d0=[0,0,0,0]
d1=[0,0,0,0]
d2=[0,0,0,0]
loci=[u0,u1,u2,u3,u4,u5,u6,u7,u8,d0,d1,d2]
#indexes to iterate through on the lines in dataset
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

def printData():
	count=0
	for locus in loci:
		print('n' + str(count) + ': A: ' + str(locus[0]/l)[0:4] + ' --U: ' + str(locus[1]/l)[0:4] + ' --G: ' + str(locus[2]/l)[0:4] + ' --C: ' + str(locus[3]/l)[0:4])
		count+=1
#returns matrix of rows of percents [a,u,g,c]
def buildPercentMatrix():
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
printData()
print('|||||||||||||||||||||||||||||||||||||||||||||||||||')
matrix=buildPercentMatrix()
def scoreCodon(sequence):
	score=0
	count=0
	for nuc in sequence:
		if nuc=='a':
			s=float(matrix[count][0])-.25
		elif nuc=='u':
			s=float(matrix[count][1])-.25
		elif nuc=='g':
			s=float(matrix[count][2])-.25
		else:
			s=float(matrix[count][3])-.25
		score+=s*s*s*s
		count+=1
	return score
def cleanLine(line):
	nucs=[]
	for i in index:
		nucs.append(line[i])
	return nucs
#scoreCodon(cleanLine(dataset[0]))
min=100
#for line in dataset:
#	temp=scoreCodon(cleanLine(line))
#	if temp>.004:
#		print(temp)
