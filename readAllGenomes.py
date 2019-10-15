import os
flatfiles=[]
for filename in os.listdir("genomes"):
	flatfiles.append(filename)

	
#-------------------------------------Begin annotation region

path='C:\\Users\\Matt\\Desktop\\Bioinfo Project\\genomes\\'
filename=flatfiles[2]
file=open(path+filename)
fileString=file.read()
#saveFile=open(filename[0:filename.index(".")] + "Start.txt" , 'w')

#ISOLATE AUG SURROUNDING NUCLEOTIDES
def space(a):
	pre=9
	reconstruct=""
	i=0
	while i<len(a):
		if i==pre or i==pre+3:
			reconstruct=reconstruct+" "
		reconstruct=reconstruct+a[i]
		i+=1
	a1=""
	for char in reconstruct:
		a1=a1+char+" "
	return a1

#AQUIRE START INDICIES
startIndex=fileString.find("CDS")
starts=[]
while startIndex>0:
	temp=fileString[startIndex+16:fileString.find(".",startIndex)]
	if temp[0]!='>' and temp[0]!='<':
		starts.append(int(temp))
	startIndex=fileString.find("CDS",startIndex+1)
print(starts)

#ISOLATE GENOME
genome=""
for char in fileString[fileString.find("ORIGIN"):len(fileString)]:
    if char=='t'or char=='g' or char=='a' or char=='c':
        genome+=char

#CHANGE T BASE TO U TO REPRESENT AS RNA
genomeRNA=""
for char in genome:
	if char=='t':
		genomeRNA+='u'
	else:
		genomeRNA+=char

#Get all start codons in genome
startCodonIndicies=[]
for i in range(len(genomeRNA)):
	#print(genome[i:i+3])
	if genomeRNA[i:i+3]=='aug':
		startCodonIndicies.append(i+1)

falsePositiveStarts=[]
for index in startCodonIndicies:
	if index in starts:
		pass
	else:
		falsePositiveStarts.append(index)

#print True start codons
print('True start codons: ')
for index in starts:
	print(space(genomeRNA[index-10:index+5]))

#print false positive start codons
print('False positives: ')
for index in falsePositiveStarts:
	print(space(genomeRNA[index-10:index+5]))
	
#test environment
print('In both true and false: ')
trueStartSpace=[]
for index in starts:
	trueStartSpace.append(space(genomeRNA[index-10:index+5]))
falseStartSpace=[]
for index in falsePositiveStarts:
	falseStartSpace.append(space(genomeRNA[index-10:index+5]))
for i in trueStartSpace:
	if i in falseStartSpace:
		print(i)

f=open('genDataSet.txt','w')
for i in trueStartSpace:
	f.write(i+',Y\n')
for i in falseStartSpace:
	f.write(i+',N\n')

#------------------------------------------------END annotation region

numNucs=15
#build consensus matrices
print("True start codon identities by position")
startIdentities={"A": [0 for i in range(numNucs)], "U": [0 for i in range(numNucs)], "G": [0 for i in range(numNucs)], "C": [0 for i in range(numNucs)]}
for i in trueStartSpace:
	temp=i.replace(" ","").upper()
	for n in range(numNucs):
		startIdentities[temp[n]][n]+=1
for key in startIdentities:
	print(key + ": " + str(startIdentities[key]))
print("\nTrue start codon frequency by position")
for key in startIdentities:
	for i in range(numNucs):
		startIdentities[key][i]=round(startIdentities[key][i]/len(trueStartSpace),2)
for key in startIdentities:
	print(key + ": " + str(startIdentities[key]))

print("\nFalse start codon identities by position")
falseStartIdentities={"A": [0 for i in range(numNucs)], "U": [0 for i in range(numNucs)], "G": [0 for i in range(numNucs)], "C": [0 for i in range(numNucs)]}
for i in falseStartSpace:
	temp=i.replace(" ","").upper()
	for n in range(numNucs):
		falseStartIdentities[temp[n]][n]+=1
for key in falseStartIdentities:
	print(key + ": " + str(falseStartIdentities[key]))
print("\nFalse start codon frequency by position")
for key in falseStartIdentities:
	for i in range(numNucs):
		falseStartIdentities[key][i]=round(falseStartIdentities[key][i]/len(falseStartSpace),2)
for key in falseStartIdentities:
	print(key + ": " + str(falseStartIdentities[key]))







	