from Bio import Entrez
import re
Entrez.email='mgheffel@ksu.edu'
with open('genomeIDs.txt', 'r') as f:
	r=f.readlines()
	print(len(r))
count=0
ids=[]
f=open('genomeIDs.txt','w')
for i in range(len(r)):
	if i<len(r)-500:
		f.write(r[i])
	else:
		ids.append(r[i])
#print(ids)
print('Gathering Genomes')
genomes=[]
genLens=[]
pattern=r'\d+ bp'
index=0
for i in ids:
	print(index)
	handle=Entrez.efetch(db='nucleotide', id=i, rettype='gb')
	record=handle.read()
	genomes.append(record)
	try:
		match=re.search(pattern,record)
		genLens.append(int(match.group(0)[0:-3]))
	except:
		genLens.append(-1)
	index+=1

print('Checking lengths')
index=0
fullGenomes=[]
while index<len(genomes):
	print(index)
	if genLens[index]>14000:
		fullGenomes.append(genomes[index])
	index+=1

print('Creating genome files')
index=0
pattern=r'LOCUS +[A-Z]+\d+'
path='C:\\Users\\Matt\\Desktop\\Marthaler\\genomes\\'
while index<len(fullGenomes):
	try:
		match=re.search(pattern, fullGenomes[index])
		filename=path + match.group(0)[12:] + '.txt'
		with open(filename,'w') as f:
			f.write(fullGenomes[index])
	except:
		pass
	index+=1
print(len(r)-500)