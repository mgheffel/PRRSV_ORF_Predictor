from Bio import Entrez
import re
Entrez.email='mgheffel@ksu.edu'
handle=Entrez.esearch(db='nucleotide', term='prrsv', idtype='acc', retmax='27989')
record=Entrez.read(handle)
handle.close()
print(record['IdList'])
ids=record['IdList']
with open('genomeIDs.txt', 'w') as f:
	for id in ids:
		f.write(id + '\n')
exit() #GOT ALL THE IDS NOW READ THROUGH THEM ELSEWHERE

print('Gathering Genomes')
genomes=[]
genLens=[]
pattern=r'\d+ bp'
index=0
ids=['CP006415']
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
	fullGenomes.append(genomes[index])
	index+=1

print('Creating genome files')
index=0
pattern=r'LOCUS +[A-Z]+\d+'
path='C:\\Users\\Matt\\Desktop\\Marthaler\\genomes\\'
while index<len(fullGenomes):
	try:
		match=re.search(pattern, fullGenomes[index])
		filename=path + match.group(0)[11:] + '.txt'
		with open('Yeast.txt','w') as f:
			f.write(fullGenomes[index])
	except:
		pass
	index+=1