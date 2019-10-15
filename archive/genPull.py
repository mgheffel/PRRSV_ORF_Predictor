from Bio import Entrez
import re
Entrez.email='mgheffel@ksu.edu'
handle=Entrez.esearch(db='nucleotide', term='prrsv', idtype='acc') #, retmax='27989'
record=Entrez.read(handle)
handle.close()
print(record['IdList'])
ids=record['IdList']

genomes=[]
genLens=[]
pattern=r'\d+ bp'
for i in ids:
	handle=Entrez.efetch(db='nucleotide', id=i, rettype='gb')
	record=handle.read()
	genomes.append(record)
	try:
		match=re.search(pattern,record)
		genLens.append(int(match.group(0)[0:-3]))
	except:
		genLens.append(-1)
	#genomes.append('\n |||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n')


index=0
toPrint=[]
while index<len(genomes):
	print(genLens[index])
	print(genLens[index]>14000)
	if genLens[index]>14000:
		print('in')
		toPrint.append(genomes[index])
		toPrint.append('|||')
		#f.write(genomes[index])
		#f.write('|||')
	index+=1
with open('EntrezRead.txt','w') as f:
	f.writelines(toPrint)
print(toPrint)

#SAVE ALL PRRSV

