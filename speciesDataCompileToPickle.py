import genomeDataExtract as gde
import os
import pickle

path=path='C:\\Users\\Matt\\Desktop\\Marthaler\\genomes\\'

allGenomes={}
for filename in os.listdir('C:\\Users\\Matt\\Desktop\\Marthaler\\genomes'):
	if filename!='clean':
		try:
			with open(path+filename,'r') as f:
				allGenomes[filename]=gde.buildDataDictionary(f.read())
		except:
			with open(path+filename,'r') as f:
				allGenomes[filename]=gde.buildDataDictionary(f.read())
			print('Error on file: ' + filename)


f=open('prrsvDataSet.pckl', 'wb')
pickle.dump(allGenomes,f)
f.close()
