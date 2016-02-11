import pandas as pd
from apriori import *
import os
import re
import csv
from json import dumps

# Opening files
df = pd.read_csv( os.path.join(os.path.dirname(__file__), 'data', 'data.tsv'), header=0, delimiter="\t", quoting=3 )
cluster1 = csv.reader(open(os.path.join(os.path.dirname(__file__), 'data', 'clu1.csv'),"rb"))
cluster2 = csv.reader(open(os.path.join(os.path.dirname(__file__), 'data', 'clu2.csv'),"rb"))
cluster3 = csv.reader(open(os.path.join(os.path.dirname(__file__), 'data', 'clu3.csv'),"rb"))
cluster4 = csv.reader(open(os.path.join(os.path.dirname(__file__), 'data', 'clu4.csv'),"rb"))

# Define position of each pattern recognized in genes

def pos(l, patrn):
	l_pos = []
	l_gene = []
	p = re.compile("%s{1}" %patrn)
	for gene in l:
		for m in p.finditer(gene):
			l_pos.append(m.start())
			l_gene.append(m.group())
	return l_pos, l_gene


# Running a priori algorithm to detect patterns and store them with their position in Json File
def apriori_algo(d, k):

	st = []
	l_pos = []
	l_gene = []
	data_pattern = []
	data_occurence = []
	data_position = []

	for line in d:
		nbr = line[0]
		stri = df["spacer_seq"][int(nbr)-1]
		st.append(stri)

	patterns = Apriori(st, 50)

	for i in patterns:
		data_pattern.append(i)
		data_occurence.append(patterns[i])
		l_pos, l_gene = pos(st,i)
		compte = {}.fromkeys(set(l_pos),0)
		for valeur in l_pos:
			compte[valeur] += 1
		data_position.append(compte)
		
	with open("data/data.json" %k , "w") as file:
		file.write(dumps({'pattern':data_pattern, 'occur' : data_occurence, 'position':data_position}, file, indent=3))
		
# This function is made to detect fuzzy data, a pattern that is present in every cluster and can not be deduced as a pattern
def clean_data():
	with open('data/data_1.json') as data_file:
		data_1 = json.load(data_file)
	with open('data/data_2.json') as data_file:    
   		data_2 = json.load(data_file)
	with open('data/data_3.json') as data_file:    
   		data_3 = json.load(data_file)
	with open('data/data_3.json') as data_file:    
   		data_4 = json.load(data_file)

	l = []

	for i in data_1["pattern"]:
		l.append(i)
	for i in data_2["pattern"]:
		l.append(i)
	for i in data_3["pattern"]:
		l.append(i)
	for i in data_4["pattern"]:
		l.append(i)

	compte = {}.fromkeys(set(l),0)

	for valeur in l:
		compte[valeur] += 1
	l_clean = []
	l_position = []
	for i in compte:
		if compte[i]<4:
			l_clean.append(i)
			try:
				if i in data_1["pattern"]:
					l_position.append(1)
			except:
				pass
			try:
				if i in data_2["pattern"]:
					l_position.append(2)
			except:
				pass
			try:
				if i in data_3["pattern"]:
					l_position.append(3)
			except:
				pass
			try:
				if i in data_4["pattern"]:
					l_position.append(4)
			except:
				pass
			l_clean.append(l_position)
	return l_clean
	
	
#running the algorithm to generate a json file for every cluster
l = []
apriori_algo(cluster1, 1)
apriori_algo(cluster2, 2)
apriori_algo(cluster3, 3)
apriori_algo(cluster4, 4)
l = clean_data()
if len(l) is 0:
	print "no matching pattern"
else:
	print l