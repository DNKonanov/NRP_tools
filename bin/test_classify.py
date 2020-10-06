from sklearn import tree
import csv
import sys
import numpy as np

## Read in features
## Read in features
features = [] 
print('TEST STRING')
with open('etmp.features.tsv','rb') as tsvin:
	tsvin = np.loadtxt(tsvin, delimiter='\t')
	print(tsvin)
	print('TEST2')
	for row in tsvin:
		features.append(row)
## Read in labels
labels = []
with open('etmp.labels.tsv','rb') as tsvin:
	tsvin = np.loadtxt(tsvin, delimiter='\t')
	print(tsvin)
	for row in tsvin:
		labels.append(row)
