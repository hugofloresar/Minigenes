#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 15:38:21 2023

@author: hugof
"""


from Bio.KEGG import REST
import pandas as pd

#### Test aa sequences and eliminate duplicates

code = '32213'

with open('Results_core/Combination/'+str(code)+'.fn') as f:
    lines = f.readlines()


newlist = [] 
freqlist = [] 
newlist.append(lines[0])
newlist.append(lines[1])
freqlist.append(1)
nn = len(lines)
### check if the aa sequences are different
for i in range(int(nn/2)-1):
    ele = lines[2+2*i+1]
    pp = len(newlist)
    for k in range(int(pp/2)):
        seq = newlist[2*k+1]
        crit = False
        if seq == ele:
            crit = True
            freqlist[k]+=1 
            break
    if not(crit):
        newlist.append(lines[2+2*i])
        newlist.append(lines[2+2*i+1])
        freqlist.append(1)    

# print(sum(freqlist))

# ccode = open('Results_core/Combination_filtered/'+str(code)+'.fn', "a")
# for i in range(len(newlist)):
#     ccode.write(newlist[i])
# ccode.close()

# cfreq = open('Results_core/Combination_filtered/'+str(code)+'_freq.fn', "a")
# for i in range(len(freqlist)):
#     cfreq.write(str(freqlist[i]))
#     cfreq.write("\n")
# cfreq.close()


'''From the KEGG Ortholog, obtain the metabolic pathways '''


with open('Results_core/Combination_filtered/'+str(code)+'_result.txt') as f:
    lines = f.readlines()


KO_list = []
for line in lines:
    information = line.split()
    if len(information)>1:
       KO_list.append(information[1]) 
    else:
       KO_list.append("No")         

import numpy as np

def count_freq(lista):
    freq0 = {}
    
    for item in lista:
        if item in freq0:
          freq0[item] += 1
        else:
          freq0[item] = 1

    sorted_list = sorted(freq0.items(), key=lambda item:item[0])
    sorted_list = np.array(sorted_list)

    return sorted_list

KO_list_freq = count_freq(KO_list)

#######################################################################################


lista_pathways_freq = []

for item in KO_list_freq:
    KO = item[0]
    if KO == 'No':
        # print('For KO= ',KO)
        # print('unclassified')
        freq = np.sum(np.array(freqlist)[np.where(np.array(KO_list)==KO)[0]])
        lista_pathways_freq.append(['Nomap','unclassified',freq])
    else:
        # print('For KO= ',KO)
        freq = np.sum(np.array(freqlist)[np.where(np.array(KO_list)==KO)[0]])
        # print('Frequency = ', freq)
        result = REST.kegg_get(KO).read()
        try:
            path = result.split('PATHWAY')[1]
            pathw = path.split('BRITE')[0]
            try:
                pathw = pathw.split('MODULE')[0]
            except:
                print(pathw)
            lista = pathw.split('\n')
            nn = len(lista)
            for i in range(nn-1):
                pathway = lista[i]
                pathway = pathway.split('  ')
                lista_pathways_freq.append([pathway[-2],pathway[-1],freq])
        except:
            path = result.split('PATHWAY')[0]
            # print('case2')
            # print('uncharacterized')
            lista_pathways_freq.append(['map00000','uncharacterized',freq])


nn = len(lista_pathways_freq)
for i in range(nn):
    item  = lista_pathways_freq[i]
    parts = item[0].split(' ')
    if len(parts)>1:
        pathw = parts[1]
        lista_pathways_freq[i][0] = pathw

# for item in lista_pathways_freq:    
#     print(item)


def summary_list(lista):
    freq0 = {}
    
    for item in lista:
        if item[0] in freq0:
          freq0[item[0]] += item[2]
        else:
          freq0[item[0]] = item[2]

    sorted_list = sorted(freq0.items(), key=lambda item:item[0])
    sorted_list = np.array(sorted_list)

    return sorted_list

sum_path = summary_list(lista_pathways_freq)


def summary_labels(lista):
    label0 = {}
    
    for item in lista:
        if item[0] in label0:
          label0[item[0]] = item[1]
        else:
          label0[item[0]] = item[1]

    sorted_list = sorted(label0.items(), key=lambda item:item[0])
    sorted_list = np.array(sorted_list)

    return sorted_list

sum_labels = summary_labels(lista_pathways_freq)


df1 = pd.DataFrame({
    "Namepath": summary_labels(lista_pathways_freq)[:,1],
    code: np.array(sum_path[:,1],dtype = int)
    }, index=sum_path[:,0])

print(df1)

