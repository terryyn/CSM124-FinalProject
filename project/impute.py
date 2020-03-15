# -*- coding: utf-8 -*-

import sys
import numpy as np
from itertools import product ###For 
import copy
from collections import defaultdict

masked_f = 'example_data_1_masked.txt'
unmasked_f = 'example_data_1.txt'

for idx, l in enumerate(open(masked_f)):
	pass
num_snps = idx + 1
num_indivs = len(l.split())

print('SNPs: {}'.format(num_snps))
print('Individuals: {}'.format(num_indivs))


X = np.zeros((num_snps, num_indivs), np.int32)
conversions = {'0': 0, '1': 1, '2': 2, '*': 3}


# Fill X array
for idx, l in enumerate(open(masked_f)):
	X[idx] = [np.int32(conversions[ch]) for ch in l.split()]




# Calculate the allele frequency of 2 for all
#afrqs = np.zeros(num_snps)
#for ii in range(num_snps):
#	unmasked = [X[ii][jj] for jj in range(num_indivs) if X[ii][jj] != 3]
#	afrqs[ii] = np.mean(unmasked) / 2.0 if len(unmasked) else 0.5


####Let's look at the first five SNPs
SNP_2look = 5

num_iterate = int(num_snps/SNP_2look)
last_SNPnum = num_snps%SNP_2look


for i in range(num_iterate+1):
    if i == num_iterate:
        sbarr = X[SNP_2look * i: SNP_2look *i + last_SNPnum]  ##Last case, fewer than 5 SNPs
    else:
        sbarr = X[SNP_2look * i: SNP_2look * (i+1)]
        if i % 1000 == 0:
            print(i)

#########This is for each SNP block
            ###We need to make dictionary  count table for each SNP block
    length = sbarr.shape[0]
    d = defaultdict()

    for ii in range(num_indivs):
        if 3 in sbarr[:, ii]: ##For masked
            mcnt = np.count_nonzero(sbarr[:, ii] == 3)
            for p in product([0, 2], repeat=mcnt):  ###every p, we should replace, we umasked value with the masked value 
                count = 0
                str_interest = copy.deepcopy(sbarr[:, ii])
                for jj in range(length):
                    if str_interest[jj] == 3:
                        str_interest[jj] = p[count]
                        count += 1
                        perm = np.array2string(str_interest)
                        if perm in d:
                            d[perm] = d[perm]+1
        
                        else:
                            d[perm] = 1.0/(2**mcnt)
        
        else: ###Count the number of unmasked SNPs ##For unmasked
            if np.array2string(sbarr[:, ii]) in d:
                d[np.array2string(sbarr[:, ii])] = d[np.array2string(sbarr[:, ii])]+1.0
        
            else:
                d[np.array2string(sbarr[:, ii])] =1.0



####Now we have dictionary of SNP patterns and their likelihoods. 
####Kind of inefficient but let's go back to all masked ones and get the most probabl one

    for ii in range(num_indivs):
        if 3 in sbarr[:, ii]: ##For masked
            mcnt = np.count_nonzero(sbarr[:, ii] == 3)        
            all_comb = []
            for p in product([0, 2], repeat=mcnt):  ###every p, we should replace, we umasked value with the masked value 
                count = 0
                str_interest = copy.deepcopy(sbarr[:, ii])
                for jj in range(length):
                    if str_interest[jj] == 3:
                        str_interest[jj] = p[count]
                        count += 1
            
                all_comb.append(str_interest)
            
        ###Let's look at each element and get values 
            getmax = []
            for gg in range(len(all_comb)):
                getmax.append(d[np.array2string(all_comb[gg])])
                max_index = np.argmax(getmax)
                sbarr[:, ii] = all_comb[max_index]

np.savetxt('exmpl1_gpred.txt', X, fmt='%d', delimiter=" ")



with open(masked_f) as mf, open(unmasked_f) as uf:
	idx = 0
	preds = 0
	correct = 0
	for m, u in zip(mf, uf):
		mvals = m.strip().split()
		uvals = u.strip().split()
		for ii, geno in enumerate(mvals):
			if geno == '*':
				preds += 1
				if str(X[idx][ii]) == uvals[ii]:
					correct += 1
		idx += 1
	print('Accuracy: {}'.format(float(correct) / preds))
