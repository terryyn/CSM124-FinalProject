#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#################Start with the new one


import sys
import numpy as np
from itertools import product ###For 
import copy
from collections import defaultdict


def add_phase(geno, dic, mcnt, homo):  ###Add phase to the dictionary for counting haplotype freq
    hap1 = []
    hap2 = []
    for jj in range(len(geno)):  ###let's make haplotype counts for each of these
           # print(sbarr[jj,ii])
            if geno[jj] == 0:  ###If it's 0, add 0 to haplotype
                hap1.append(0)
                hap2.append(0)
            if geno[jj] == 1: ###If it's 1, add 1 and 0 to haplotype
                hap1.append(1)
                hap2.append(0)
            if geno[jj] == 2: ###If it's 2, add 2 to haplotype
                hap1.append(1)
                hap2.append(1)
       
    str1 = np.array2string(np.array(hap1))
    str2 = np.array2string(np.array(hap2))
    hap_return = [np.array(hap1), np.array(hap2)]
    hap_all = [str1, str2]
       # print(hap_all)
       # print(sbarr_hap[0:SNP_2look, ii*2:(ii*2 +2)])  ###Need to add 2 to get two columns
    for i in range(len(hap_all)): ###For each element of hap_all
           if hap_all[i] in dic:
               if homo == 0:
                   dic[hap_all[i]] = dic[hap_all[i]] + 1.0 ###This part can be defintely be better. To add 1 for complete phased and 
               else:
                   dic[hap_all[i]] =dic[hap_all[i]] + 1.0/(2.0**mcnt)
           else:                                      ###For incomplete haplotypes, I should add partial probabilities
                dic[hap_all[i]] = 1.0/(2.0**mcnt)  ###
    return hap_return

def get_phase(geno, dic): ###Phase genotypes and return haplotypes
    hap1 = []
    hap2 = []
    for jj in range(len(geno)):  ###let's make haplotype counts for each of these
           # print(sbarr[jj,ii])
            if geno[jj] == 0:  ###If it's 0, add 0 to haplotype
                hap1.append(0)
                hap2.append(0)
            if geno[jj] == 1: ###If it's 1, add 1 and 0 to haplotype
                hap1.append(1)
                hap2.append(0)
            if geno[jj] == 2: ###If it's 2, add 2 to haplotype
                hap1.append(1)
                hap2.append(1)
    hap_return = [np.array(hap1), np.array(hap2)]

    return hap_return

##Get imputated genotypes
unmasked_f = 'exmpl_gpred.txt'

for idx, l in enumerate(open(unmasked_f)):
	pass
num_snps = idx + 1
num_indivs = len(l.split())

print('SNPs: {}'.format(num_snps))
print('Individuals: {}'.format(num_indivs))


###Allocate X
X = np.zeros((num_snps, num_indivs), np.int32)
# Fill X array
conversions = {'0': 0, '1': 1, '2': 2}
for idx, l in enumerate(open(unmasked_f)):
	X[idx] = [np.int32(conversions[ch]) for ch in l.split()]

###Let's make an array or haplotypes, simply multiply by 2
Hap = np.zeros((num_snps, num_indivs *2), np.int32)

SNP_2look = 8

num_iterate = int(num_snps/SNP_2look)
last_SNPnum = num_snps%SNP_2look


for i in range(num_iterate+1):
    if i%50 == 0:
        print("I am at ")
        print(i)
    if i == num_iterate:
        if last_SNPnum == 0: ###For a divisble case
            break
        sbarr = X[SNP_2look * i: SNP_2look *i + last_SNPnum]  ##Last case, fewer than 5 SNPs
        sbarr_hap = Hap[SNP_2look * i: SNP_2look *i + last_SNPnum]
        
        
    else:
        sbarr = X[SNP_2look * i: SNP_2look * (i+1)]  ##Look at SNP blocks
        sbarr_hap = Hap[SNP_2look * i: SNP_2look * (i+1)]
       # print("subarray")
       # print(sbarr)
       # print("subarray for haplotype")
       # print(sbarr_hap)
        if i % 1000 == 0:
            print(i)
    d = defaultdict()
    length = sbarr.shape[0]

    for ii in range(num_indivs):
       # print("genotype") 
       # print(sbarr[:, ii])
        if np.count_nonzero(sbarr[:, ii] == 1) > 1: ##For genotypes with at least one heterozygous locus
#            print("This is the genotype")
 #           print(sbarr[:, ii])
            mcnt = np.count_nonzero(sbarr[:, ii] == 1)  ###Count the occurence of 1
            for p in product([0, 1], repeat=mcnt):  ###every p, we should replace, we umasked value with the masked value 
                    count = 0
                    str_interest = copy.deepcopy(sbarr[:, ii])
                    for jj in range(length):
                        if str_interest[jj] == 1:
                            str_interest[jj] = p[count]
                            count += 1
  #              print("input to add_phase")
   #             print(str_interest)
                    add_phase(str_interest, d, mcnt, 1) ###Add phase for each possible ones
                        
                        ####For each permutated possibility, phase and add it to the table
                        
                        
       # print("Done with this genotype")                
       
    
        else:  ##For genotype blocks with homozygous loci of only one heterozygous locus
        #print("array to string")
        #print(hap1)
        #print(hap2)
        
        ##Add these two haplotypes to dictionary and count them
            temp = add_phase(sbarr[:,ii], d, 0, 0)
       # print(temp[1])
        #print(sbarr_hap[0:SNP_2look, ii *2])
            sbarr_hap[:, ii *2] = temp[0]
            sbarr_hap[:, ii *2 + 1] = temp[1]
            
    
   # print("After adding stuff for homo")
   # print(sbarr_hap)
    

###Now do the same thing for heterozygous-containing genotypes
    for ii in range(num_indivs):
            if np.count_nonzero(sbarr[:, ii] == 1) > 1: ##For masked
                mcnt = np.count_nonzero(sbarr[:, ii] == 1)        
                all_comb = []
                for p in product([0, 1], repeat=mcnt):  ###every p, we should replace, we umasked value with the masked value 
                    count = 0
                    str_interest = copy.deepcopy(sbarr[:, ii])
                    for jj in range(length):
                        if str_interest[jj] == 1:
                            str_interest[jj] = p[count]
                            count += 1
                    temp = get_phase(str_interest, d)
               # print(temp[0])
                  #  all_comb.append(temp[0])
                   # all_comb.append(temp[1])
                    all_comb.append(temp)
           # print(all_comb)
        ###Let's look at each element and get values 
                getmax = []
                for gg in range(len(all_comb)):
           #     print(all_comb[gg])
                    getmax.append(d[np.array2string(all_comb[gg][0])] * d[np.array2string(all_comb[gg][1])])
                    max_index = np.argmax(getmax)
            #print("for genotype ")
           # print(sbarr[:, ii])
           # print(getmax)
           # print(max_index)
            #print(sbarr[:, ii] - all_comb[max_index])

                sbarr_hap[:, ii *2] = all_comb[max_index][0]   ###Find the most common haplotype
                sbarr_hap[:, ii *2 + 1] = all_comb[max_index][1] ###Find its haplotype pair
                
np.savetxt('exmpl_phased.txt', Hap, fmt='%d', delimiter=" ")
