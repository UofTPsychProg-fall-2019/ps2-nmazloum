#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scene-cat problem set for PSY 1210 - Fall 2018

@author: Michael Mack
"""

#%% import block 
import numpy as np
import scipy as sp
import scipy.stats
import os
import shutil

#%%
# copy files from testing room folders to raw data, rename files to include
# testing room letter in the filename
#
testingrooms = ['A','B','C']
for room in testingrooms:
    path="/Users/seyeeh-negarmazloomfarzaghi/Desktop/ps2-nmazloum/"+"testingroom"+room+"/experiment_data.csv"
    newpath="/Users/seyeeh-negarmazloomfarzaghi/Desktop/ps2-nmazloum/"+"rawdata/"+"experiment_data_"+room+".csv"
    shutil.copyfile(path,newpath)


#%%
# read in all the data files in rawdata directory using a for loop
# columns: subject, stimulus, pairing, accuracy, median RT
#
data = np.empty((0,5))
for room in testingrooms:
    #file = open(path,"r")
    path2 = "/Users/seyeeh-negarmazloomfarzaghi/Desktop/ps2-nmazloum/"+"rawdata/"+"experiment_data_"+room+".csv"
    tmp = sp.loadtxt(path2,delimiter=',')
    data = np.vstack([data,tmp])


#%%
# calculate overall average accuracy and average median RT
#
subject = data [:,0]
stimulus = data[:,1]
pairing = data[:,2]
accuracy = data[:,3]
medianrt = data[:,4]

acc_avg = np.mean(accuracy*100)
mrt_avg=np.mean(medianrt)


#%%
# calculate averages (accuracy & RT) split by stimulus using a for loop and an 
# if statement. (i.e., loop through the data to make a sum for each condition, 
# then divide by the number of data points going into the sum)

sum1=0
sum2=0
sum3=0
sum4=0
wordscount=0
facecount=0
for i in range(len(subject)):
    
    if stimulus [i]==1:
        
        wordscount += 1
        #wordscount = wordscount + 1
        
        sum1 = sum1 + accuracy[i]
        sum2 = sum2 + medianrt[i]
    else:
        
        facecount += 1
        sum3 = sum3 + accuracy[i]
        sum4 = sum4 + medianrt[i]
        
wordsACCavg = (sum1 / wordscount) *100
faceACCavg = (sum3 / facecount) * 100
wordsAMRavg = sum2 / wordscount
faceAMRavg = sum4 / facecount
        

# words: 88.6%, 489.4ms   faces: 94.4%, 465.3ms


#%%
# calculate averages (accuracy & RT) split by congruency using indexing, 
# slicing, and numpy's mean function 
# wp - white/pleasant, bp - black/pleasant
# (hint: only one line of code is needed per average)
#
wp1 = list()
wp2 = list()
bp1 = list()
bp2 = list()

for i in range(len(subject)):
    
    if pairing[i] == 1:
        wp1.append(accuracy[i])
        wp2.append(medianrt[i])
    else:
        bp1.append(accuracy[i])
        bp2.append(medianrt[i])
        
acc_wp = np.mean(wp1)*100
acc_bp = np.mean(bp1)*100
mrt_wp= np.mean(wp2)
mrt_bp = np.mean(bp2)

acc_wp = ...  # 94.0%
acc_bp = ...  # 88.9%
mrt_wp = ...  # 469.6ms
mrt_bp = ...  # 485.1ms


#%% 
# calculate average median RT for each of the four conditions
# use for loops, indexing/slicing, or both!
# (hint: might be easier to slice data into separate words and faces datasets)
#
wwp = list()
wbp = list()
fwp = list()
fbp = list()

for i in range(len(subject)):
    
    if pairing [i] == 1 and stimulus[i] ==1:
        wwp.append (medianrt[i])
        
    elif pairing [i] == 2 and stimulus[i] == 1:
        wbp.append (medianrt[i])
        
    elif pairing[i] == 1 and stimulus[i] == 2:
        fwp.append (medianrt[i])
        
    elif pairing [i] == 2 and stimulus [i] == 2:
        fbp.append (medianrt[i])
        
mrt_wwp = np.mean(wwp)
mrt_wbp = np.mean(wbp)

mrt_fwp = np.mean(fwp)
mrt_fbp = np.mean(fbp)


#converting to array for ttest

wwp_array = np.asarray(wwp)
wbp_array = np.asarray(wbp)

fwp_array = np.asarray(fwp)
fbp_array = np.asarray(fbp)



        


# words - white/pleasant: 478.4ms
# words - black/pleasant: 500.3ms
# faces - white/pleasant: 460.8ms
# faces - black/pleasant: 469.9ms


#%%        
# compare pairing conditions' effect on RT within stimulus using scipy's 
# paired-sample t-test: scipy.stats.ttest_rel()
#
import scipy.stats

t1 = scipy.stats.ttest_rel(wwp_array, wbp_array)
t2 = scipy.stats.ttest_rel(fwp_array, fbp_array)



...

# words: t=-5.36, p=2.19e-5
# faces: t=-2.84, p=0.0096


#%%
# print out averages and t-test results
# (hint: use the ''.format() method to create formatted strings)
#

# you need to print this out for all parts of question 3. have it labelled clearly 
print('\nOVERALL: {:.1f}%, {:.1f} ms'.format(acc_avg,mrt_avg))

print('\nAVERAGES_STIMULUS: {:.1f}%, {:.1f}%, {:.1f} ms, {:.1f} ms'.format(wordsACCavg,faceACCavg,wordsAMRavg,faceAMRavg))

print('\nAVERAGES_PAIRING: {:.1f}%, {:.1f}%, {:.1f} ms, {:.1f} ms'.format(acc_wp,acc_bp, mrt_wp,mrt_bp))

print('\nAVERAGES_FOUR_CONDITIONS:  {:.1f} ms, {:.1f} ms, {:.1f} ms, {:.1f} ms'.format(mrt_wwp,mrt_wbp,mrt_fwp,mrt_fbp))


print("t test [wwp_array,wbp_array]", t1)
print("t test [fwp_array, fbp_array]", t2)


#