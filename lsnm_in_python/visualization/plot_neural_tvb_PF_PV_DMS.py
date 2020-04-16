# ============================================================================
#
#                            PUBLIC DOMAIN NOTICE
#
#       National Institute on Deafness and Other Communication Disorders
#
# This software/database is a "United States Government Work" under the 
# terms of the United States Copyright Act. It was written as part of 
# the author's official duties as a United States Government employee and 
# thus cannot be copyrighted. This software/database is freely available 
# to the public for use. The NIDCD and the U.S. Government have not placed 
# any restriction on its use or reproduction. 
#
# Although all reasonable efforts have been taken to ensure the accuracy 
# and reliability of the software and data, the NIDCD and the U.S. Government 
# do not and cannot warrant the performance or results that may be obtained 
# by using this software or data. The NIDCD and the U.S. Government disclaim 
# all warranties, express or implied, including warranties of performance, 
# merchantability or fitness for any particular purpose.
#
# Please cite the author in any work or product based on this material.
# 
# ==========================================================================



# ***************************************************************************
#
#   Large-Scale Neural Modeling software (LSNM)
#
#   Section on Brain Imaging and Modeling
#   Voice, Speech and Language Branch
#   National Institute on Deafness and Other Communication Disorders
#   National Institutes of Health
#
#   This file (plot_neural_tvb_PF_PV_DMS.py) was created on December 12, 2017.
#
#
#   Author: Antonio Ulloa. Last updated by Antonio Ulloa on December 12, 2017  
# **************************************************************************/

# plot_neural_tvb_PF_PV_DMS.py
#
# Plot neuronal activity data files of TVB nodes during PF, PV, and DMS conditions

# what TVB nodes are the TVB host nodes connected to ?
# the nodes below were taken from an output given by the script:
# "display_hagmanns_brain_connectivity.py"
v1_cxn = [330 , 331 , 334 , 335 , 337 , 339 , 340 , 341 , 342 , 343 , 344 , 346 , 347 , 348 , 349 , 350 , 351 , 352 , 353 ,
          354 , 363 , 373 , 374 , 375 , 376 , 377 , 382 , 838 , 839 , 840 , 842 , 843 , 847 , 848 , 873 , 874 , 877]

v4_cxn = [282 , 297 , 361 , 365 , 368 , 377 , 378 , 379 , 381 , 391 , 392 , 397 , 398 , 401 , 402 , 423 , 424]

it_cxn = [380 , 398 , 401 , 402 , 404 , 412 , 415]

fs_cxn = [0 , 39 , 40 , 42 , 43 , 44 , 45 , 48 , 49 , 50 , 51 , 52 , 53 , 73 , 102 , 103 , 106 , 107 , 114 , 450]

d1_cxn = [41 , 42 , 44 , 48 , 52 , 53 , 54 , 55 , 66 , 67 , 71 , 72 , 73 , 75 , 76 , 77 , 129 , 130 , 445]
d2_cxn = [1 , 10 , 19 , 20 , 21 , 22 , 24 , 39 , 40 , 42 , 43 , 44 , 45 , 46 , 52 , 53 , 54 , 55 , 63 , 64 , 65 , 66 , 67 ,
          70 , 72 , 73 , 74 , 75 , 76 , 77 , 85 , 91 , 92 , 94 , 103 , 143 , 145 , 146 , 150 , 239 , 249 , 395 , 429 , 432 ,
          443 , 444 , 445 , 446 , 447 , 448 , 452 , 464 , 469 , 471 , 472 , 486 , 492 , 493 , 494 , 496]

fr_cxn = [99 , 100 , 108 , 109 , 110 , 112 , 113 , 114 , 115 , 116 , 117 , 118 , 119 , 121 , 122 , 123 , 124 , 167 , 173 ,
          183 , 184 , 192 , 193 , 194 , 195 , 196 , 314 , 315 , 612 , 614 , 615 , 616 , 617 , 618 , 622 , 623 , 624 , 626 ,
          686 , 694 , 695 , 696]


import numpy as np
import matplotlib.pyplot as plt

from scipy import signal

# load file containing TVB nodes electrical activity
pf = np.load('subject_12/output.Fixation_dot_incl_PreSMA_3.0_0.15/tvb_neuronal.npy')
pv = np.load('subject_12/output.PV_incl_PreSMA_w_Fixation_3.0_0.15/tvb_neuronal.npy')
dms = np.load('subject_12/output.DMSTask_incl_PreSMA_w_Fixation_3.0_0.15/tvb_neuronal.npy')
pf_syn   = np.load('subject_12/output.Fixation_dot_incl_PreSMA_3.0_0.15/synaptic_in_998_ROIs.npy')
pv_syn  = np.load('subject_12/output.PV_incl_PreSMA_w_Fixation_3.0_0.15/synaptic_in_998_ROIs.npy')
dms_syn = np.load('subject_12/output.DMSTask_incl_PreSMA_w_Fixation_3.0_0.15/synaptic_in_998_ROIs.npy')
pf_bold   = np.load('subject_12/output.Fixation_dot_incl_PreSMA_3.0_0.15/bold_balloon_998_regions_3T_0.25Hz.npy')
pv_bold  = np.load('subject_12/output.PV_incl_PreSMA_w_Fixation_3.0_0.15/bold_balloon_998_regions_3T_0.25Hz.npy')
dms_bold = np.load('subject_12/output.DMSTask_incl_PreSMA_w_Fixation_3.0_0.15/bold_balloon_998_regions_3T_0.25Hz.npy')

print 'Shape of TVB neuronal activity files: ', pf.shape
print 'Shape of TVB synaptic activity files: ', pf_syn.shape
print 'Shape of TVB BOLD activity files: ', pf_bold.shape


# define neural activity window to look at (typically the duration of one trial)
start = 3750
end = 3850

# the following nodes have connections with host nodes
# ... as determined by script "display_hagmann_brain_connectivity"
pf_v1 = np.mean(pf[start:end, 0, v1_cxn, 0], axis=1)
#pf_v1 = np.mean(pf_bold, axis=0)
#pf_v1 = np.mean(pf_bold[v1_cxn], axis=0)
#pf_v1 = pf_syn[331, start:end]
pf_v4 = np.mean(pf[start:end, 0, v4_cxn, 0], axis=1)
pf_it = np.mean(pf[start:end, 0, it_cxn, 0], axis=1)
pf_fs = np.mean(pf[start:end, 0, fs_cxn, 0], axis=1)
pf_d1 = np.mean(pf[start:end, 0, d1_cxn, 0], axis=1)
pf_d2 = np.mean(pf[start:end, 0, d2_cxn, 0], axis=1)
pf_fr = np.mean(pf[start:end, 0, fr_cxn, 0], axis=1)

pv_v1 = np.mean(pv[start:end, 0, v1_cxn, 0], axis=1)
#pv_v1 = np.mean(pv_bold, axis=0)
#pv_v1 = np.mean(pv_bold[v1_cxn], axis=0)
#pv_v1 = pv_syn[331, start:end]
pv_v4 = np.mean(pv[start:end, 0, v4_cxn, 0], axis=1)
pv_it = np.mean(pv[start:end, 0, it_cxn, 0], axis=1)
pv_fs = np.mean(pv[start:end, 0, fs_cxn, 0], axis=1)
pv_d1 = np.mean(pv[start:end, 0, d1_cxn, 0], axis=1)
pv_d2 = np.mean(pv[start:end, 0, d2_cxn, 0], axis=1)
pv_fr = np.mean(pv[start:end, 0, fr_cxn, 0], axis=1)

dms_v1 = np.mean(dms[start:end, 0, v1_cxn, 0], axis=1)
#dms_v1 = np.mean(dms_bold, axis=0)
#dms_v1 = np.mean(dms_bold[v1_cxn], axis=0)
#dms_v1 = dms_syn[331, start:end]
dms_v4 = np.mean(dms[start:end, 0, v4_cxn, 0], axis=1)
dms_it = np.mean(dms[start:end, 0, it_cxn, 0], axis=1)
dms_fs = np.mean(dms[start:end, 0, fs_cxn, 0], axis=1)
dms_d1 = np.mean(dms[start:end, 0, d1_cxn, 0], axis=1)
dms_d2 = np.mean(dms[start:end, 0, d2_cxn, 0], axis=1)
dms_fr = np.mean(dms[start:end, 0, fr_cxn, 0], axis=1)

# Extract number of timesteps from one of the matrices
timesteps = pf_v1.shape[0]

ts_to_plot = end - start

# what was the duration of simulation in real time (in ms)?
x_lim = (end - start) * 10. * 5. / 1000.

print 'Simulation timesteps to plot = ', timesteps
print 'Time to plot (in seconds) = ', x_lim

print 'Shape of synaptic arrays: ', pf_v1.shape

# Contruct a numpy array of timesteps (data points provided in data file)
real_time = np.linspace(0, (ts_to_plot-1)*50./1000., num=ts_to_plot)

# increase font size
plt.rcParams.update({'font.size': 15})

# Set up plot
plt.figure()

# Plot V1 module
ax = plt.subplot(7,1,7)
ax.plot(real_time, pf_v1, color='blue', linewidth=1)
ax.plot(real_time, pv_v1, color='green', linewidth=1)
ax.plot(real_time, dms_v1, color='red', linewidth=1)
ax.set_yticks([])
#ax.set_xlim([0, x_lim])
#ax.set_ylim([0, 1])
#ax.set_title('SIMULATED ELECTRICAL ACTIVITY, HAGMANNS BRAIN')
plt.ylabel('V1', rotation='horizontal', horizontalalignment='right')

ax = plt.subplot(7,1,6)
ax.plot(real_time, pf_v4, color='blue', linewidth=1)
ax.plot(real_time, pv_v4, color='green', linewidth=1)
ax.plot(real_time, dms_v4, color='red', linewidth=1)
ax.set_yticks([])
ax.set_xticks([])
ax.set_xlim(0, x_lim)
#ax.set_ylim([0, 1])
plt.ylabel('V4', rotation='horizontal', horizontalalignment='right')

ax = plt.subplot(7,1,5)
ax.plot(real_time, pf_it, color='blue', linewidth=1)
ax.plot(real_time, pv_it, color='green', linewidth=1)
ax.plot(real_time, dms_it, color='red', linewidth=1)
ax.set_yticks([])
ax.set_xticks([])
ax.set_xlim(0, x_lim)
#ax.set_ylim([0, 1])
plt.ylabel('IT', rotation='horizontal', horizontalalignment='right')

ax = plt.subplot(7,1,4)
ax.plot(real_time, pf_fs, color='blue', linewidth=1)
ax.plot(real_time, pv_fs, color='green', linewidth=1)
ax.plot(real_time, dms_fs, color='red', linewidth=1)
ax.set_yticks([])
ax.set_xticks([])
ax.set_xlim(0, x_lim)
#ax.set_ylim([0, 1])
plt.ylabel('FS', rotation='horizontal', horizontalalignment='right')

ax = plt.subplot(7,1,3)
ax.plot(real_time, pf_d1, color='blue', linewidth=1)
ax.plot(real_time, pv_d1, color='green', linewidth=1)
ax.plot(real_time, dms_d1, color='red', linewidth=1)
ax.set_yticks([])
ax.set_xticks([])
ax.set_xlim(0, x_lim)
#ax.set_ylim([0, 1])
plt.ylabel('D1', rotation='horizontal', horizontalalignment='right')

ax = plt.subplot(7,1,2)
ax.plot(real_time, pf_d2, color='blue', linewidth=1)
ax.plot(real_time, pv_d2, color='green', linewidth=1)
ax.plot(real_time, dms_d2, color='red', linewidth=1)
ax.set_yticks([])
ax.set_xticks([])
ax.set_xlim(0, x_lim)
#ax.set_ylim([0, 1])
plt.ylabel('D2', rotation='horizontal', horizontalalignment='right')

ax = plt.subplot(7,1,1)
ax.plot(real_time, pf_fr, color='blue', linewidth=1)
ax.plot(real_time, pv_fr, color='green', linewidth=1)
ax.plot(real_time, dms_fr, color='red', linewidth=1)
ax.set_yticks([])
ax.set_xticks([])
ax.set_xlim(0, x_lim)
#ax.set_ylim([0, 1])
plt.ylabel('FR', rotation='horizontal', horizontalalignment='right')

# Show the plot on the screen
plt.show()

