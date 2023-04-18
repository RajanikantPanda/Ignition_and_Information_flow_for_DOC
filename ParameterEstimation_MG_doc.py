#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: matgilson
"""

import os
import numpy as np
import scipy.io as sio
import scipy.signal as spsg
import scipy.stats as stt
import WBLECmodel
import pyMOU.pymou as pymou
import matplotlib.pyplot as plt


#%% create a folder for saving the data and model parameters

res_dir = 'model_param/'
if not os.path.exists(res_dir):
    print('create directory:',res_dir)
    os.makedirs(res_dir)

show_graphs = True

#%% general properties
# load fMRI data
ts_tmp_cnt = sio.loadmat('data_controls_all')['ts_all']
ts_tmp_mcs = sio.loadmat('data_MCS_all_2_fMRI')['ts_all']
ts_tmp_uws = sio.loadmat('data_UWS_all_2_fMRI')['ts_all']
ts_tmp= np.concatenate((ts_tmp_cnt, ts_tmp_mcs, ts_tmp_uws), axis=0)  
print(ts_tmp.shape)

# n_sub= number of subjects (controls + patients)
# N= number of ROIs (brain regions)
# T=number of time points (in TR) of the fMRI/BOLD
n_sub, N,T = np.shape(ts_tmp)

# time shifts for spatiotemporal covariances in our FC: 0 and 1 TR
v_tau = np.arange(2,dtype=float)
n_tau = v_tau.size


#%% functional data
# filtering between 0.01 and 0.1 Hz
n_order = 3
TR = 2.
Nyquist_freq = 0.5 / TR
low_f = 0.01 / Nyquist_freq
high_f = 0.1 / Nyquist_freq
b,a = spsg.iirfilter(n_order,[low_f,high_f],btype='bandpass',ftype='butter')

# calculate FC
ts_emp = np.zeros([n_sub,N,T])
FC_emp = np.zeros([n_sub,n_tau,N,N])
for i_sub in range(n_sub):
    # appply filter
    ts_emp[i_sub,:,:] = spsg.filtfilt(b,a,ts_tmp[i_sub,:,:],axis=1)
    # center the time series
    ts_emp[i_sub,:,:] -= np.outer(ts_emp[i_sub,:,:].mean(1),np.ones(T))

for i_sub in range(n_sub):
    # calculate BOLD covariances for all time shifts
    for i_tau in range(n_tau):
        FC_emp[i_sub,i_tau,:,:] = np.tensordot(ts_emp[i_sub,:,0:T-n_tau+1],ts_emp[i_sub,:,i_tau:T-n_tau+1+i_tau],axes=(1,1)) / float((T-n_tau))

# renormalization of FC (homogeneous for all sessions)
renorm_factor = FC_emp[:,0,:,:].mean()
FC_emp *= 0.5 / renorm_factor
print('mean FC value (most of the distribution should be between 0 and 1):',FC_emp.mean())

# same renormalization for time series
ts_emp *= 0.5 / renorm_factor


#%% structural data

# load DTI data
SC_anat = sio.loadmat('data_controls_all.mat')['sc_all'].mean(0)
print(SC_anat.shape)

# mask for existing connections in EC
mask_EC = np.zeros([N,N],dtype=bool)
# limit DTI value to determine SC (only connections with larger values are tuned)
lim_SC = 0.0067 # SC density value
# threshold DTI to obtain EC mask
mask_EC[SC_anat>lim_SC] = True
for i in range(N):
    # no self connection
    mask_EC[i,i] = False
    # additional interhemispheric connections (known to be missed by DTI)
#    mask_EC[i,i-int(N/2)] = True

print('EC density:',mask_EC.sum()/float(N*(N-1)))

# diagonal mask for input noise matrix (here, no input cross-correlation)
mask_Sigma = np.eye(N,dtype=bool)


#%% test figures for ROI ordering

if show_graphs:
    plt.figure()
    plt.imshow(FC_emp.mean(0)[0,:,:])
    plt.title('FC0')
    
    plt.figure()
    plt.imshow(SC_anat)
    plt.title('SC')
    
    plt.figure()
    plt.imshow(mask_EC)
    plt.title('mask_EC')


#%%  model optimization > EC connectivity

# old version
J_mod = np.zeros([n_sub,N,N]) # Jacobian (off-diagonal elements = EC)
Sigma_mod = np.zeros([n_sub,N,N]) # local variance

# new version
model = pymou.MOU()
J_mod2 = np.zeros([n_sub,N,N]) # Jacobian (off-diagonal elements = EC)
Sigma_mod2 = np.zeros([n_sub,N,N]) # local variance

for i_sub in range(n_sub):
    print('sub',i_sub)

    # OLD VERSION
    print('OLD')

    # objective FC matrices (empirical)
    FC0_obj = FC_emp[i_sub,0,:,:]
    FC1_obj = FC_emp[i_sub,1,:,:]
    # autocovariance for time shifts in v_tau; with lower bound to avoid negative values (cf. log)
    ac_tmp = np.maximum(FC_emp[i_sub,:,:,:].diagonal(axis1=1,axis2=2),1e-6)
    # time constant for BOLD autocovariances (calculated with lag from 0 to 2 TRs)
    # it relates to inverse of negative slope of autocovariance over all ROIs
    tau_x = -1. / np.polyfit(v_tau,np.log(ac_tmp).mean(1),1)[0]
    # this would be th eversion wuth 1 value of tau_x for each ROI
    #tau_x_n = -1. / np.polyfit(v_tau,np.log(ac_tmp),1)[0]
    # perform model inversion
    (J_mod_tmp,Sigma_mod_tmp) = WBLECmodel.optimize(FC0_obj,FC1_obj,tau_x,mask_EC,mask_Sigma)
    tau_save=J_mod_tmp[0][0]
    # store results
    J_mod[i_sub,:,:] = J_mod_tmp  #j_mod effectivie connectivity
    Sigma_mod[i_sub,:,:] = Sigma_mod_tmp
    
    # NEW VERSION
    print('NEW')
    #model.fit(ts_emp,i_tau_opt=1, algo_version='PCB2016') 
    model.fit(ts_emp[i_sub,:,:].T, mask_C=mask_EC, mask_Sigma=mask_Sigma, epsilon_C=0.0005, epsilon_Sigma=0.05)
    J_mod2[i_sub,:,:] = model.J
    Sigma_mod2[i_sub,:,:] = model.Sigma
    print('steps:', model.d_fit['iterations'], 'error (distance):', model.d_fit['distance'],'; Pearson:', model.d_fit['correlation'])


# calculate error and Pearson between estimated EC
err_J = np.zeros([n_sub])
Pearson_J = np.zeros([n_sub])
for i_sub in range(n_sub):
    err_J[i_sub] = np.linalg.norm(J_mod[i_sub,mask_EC]-J_mod2[i_sub,mask_EC]) / np.linalg.norm(J_mod[i_sub,mask_EC])
    Pearson_J[i_sub] = stt.pearsonr(J_mod[i_sub,mask_EC], J_mod2[i_sub,mask_EC])[0]
        
#%% compare results

i_sub = 0
plt.figure()
plt.scatter(J_mod[i_sub,mask_EC], J_mod2[i_sub,mask_EC])
plt.xlabel('old EC')
plt.ylabel('new EC')
plt.title('subject '+str(i_sub))


plt.figure()
plt.violinplot(err_J, positions=[0])
plt.violinplot(Pearson_J, positions=[1])
plt.axis(ymin=0, ymax=1)
plt.xticks([0,1],['error','Pearson'])
plt.show()
    
#%% save results
np.save(res_dir+'FC_emp.npy',FC_emp) # empirical spatiotemporal FC
np.save(res_dir+'mask_EC.npy',mask_EC) # mask of optimized connections
np.save(res_dir+'mask_Sigma.npy',mask_Sigma) # mask of optimized Sigma elements

np.save(res_dir+'J_mod.npy',J_mod2) # estimated Jacobian matrices (EC + inverse time constant on diagonal)
np.save(res_dir+'Sigma_mod.npy',Sigma_mod2) # estimated Sigma matrices

np.save(res_dir+'J_mod_old.npy',J_mod) # estimated Jacobian matrices (EC + inverse time constant on diagonal)
np.save(res_dir+'Sigma_mod_old.npy',Sigma_mod) # estimated Sigma matrices

sio.savemat(res_dir+'J_mod.mat', {'J_mod2':J_mod2})
sio.savemat(res_dir+'Sigma_mod.mat', {'Sigma_mod2':Sigma_mod2})
sio.savemat(res_dir+'FC_emp.mat', {'FC_emp':FC_emp})







