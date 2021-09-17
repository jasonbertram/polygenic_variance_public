# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 11:04:41 2019

@author: jason
"""

from processing_functions import bootstrap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import os

def convert_frac(string):
    temp=string.split('/')
    return int(temp[0])/int(temp[1])

af_df=pd.read_csv('D:\\snp_data\\barghi\\frequencies.tab',sep='\t')
af_df[af_df.columns[3:]]=af_df[af_df.columns[3:]].applymap(convert_frac)
save_path='D:\\snp_data\\barghi\\boot\\'


#%%
#bootstrapping

#width of p intervals
delt=0.025
pvals=[[p, p+delt] for p in np.arange(0.5,1,delt)]

for af_cols in [['R'+str(rep+1)+'_F'+str(gen) for gen in range(0,61,10)] for rep in range(10)]:
    bootstrap(af_df,pvals,1e6,1000,'chr','pos',af_cols,save_path)



#%%

delt=0.025
pvals=[[p, p+delt] for p in np.arange(0.5,1,delt)]

boot_var_dict={}
boot_mean_dict={}
for filename in os.listdir(save_path):
    if filename[:3]=='var':
        boot_var_dict[filename[4:-4]]=np.loadtxt(save_path+filename)
    elif filename[:4]=='mean':
        boot_mean_dict[filename[5:-4]]=np.loadtxt(save_path+filename)



#%%
#C difference between 0.5 and 0.9 AF with 95% block bootstrap ci
#All time differences


for rep in range(1,11):
    plt.figure()
    ax = plt.gca()
    for gen in range(0,51,10):
        offs=(gen-30)/20
        Cdiff_bootstrap=[]        
        for nextgen in range(gen+10,61,10):
            C_candidates=boot_var_dict[str(rep)+'_'+str(gen)+'_'+str(nextgen)]/((pvals+delt/2)*(1-pvals-delt/2))
            Cdiff_bootstrap.append(C_candidates[:,0]-C_candidates[:,-5])
            
        
        color = next(ax._get_lines.prop_cycler)['color']
        Cdiff_mean=np.mean(Cdiff_bootstrap,1)
        plt.errorbar(np.arange(gen+10,61,10)+offs,Cdiff_mean,
                    yerr=[Cdiff_mean-np.percentile(Cdiff_bootstrap,q=2.5,axis=1),np.percentile(Cdiff_bootstrap,q=97.5,axis=1)-Cdiff_mean]
                    ,color=color)
        plt.plot(np.arange(gen+10,61,10)+offs,Cdiff_mean,'.',markersize=15,color=color,label=str(gen))
        
        plt.xlabel('Generation',fontsize=14)
        plt.ylabel(r'$C(0.5)-C(0.9)$',fontsize=14)
        plt.legend(title="Reference generation",loc="upper left", ncol=2)    
        plt.savefig("Cdiff_"+str(rep)+".pdf")

#%%
#mean s

for rep in range(10):
    #plt.figure()
    #ax = plt.gca()
    
    af_cols=['R'+str(rep+1)+'_F'+str(gen) for gen in range(0,61,10)]
    
    for ind,gen in enumerate(af_cols[:-1]):    
        #for nextgen in af_cols[ind+1:]:
        nextgen=af_cols[ind+1]
    
        means=boot_mean_dict[str(gen)+'_'+str(nextgen)]#/(np.mean(pvals,1)*(1-np.mean(pvals,1)))
        
        plt.figure()
        plt.fill_between(np.mean(pvals,1),np.percentile(means,q=5,axis=0),
                         np.percentile(means,q=95,axis=0),color='C2',alpha=0.5)
                    
        plt.plot(np.mean(pvals,1),np.mean(means,0),label=str(rep)+":"+str(gen)+'->'+str(nextgen))
        
        plt.plot([0.5,1],[0,0])
        
        plt.ylim([-0.05,0.05])
        plt.legend()  
            
        
        #color = next(ax._get_lines.prop_cycler)['color']
        #Cdiff_mean=np.mean(Cdiff_bootstrap,1)
        #plt.errorbar(np.arange(gen+10,61,10)+offs,Cdiff_mean,
        #            yerr=[Cdiff_mean-np.percentile(Cdiff_bootstrap,q=2.5,axis=1),np.percentile(Cdiff_bootstrap,q=97.5,axis=1)-Cdiff_mean]
        #             ,color=color)
        # plt.plot(np.arange(gen+10,61,10)+offs,Cdiff_mean,'.',markersize=15,color=color,label=str(gen)+'->'+str(nextgen))
        
        # plt.xlabel('Generation',fontsize=14)
        # plt.ylabel(r'Excess $\Delta p$ variance',fontsize=14)
        # plt.legend(title="Reference generation",loc="upper left", ncol=2)    
        # #plt.savefig("Cdiff_"+str(rep)+".pdf")

#%%
#mean delta p

af_df=af_df.sort_values(by='PA_7_2011')

deltp=pd.DataFrame(af_df['R1_F10']-af_df['R1_F0'],columns=['deltp'])
deltp['R1_F0']=af_df['R1_F0']
deltp=deltp.sort_values(by='R1_F0')
plt.scatter(deltp['R1_F0'],deltp['deltp'],s=.1)
plt.plot([0,1],[0,0],'k')
plt.plot(deltp['R1_F0'],deltp['deltp'].rolling(50000).mean(),'k')
plt.ylim([-0.05,0.05])

