# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 11:04:41 2019

@author: jason
"""

from processing_functions import bootstrap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os


af_df=pd.read_csv('D:\\snp_data\\kelly_hughes\\kelly_hughes_2019.csv',sep=',')

#width of p intervals
delt=0.025
pvals=[[p, p+delt] for p in np.arange(0.5,1,delt)]
save_path='D:\\snp_data\\kelly_hughes\\boot\\'


#%%
#bootstrapping

for af_cols in [['pRA0','pRA7'],['pRB0','pRB7'],['pRC0','pRC7']]:
    bootstrap(af_df,pvals,1e6,1000,'chom','snp',af_cols,save_path)#,identifier=af_cols[0][:-1])


#%%
#read bootstrap

boot_var_dict={}
boot_mean_dict={}
for filename in os.listdir(save_path):
    if filename[:3]=='var':
        boot_var_dict[filename[4:-4]]=np.loadtxt(save_path+filename)
    elif filename[:4]=='mean':
        boot_mean_dict[filename[5:-4]]=np.loadtxt(save_path+filename)

#%%

for af_cols in [['pRA0','pRA7'],['pRB0','pRB7'],['pRC0','pRC7']]:
    
    plt.figure()
    ax=plt.gca()
    
    for ind,gen in enumerate(af_cols):
        Cdiff_bootstrap=[]        
        for nextgen in af_cols[ind+1:]:
            #offs=(gen-5000)/10
            offs=0
    
            C_candidates=boot_var_dict[str(gen)+'_'+str(nextgen)]/(np.mean(pvals,1)*(1-np.mean(pvals,1)))
            Cdiff_bootstrap.append(C_candidates[:,0]-C_candidates[:,-5])
            
            #C_candidates=np.transpose(np.transpose(C_candidates)-C_candidates[:,-5])
    
            plt.figure()
            plt.fill_between(np.mean(pvals,1),np.percentile(C_candidates,q=2.5,axis=0),
                            np.percentile(C_candidates,q=97.5,axis=0),color='C2',alpha=0.5)
                        
            plt.plot(np.mean(pvals,1),np.mean(C_candidates,0),label=str(gen)+'->'+str(nextgen))
            
            plt.legend(title="Reference generation",loc="upper left", ncol=2)
            
            #plt.ylim([-0.025,0.025])
    
        
        # color = next(ax._get_lines.prop_cycler)['color']
        # Cdiff_mean=np.mean(Cdiff_bootstrap,1)
        # plt.errorbar(np.arange(ind+1,len(af_cols))+offs,Cdiff_mean,
        #             yerr=[Cdiff_mean-np.percentile(Cdiff_bootstrap,q=2.5,axis=1),np.percentile(Cdiff_bootstrap,q=97.5,axis=1)-Cdiff_mean]
        #             ,color=color)
        # plt.plot(np.arange(ind+1,len(af_cols))+offs,Cdiff_mean,'.',markersize=15,color=color,label=str(gen))
        
        # plt.hlines(0,4995,5010)
        
        # plt.xlabel('Generation',fontsize=14)
        # plt.ylabel(r'Excess $\Delta p$ variance',fontsize=14)
        # plt.ylim([-0.004,0.004])
        #plt.legend(title="Reference generation",loc="upper left", ncol=2)    
       #plt.savefig("Cdiff_"+str(rep)+".pdf")


#%%
#mean s

for rep in range(3):
    
    af_cols=[['pRA0','pRA7'],['pRB0','pRB7'],['pRC0','pRC7']][rep]
    
    for ind,gen in enumerate(af_cols[:-1]):    
        #for nextgen in af_cols[ind+1:]:
        nextgen=af_cols[ind+1]
    
        means=boot_mean_dict[str(gen)+'_'+str(nextgen)]#/(np.mean(pvals,1)*(1-np.mean(pvals,1)))
        
        
        plt.fill_between(np.mean(pvals,1),np.percentile(means,q=5,axis=0),
                         np.percentile(means,q=95,axis=0),color='C2',alpha=0.5)
                    
        plt.plot(np.mean(pvals,1),np.mean(means,0),label=['A','B','C'][rep])
        
        plt.plot([0.5,1],[0,0],'k')
        
        plt.ylim([-0.02,0.02])
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