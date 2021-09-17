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

def convert_frac(string):
    temp=string.split(':')
    if int(temp[1])==0:
        return(np.nan)
    else:
        return int(temp[0])/int(temp[1])


af_df=pd.read_csv('D:\\snp_data\\bergland\\6d_v7.3_output.vcf',sep='\t',skiprows=24)
#compute allele frequencies for PA
af_df[af_df.columns[16:]]=af_df[af_df.columns[16:]].applymap(convert_frac)
#Filter for SNPs tagged as used by Bergland et al.
af_df=af_df[list(map(lambda x: x.split(';')[2][-4:]=='TRUE', af_df['INFO']))]

#%%
#bootstrapping

#width of p intervals
delt=0.05
pvals=[[p, p+delt] for p in np.arange(0.5,1,delt)]
af_cols=af_df.columns[16:]
save_path='D:\\snp_data\\bergland\\boot\\'
bootstrap(af_df,pvals,1e6,1000,'#CHR','POS',af_df.columns[16:],save_path)

#%%
#read bootstrap

# save_path='D:\\snp_data\\bergland\\boot\\'
# af_cols=af_df.columns[16:]
# boot_var_dict={}
# boot_mean_dict={}
# for filename in os.listdir(save_path):
#     if filename[:3]=='var':
#         boot_var_dict[filename[4:-4]]=np.loadtxt(save_path+filename)
#     elif filename[:4]=='mean':
#         boot_mean_dict[filename[5:-4]]=np.loadtxt(save_path+filename)

#%%

# plt.figure()
# ax=plt.gca()

# for ind,gen in enumerate(af_cols):
#     Cdiff_bootstrap=[]        
#     for nextgen in af_cols[ind+1:]:
#         #offs=(gen-5000)/10
#         offs=0

#         C_candidates=boot_var_dict[str(gen)+'_'+str(nextgen)]/(np.mean(pvals,1)*(1-np.mean(pvals,1)))
#         Cdiff_bootstrap.append(C_candidates[:,0]-C_candidates[:,-5])
        
#         #C_candidates=np.transpose(np.transpose(C_candidates)-C_candidates[:,-5])

#         # plt.figure()
#         # plt.fill_between(np.mean(pvals,1),np.percentile(C_candidates,q=2.5,axis=0),
#         #                 np.percentile(C_candidates,q=97.5,axis=0),color='C2',alpha=0.5)
                    
#         # plt.plot(np.mean(pvals,1),np.mean(C_candidates,0),label=str(gen)+'->'+str(nextgen))
        
#         # plt.legend(title="Reference generation",loc="upper left", ncol=2)
        
#         # plt.ylim([0,0.1])

    
#     color = next(ax._get_lines.prop_cycler)['color']
#     Cdiff_mean=np.mean(Cdiff_bootstrap,1)
#     plt.errorbar(np.arange(ind+1,len(af_cols))+offs,Cdiff_mean,
#                 yerr=[Cdiff_mean-np.percentile(Cdiff_bootstrap,q=2.5,axis=1),np.percentile(Cdiff_bootstrap,q=97.5,axis=1)-Cdiff_mean]
#                 ,color=color)
#     plt.plot(np.arange(ind+1,len(af_cols))+offs,Cdiff_mean,'.',markersize=15,color=color,label=str(gen))
    
#     # plt.hlines(0,4995,5010)
    
#     # plt.xlabel('Generation',fontsize=14)
#     # plt.ylabel(r'Excess $\Delta p$ variance',fontsize=14)
#     # plt.ylim([-0.004,0.004])
#     #plt.legend(title="Reference generation",loc="upper left", ncol=2)    
#    #plt.savefig("Cdiff_"+str(rep)+".pdf")

#%%
#mean delta p

# af_df=af_df.sort_values(by='PA_7_2011')

# deltp=af_df['PA_11_2011']-af_df['PA_7_2011']
# plt.scatter(af_df['PA_7_2011'],deltp,s=.1)
# plt.plot([0,1],[0,0],'k')
# plt.plot(af_df['PA_7_2011'],deltp.rolling(2000).mean(),'k')

