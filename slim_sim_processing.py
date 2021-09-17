# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 11:44:29 2020

@author: jason
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
from tqdm import tqdm
import pandas as pd
from processing_functions import bootstrap
import pickle

#%%
#multiple time points AF only
#dataframe construction

N=1000


for directory in os.listdir('D:\\slim\\'):
    path ='D:\\slim\\'+ directory+'\\'
    
    print(path)

    seeds=set([_[-13:] for _ in os.listdir(path) if os.path.isfile(os.path.join(path,_)) and _[-4:] not in ['slim','.zip']])
    #arg=[ind for ind,seed in enumerate(seeds) if [_ for _ in os.listdir('.') if _[-13:]==seed][0][-3:]=='870'][0]
    for filenames in tqdm([[_ for _ in os.listdir(path) if _[-13:]==seed] for seed in seeds]):
        #print(filenames)
        #filenames=[_ for _ in os.listdir('.') if os.path.isfile(os.path.join('.',_))]
    
        af_df=pd.DataFrame({})
        
        for filename in filenames:
            gen=filename[-19:-14]
            temp_df=pd.read_csv(path+filename, header=None, skiprows=2,
                        sep=' ', engine='python',
                       usecols=[0,1,2,3,4,8], 
                       names=[gen+'_sid','bid','typ','pos',gen+'_a',gen],
                       skipfooter=2*N+1)#include when genotypes output
            
            #convert from allele counts to frequency
            temp_df[gen]=temp_df[gen]/(2*N)
            
            if gen=='10000':
                mut_df_rep=temp_df
            else:
                mut_df_rep=pd.merge(mut_df_rep,temp_df,on=['bid','typ','pos'],how='left')
            
        af_df=pd.concat([af_df,mut_df_rep],ignore_index=True)
        
        
        af_df.to_csv(path+'af_df_'+filename[-13:]+'.zip', index=False)
    
    
    print(directory,len(af_df))
    #drop fixed mutations
    #af_df=af_df.dropna()
    
    
    
    # #bootstrapping
    
    # #delt=0.025
    # #pvals=[[p, p+delt] for p in np.arange(0.5,1,delt)]
    
    # pvals=[[0.5,0.6],[0.875,0.925]]
    
    # win_size=1e6
    # num_wins=int(1e8/win_size)
    # B=1000
    
    
    # for filename in tqdm([_ for _ in os.listdir(path) if _[-3:]=='zip']):
    #     seed=filename[-17:-4]
    #     af_df=pd.read_csv(path+filename)
        
    #     bootstrap(af_df,pvals,win_size,B,'','pos',['10000','10010'],path+'\\boot\\',seed)
        




#%%
#calculate total s in first generation

total_s={}

path='D:\\slim\\background_s5em2_u1em8\\'
#path='D:\\slim\\positive_s2em2_u1em9\\'

seeds=list(set([_[-13:] for _ in os.listdir(path) if os.path.isfile(os.path.join(path,_)) and _[-4:] not in ['slim','.zip','_s.p']]))
#seed=seeds[1]
#print(seed)

for seed in seeds:
    filenames=[_ for _ in os.listdir(path) if _[-13:]==seed]
    
    af_df=pd.read_csv(path+'af_df_'+seed+'.zip')
    
    non_neutral_query=(af_df['typ']=='m2') #& ((af_df[gen] > 0.9) | (af_df[gen] < 0.1))
    #query_50=(af_df[gen] < 0.6) & (af_df[gen] > 0.4) & (af_df['typ']=='m1')
    #query_95=((af_df[gen] < 0.96) & (af_df[gen] > 0.94)) | ((af_df[gen] < 0.06) & (af_df[gen] > 0.04)) & (af_df['typ']=='m1')
    #query_single=(af_df[gen] == 5e-4)
    #query=non_neutral_query & query_single #| query_50 | query_95
    query=(af_df['bid']==af_df['bid']) #Dummy all True query
            
    #Specify Gaussian stabilizing selection
    stabilizing=False
    #not used when not stabilizing
    opt=0.
    sigma2=0.5
    
    # def ss_fitness(breed_val,opt,sigma2):
    #     return np.exp(-(breed_val-opt)**2/(2*sigma2))
    
    
    gen='10000'
    nextgen='10010'
        
    relevant_mutations=set(af_df[query][gen+'_sid'])
    
    filename=[filename for filename in filenames if filename[:5]==gen][0]
    
    #haplotypes
    gen_dict={}
    with open(path+filename) as f:
        for line in f:
            if line=='Genomes:\n':
                break
        
        for i,line in tqdm(enumerate(f)):
            gen_dict[i]=set(map(int,line.split()[2:])) & relevant_mutations
    
    breedval_map=dict(zip(af_df[query][gen+'_sid'],af_df[query][gen+'_a']))
    
    #individual fitnesses
    fit_dict={}
    for key in tqdm(range(0,2*N,2)):       
        if stabilizing:
            #individual fitness from Gaussian stabilizing selection
            breed_val=np.sum([breedval_map[_] for _ in gen_dict[key]]+[breedval_map[_] for _ in gen_dict[key+1]])
            fit_dict[int(key/2)]=ss_fitness(breed_val, opt, sigma2)
        else:
            #individual fitness from locus-level selection coefficients
            all_loci=(gen_dict[key] | gen_dict[key+1])
            homozygous_loci=(gen_dict[key] & gen_dict[key+1])
            fit_dict[int(key/2)]=np.prod([1+0.5*(1+(_ in homozygous_loci))*breedval_map[_] for _ in all_loci])
            
    
    mean_fit=np.mean(list(fit_dict.values()))
    
    total_s[seed]=np.zeros([len(af_df),3])
    
    #total selection coefficients
    for ind in tqdm(af_df[query].index):
        mut=af_df.loc[ind,gen+'_sid']
        mean_fit_derived=np.mean([fit_dict[int(_/2)] for _ in gen_dict.keys() if mut in gen_dict[_]])
        mean_fit_ancestral= (mean_fit - mean_fit_derived*af_df.loc[ind,gen])/(1-af_df.loc[ind,gen])
        
        #ensure that the major allele is focal for sigma^2(s) calculation 
        #to be consistent with empirical sigma^2(s) calculation
        sign=-(-1)**(af_df.loc[ind,gen]>0.5)
        
        total_s[seed][ind,0]=af_df.loc[ind,gen]
        total_s[seed][ind,1]=af_df.loc[ind,nextgen]
        total_s[seed][ind,2]=sign*(mean_fit_derived-mean_fit_ancestral)/mean_fit
        #total_s.loc[ind,gen+'_tots']=sign*(mean_fit_derived-mean_fit_ancestral)/mean_fit


with open(path+'total_s.p', 'wb') as file:
    pickle.dump(total_s, file)

#%% 

def var_s_p(svals,pvals):
    return np.array([np.var([_[2] for _ in svals if ((_[0]<=p[1]) & (_[0]>=p[0])) | ((_[0]>=1-p[1]) & (_[0]<=1-p[0]))]) for p in pvals])


def var_s_p_lb(svals,pvals):
    
    p=pvals[1]
    var_p_star=np.var([_[1] - _[0] for _ in svals if ((_[0]<=p[1]) & (_[0]>=p[0])) | ((_[0]>=1-p[1]) & (_[0]<=1-p[0]))])
    C_p_star=var_p_star/(np.mean(p)*(1-np.mean(p)))
    
    p=pvals[0]
    var_p=np.var([_[1] - _[0] for _ in svals if ((_[0]<=p[1]) & (_[0]>=p[0])) | ((_[0]>=1-p[1]) & (_[0]<=1-p[0]))])
    C_p=var_p/(np.mean(p)*(1-np.mean(p)))
    
    return (C_p-C_p_star)/(np.mean(p)*(1-np.mean(p)))




#path='D:\\slim\\positive_s2em2_u1em9\\'
path='D:\\slim\\background_s5em2_u1em8\\'

with open(path+'total_s.p', 'rb') as file:
    total_s = pickle.load(file)


delt=0.05
pvals=[[p, p+delt] for p in np.arange(0.5,1,delt)]
var_s_p_all=np.array([var_s_p(total_s[seed],pvals) for seed in total_s.keys()])


pvals=[[0.5,0.6],[0.875,0.925]]
var_s_p_lb_all=np.array([var_s_p_lb(total_s[seed],pvals) for seed in total_s.keys()])


ps=np.array([p+delt/2 for p in np.arange(0.5,1,delt)])
plt.boxplot(100*var_s_p_all*ps*(1-ps))
plt.ylim([0.004,0.035])


plt.scatter(0.01*np.array(var_s_p_lb_all),np.array(var_s_p_all[:,0]))
plt.plot([0,0.001],[0,0.001])


#var_s_p_dict={seed:var_s_p(total_s[seed],pvals) for seed in total_s.keys()}

#af_df.to_csv('totals_'+filename[-13:]+'.zip', index=False,compression='infer')

#tots_cum=np.cumsum(af_df[[str(_)+'_tots' for _ in range(4995,5005)]],1)

#plt.figure()
#plt.plot(range(4995,5005),[np.sqrt(np.var(af_df[query_50][str(gen)+'_tots'])) for gen in range(4995,5005)])
#plt.plot(range(4995,5005),[np.sqrt(np.var(af_df[query_95][str(gen)+'_tots'])) for gen in range(4995,5005)])
#plt.plot(range(4995,5005),[0.25*np.var(af_df[query_50][str(gen)+'_tots'])-0.95*0.05*np.var(af_df[query_95][str(gen)+'_tots']) for gen in range(4995,5005)])

# plt.figure()
# plt.plot(range(4995,5005),[np.sqrt(np.var(tots_cum[query_50][str(gen)+'_tots'])) for gen in range(4995,5005)])
# plt.plot(range(4995,5005),[np.sqrt(np.var(tots_cum[query_95][str(gen)+'_tots'])) for gen in range(4995,5005)])

# plt.figure()
# plt.plot(range(4995,5005),[0.25*np.var(tots_cum[query_50][str(gen)+'_tots'])-0.95*0.05*np.var(tots_cum[query_95][str(gen)+'_tots']) for gen in range(4995,5005)])


for _ in hist_s_p(total_s[seed],pvals):
    plt.figure()
    plt.plot(_[1][:-1],_[0])
    plt.xlim([-1.5,1])


#%% LD

non_neutral_present=np.array([[(af_df.loc[ind,gen+'_sid'] in gen_dict[_]) for _ in gen_dict] for ind in af_df[non_neutral_query].index])

ind=af_df[query_50].index[0]
mut_present=np.array([(af_df.loc[ind,gen+'_sid'] in gen_dict[_]) for _ in gen_dict])

jointp=np.array([np.sum(mut_present & non_neutral_present[_])/(2*N) for _ in non_neutral_present])

D=jointp-af_df['4995'].loc[ind]*np.array(af_df[non_neutral_query]['4995'])


#%%
#mean s

for gen in range(4995,5005):
    #offs=(gen-30)/20
    #means_bootstrap=[]        
    for nextgen in range(gen+1,gen+2):
        means=boot_mean_dict[str(gen)+'_'+str(nextgen)]/np.mean(pvals,1)
        
        plt.figure()
        plt.fill_between(np.mean(pvals,1),np.percentile(means,q=5,axis=0),
                         np.percentile(means,q=95,axis=0),color='C2',alpha=0.5)
                    
        plt.plot(np.mean(pvals,1),np.mean(means,0),label=str(rep)+":"+str(gen)+'->'+str(nextgen))
        
        #plt.ylim([-0.1,0.1])
        plt.legend()  
        
#%%
#mean delta p

N=1000
#path='D:\\slim\\neutral\\'
path='D:\\slim\\background_s1em2_u1em8\\'
#path='D:\\slim\\positive_s2em2_u1em9\\'

from processing_functions import mean_deltp,mean_deltp_folded

delt=0.02
pvals=np.array([[p, p+delt] for p in np.arange(0.,1-4*delt/5,delt/5)])

mean_deltp_all=[]

seeds=list(set([_[-13:] for _ in os.listdir(path) if os.path.isfile(os.path.join(path,_)) and _[-4:] not in ['slim','.zip','_s.p','.swp']]))

for filenames in tqdm([[_ for _ in os.listdir(path) if _[-13:]==seed] for seed in seeds]):
    
    af_df=pd.DataFrame({})
    
    for filename in filenames:
        gen=filename[-19:-14]
        temp_df=pd.read_csv(path+filename, header=None, skiprows=2,
                    sep=' ', engine='python',
                   usecols=[0,1,2,3,4,8], 
                   names=[gen+'_sid','bid','typ','pos',gen+'_a',gen],
                   skipfooter=2*N+1)#include when genotypes output
        
        #convert from allele counts to frequency
        temp_df[gen]=temp_df[gen]/(2*N)
        
        if gen=='10000':
            mut_df_rep=temp_df
        else:
            mut_df_rep=pd.merge(mut_df_rep,temp_df,on=['bid','typ','pos'],how='left')
    af_df=pd.concat([af_df,mut_df_rep],ignore_index=True)
    
    af_df.loc[af_df['10010'].isna(),'10010']=np.round(af_df[af_df['10010'].isna()]['10000'])
    #af_df=af_df.sort_values(by='10000')
    
    neutral_query=(af_df['typ']=='m1')
    
    deltp=af_df['10010']-af_df['10000']
    #plt.scatter(af_df['10000'],deltp,s=0.1)
    
    win_mean_deltp=np.array([mean_deltp_folded(pval,'10000',af_df[neutral_query],deltp[neutral_query])/(np.mean(pval)*(1-np.mean(pval))) for pval in pvals])
    #win_mean_deltp=np.array([mean_deltp_folded(pval,'10000',af_df[neutral_query],deltp[neutral_query]) for pval in pvals])
    
    mean_deltp_all.append(win_mean_deltp)
    
    #plt.plot(pvals[:,0],win_mean_deltp,'b',alpha=0.15)
    
    #period=1000
    #plt.plot(af_df['10000'][period:],deltp.rolling(period).mean()[period:],'b',alpha=0.15)

mean_deltp_all=np.array(mean_deltp_all)   
#plt.boxplot(mean_deltp_all,positions=pvals[:,0],widths=0.002) 
plt.plot(pvals[:,0],np.mean(mean_deltp_all,0))
plt.xlim([0.,1])
plt.ylim([-0.1,0.1])
plt.plot([0,1],[0,0],'b')