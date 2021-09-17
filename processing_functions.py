# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 09:46:05 2020

@author: jason
"""

import numpy as np

def win_contributions(pwin,gen,af_df,deltaf_df):
    query_major=((af_df[gen]<=pwin[1]) & (af_df[gen]>=pwin[0]))
    query_minor=((af_df[gen]>=1-pwin[1]) & (af_df[gen]<=1-pwin[0]))
    query=query_minor | query_major
    return np.array([np.sum(deltaf_df[query]**2),np.sum(deltaf_df*query_major-deltaf_df*query_minor),np.sum(query)])


def mean_deltp(pwin,gen,af_df,deltaf_df):
    query=((af_df[gen]<=pwin[1]) & (af_df[gen]>=pwin[0]))
    return np.mean(deltaf_df[query])

def mean_deltp_folded(pwin,gen,af_df,deltaf_df):
    query_major=((af_df[gen]<=pwin[1]) & (af_df[gen]>=pwin[0]))
    query_minor=((af_df[gen]>=1-pwin[1]) & (af_df[gen]<=1-pwin[0]))
    query=query_minor | query_major
    return np.sum(deltaf_df*query_major-deltaf_df*query_minor)/np.sum(query)



def bootstrap(af_df,pvals,win_size,B,chr_col,pos_col,af_cols,save_path,rep_string):
    
    #queries to filter dataframe by chromosome
    if chr_col=='':
        chrom_queries=[(af_df[pos_col]>-1)]
    else:
        chrom_queries=[(af_df[chr_col]==_) for _ in af_df[chr_col].unique()]
    
    #number of windows in each chromosome
    num_wins=[int(np.max(af_df[pos_col][_])/win_size) for _ in chrom_queries]

    for ind,gen in enumerate(af_cols):
        for nextgen in af_cols[ind+1:]:
            deltaf=af_df[nextgen]-af_df[gen]
            
            win_vals=[]
            for i in range(len(num_wins)):
                for win in range(num_wins[i]):
                    win_query=(af_df[pos_col][chrom_queries[i]]>win*win_size) & (af_df[pos_col][chrom_queries[i]]<(win+1)*win_size)
            
                    win_vals.append([win_contributions(p,gen,af_df[chrom_queries[i]][win_query],deltaf[chrom_queries[i]][win_query]) for p in pvals])
            
            win_vals=np.array(win_vals)
            
            bootstrap_sample_ind=np.random.randint(len(win_vals),size=[B,len(win_vals)])
            
            bootstrap_sample=np.array([np.sum(win_vals[_,:,0],0)/np.sum(win_vals[_,:,2],0)
                                      -(np.sum(win_vals[_,:,1])/np.sum(win_vals[_,:,2]))**2 
                                      for _ in bootstrap_sample_ind])
            
            #variances
            np.savetxt(save_path+'var'+rep_string+'_'+str(gen)+'_'+str(nextgen)+'.tab',bootstrap_sample)
            
            #means
            bootstrap_sample=np.array([np.sum(win_vals[_,:,1],0)/np.sum(win_vals[_,:,2],0) for _ in bootstrap_sample_ind])
            np.savetxt(save_path+'mean'+rep_string+'_'+str(gen)+'_'+str(nextgen)+'.tab',bootstrap_sample)