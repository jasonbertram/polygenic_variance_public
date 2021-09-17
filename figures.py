# -*- coding: utf-8 -*-
"""
Created on Sat Dec 19 14:44:32 2020

@author: jason
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import os
import seaborn as sns
import pandas as pd
import pickle
import matplotlib.patches as patches

#save_path=path to bootstrap files

save_path='D:\\snp_data\\barghi\\boot\\'
boot_var_dict_barghi={}
boot_mean_dict_barghi={}
for filename in os.listdir(save_path):
    if filename[:3]=='var':
        boot_var_dict_barghi[filename[4:-4]]=np.loadtxt(save_path+filename)
    elif filename[:4]=='mean':
        boot_mean_dict_barghi[filename[5:-4]]=np.loadtxt(save_path+filename)

save_path='D:\\snp_data\\kelly_hughes\\boot\\'
boot_var_dict_kelly={}
boot_mean_dict_kelly={}
for filename in os.listdir(save_path):
    if filename[:3]=='var':
        boot_var_dict_kelly[filename[4:-4]]=np.loadtxt(save_path+filename)
    elif filename[:4]=='mean':
        boot_mean_dict_kelly[filename[5:-4]]=np.loadtxt(save_path+filename)
        

save_path='D:\\snp_data\\bergland\\boot\\'  
boot_var_dict_bergland={}
boot_mean_dict_bergland={}
for filename in os.listdir(save_path):
    if filename[:3]=='var':
        boot_var_dict_bergland[filename[4:-4]]=np.loadtxt(save_path+filename)
    elif filename[:4]=='mean':
        boot_mean_dict_bergland[filename[5:-4]]=np.loadtxt(save_path+filename)
        

#color-blind friendly schemes adapted from https://personal.sron.nl/~pault/
#mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=['#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE',
#                    '#AA3377', '#BBBBBB']) 
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=['#CC6677', '#332288', '#DDCC77', '#117733', '#88CCEE',
                    '#882255', '#44AA99', '#999933', '#AA4499', '#BBBBBB']) 

        
#%%
#Barghi et al. and Kelly et al. C with 95% block bootstrap ci

#color-blind friendly schemes adapted from https://personal.sron.nl/~pault/
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=['#CC6677', '#332288', '#DDCC77', '#117733', '#88CCEE',
                    '#882255', '#44AA99', '#999933', '#AA4499', '#BBBBBB']) 

delt=0.025
pvals=[[p, p+delt] for p in np.arange(0.5,1,delt)]

fig = plt.figure(figsize=[3.,5],constrained_layout=True,dpi=300)
gs = gridspec.GridSpec(nrows=2, ncols=1, figure=fig)


ax1=fig.add_subplot(gs[0])
ax1.ticklabel_format(axis='y',style='sci',scilimits=(0,0)) 

for rep in range(10):
    af_cols=['R'+str(rep+1)+'_F'+str(gen) for gen in range(0,61,10)]
       
    color = next(ax1._get_lines.prop_cycler)['color']
    
    gen=af_cols[0]
    nextgen=af_cols[1]       

    offs=.01-rep/500

    C_candidates=(boot_var_dict_barghi[str(gen)+'_'+str(nextgen)]/(np.mean(pvals,1)*(1-np.mean(pvals,1))))
    C_mean=np.mean(C_candidates,axis=0)

    ax1.errorbar(np.mean(pvals,1)+offs,C_mean-np.min(C_mean),yerr=[C_mean-np.percentile(C_candidates,q=2.5,axis=0),
                                               np.percentile(C_candidates,q=97.5,axis=0)-C_mean],color=color,linewidth=1.)
                 
    ax1.plot(np.mean(pvals,1)+offs,C_mean-np.min(C_mean),'.',markersize=8,color=color,label=str(rep+1))
    
    ax1.set_ylabel(r'Variance Coefficient $C_t$',fontsize=10)
    
    ax1.set_title(r'Barghi et al. (1st iteration)',fontsize=9)
    #ax1.set_ylim([-0.0015,0.007])
    
    #print(np.min(C_mean))
    
    ax1.annotate(r'A',[0.95,0.925],xycoords='axes fraction',fontsize=10)
    
    

ax1.legend(title="Replicate",loc="upper center", columnspacing=0, labelspacing=0, 
          handletextpad=0, ncol=5,title_fontsize=9, fontsize=9, borderaxespad=0, frameon=False)
    
ax2=fig.add_subplot(gs[1])
ax2.ticklabel_format(axis='y',style='sci',scilimits=(0,0))

for rep in range(3):
    
    af_cols=[['pRA0','pRA7'],['pRB0','pRB7'],['pRC0','pRC7']][rep]
    
    color = next(ax2._get_lines.prop_cycler)['color']
    
    gen=af_cols[0]
    nextgen=af_cols[1]       

    offs=.01-rep/500
    
    C_candidates=boot_var_dict_kelly[str(gen)+'_'+str(nextgen)]/(np.mean(pvals,1)*(1-np.mean(pvals,1)))
    C_mean=np.mean(C_candidates,axis=0)
            
    ax2.errorbar(np.mean(pvals,1)+offs,C_mean-np.min(C_mean),yerr=[C_mean-np.percentile(C_candidates,q=2.5,axis=0),
                                               np.percentile(C_candidates,q=97.5,axis=0)-C_mean],color=color,linewidth=1.)
                 
    ax2.plot(np.mean(pvals,1)+offs,C_mean-np.min(C_mean),'.',markersize=8,color=color,label=['A','B','C'][rep])
    
    
    ax2.set_xlabel(r'Major Allele Frequency',fontsize=10)
    ax2.set_ylabel(r'Variance Coefficient $C_t$',fontsize=10)
      
    ax2.set_title(r'Kelly & Hughes',fontsize=10)
    ax2.set_ylim([-0.0015,0.007])
    
    ax2.annotate('B',[0.95,0.925],xycoords='axes fraction',fontsize=10)
 

ax2.legend(title="Replicate",loc="upper center", columnspacing=0, labelspacing=0, 
           handletextpad=0, ncol=5,title_fontsize=9, fontsize=9, borderaxespad=0, frameon=False)
    

#plt.savefig("Cvsp.pdf",bbox_inches='tight')
    
#%%
#Barghi et al. and Bergland et al. C difference with 95% block bootstrap ci

delt=0.025
pvals=[[p, p+delt] for p in np.arange(0.5,1,delt)]

fig = plt.figure(figsize=[3.,5.],constrained_layout=True,dpi=300)
gs = gridspec.GridSpec(nrows=2, ncols=1, figure=fig)

ax1=fig.add_subplot((gs[0]))

ax1.ticklabel_format(axis='y',style='sci',scilimits=(0,0)) 

rep=0
af_cols=['R'+str(rep+1)+'_F'+str(gen) for gen in range(0,61,10)]

for ind,gen in enumerate(af_cols[:-1]):
    color = next(ax1._get_lines.prop_cycler)['color']
    Cdiff_bootstrap=[]        
    for nextgen in af_cols[ind+1:]:

        C_candidates=boot_var_dict_barghi[str(gen)+'_'+str(nextgen)]/(np.mean(pvals,1)*(1-np.mean(pvals,1)))
        Cdiff_bootstrap.append(C_candidates[:,0]-C_candidates[:,-4])
        
    Cdiff_mean=np.mean(Cdiff_bootstrap,1)    
    offs=0.5-ind/6
    
    
    ax1.errorbar(np.arange((ind+1)*10,61,10)+offs,Cdiff_mean,yerr=[Cdiff_mean-np.percentile(Cdiff_bootstrap,q=2.5,axis=1),
                                               np.percentile(Cdiff_bootstrap,q=97.5,axis=1)-Cdiff_mean],color=color,linewidth=1.)
                 
    ax1.plot(np.arange((ind+1)*10,61,10)+offs,Cdiff_mean,'.',markersize=8,color=color,label=str(ind*10))
    
    ax1.set_xticks(range(10,61,10))
    
    ax1.set_xlabel(r'Final Generation',fontsize=10)
    ax1.set_ylabel(r'Excess Variance',fontsize=10)
    
    ax1.set_title(r'Barghi et al. (Rep. 1)',fontsize=10)
    
    ax1.annotate('A',[0.95,0.925],xycoords='axes fraction',fontsize=10)


ax1.legend(title="Initial Generation",loc="upper left", columnspacing=0, labelspacing=0, 
           handletextpad=0, ncol=3,title_fontsize=9, fontsize=9, borderaxespad=0, frameon=False)


delt=0.05
pvals=[[p, p+delt] for p in np.arange(0.5,1,delt)]

ax2=fig.add_subplot((gs[1]))

ax2.ticklabel_format(axis='y',style='sci',scilimits=(0,0)) 

af_cols=['PA_7_2009', 'PA_11_2009', 'PA_7_2010', 'PA_11_2010', 'PA_7_2011',
       'PA_10_2011', 'PA_11_2011']

season_labels=['S09','F09','S10','F10','S11','F11','LF11']

for ind,gen in enumerate(af_cols[:-1]):
    color = next(ax2._get_lines.prop_cycler)['color']
    Cdiff_bootstrap=[]        
    for nextgen in af_cols[ind+1:]:

        C_candidates=boot_var_dict_bergland[str(gen)+'_'+str(nextgen)]/(np.mean(pvals,1)*(1-np.mean(pvals,1)))
        Cdiff_bootstrap.append(C_candidates[:,0]-C_candidates[:,-5])
        
    Cdiff_mean=np.mean(Cdiff_bootstrap,1)    
    offs=0.5-ind/6
    
    
    ax2.errorbar(np.arange((ind+1)*10,61,10)+offs,Cdiff_mean,yerr=[Cdiff_mean-np.percentile(Cdiff_bootstrap,q=2.5,axis=1),
                                               np.percentile(Cdiff_bootstrap,q=97.5,axis=1)-Cdiff_mean],color=color,linewidth=1.)
                 
    ax2.plot(np.arange((ind+1)*10,61,10)+offs,Cdiff_mean,'.',markersize=8,color=color,label=season_labels[ind])
    
    
    ax2.set_xlabel(r'Final Generation',fontsize=10)
    ax2.set_ylabel(r'Excess Variance',fontsize=10)
    
    ax2.set_xticks(range(10,61,10))
    ax2.set_xticklabels(season_labels[1:])
    
    ax2.set_title(r'Bergland et al. (PA)',fontsize=10)
    
    ax2.set_ylim([0,1.2e-2])
    
    ax2.annotate('B',[0.95,0.925],xycoords='axes fraction',fontsize=10)
    
ax2.legend(title="Initial Generation",loc="upper left", columnspacing=0, labelspacing=0, 
           handletextpad=0, ncol=3,title_fontsize=9, fontsize=9, borderaxespad=0, frameon=False)

#plt.savefig("Cdiff_cumulative.eps",bbox_inches='tight')


#%%
#Barghi et al. and Kelly et al. sigma(s|p) lower bound with 95% block bootstrap ci

delt=0.025
pvals=[[p, p+delt] for p in np.arange(0.5,1,delt)]

fig = plt.figure(figsize=[3.,5.],constrained_layout=True,dpi=300)
gs = gridspec.GridSpec(nrows=3, ncols=1, figure=fig)

ax1=fig.add_subplot((gs[0]))

ax1.ticklabel_format(axis='y',style='sci',scilimits=(0,0)) 

t=10

for rep in range(10):
    color = next(ax1._get_lines.prop_cycler)['color']
    
    af_cols=['R'+str(rep+1)+'_F'+str(gen) for gen in range(0,61,10)]
    
    gen=af_cols[0]
    nextgen=af_cols[1]
         
    C_candidates=boot_var_dict_barghi[str(gen)+'_'+str(nextgen)]/(np.mean(pvals,1)*(1-np.mean(pvals,1)))
    Cdiff_bootstrap=(C_candidates[:,0]-C_candidates[:,-4])/(t**2*(np.mean(pvals[0])*(1-np.mean(pvals[0]))))
    
    
    Cdiff_mean=np.mean(Cdiff_bootstrap)
    
    
    ax1.errorbar(rep+1,Cdiff_mean,yerr=[[Cdiff_mean-np.percentile(Cdiff_bootstrap,q=2.5)],
                                               [np.percentile(Cdiff_bootstrap,q=97.5)-Cdiff_mean]],color=color,linewidth=1.)
                 
    ax1.plot(rep+1,Cdiff_mean,'.',markersize=8,color=color)
       
    ax1.set_xticks(range(1,11))
    
    ax1.set_title(r'Barghi et al. (1st iteration)',pad=-10,y=1,fontsize=9)
    
    ax1.set_ylim([0,2e-4])
    
    ax1.annotate('A',[0.95,0.85],xycoords='axes fraction',fontsize=10)
    
    ax1.set_xlabel(r'Replicate',fontsize=9)



ax2=fig.add_subplot(gs[1])
ax2.ticklabel_format(axis='y',style='sci',scilimits=(0,0))

t=15

for rep in range(3):
    
    af_cols=[['pRA0','pRA7'],['pRB0','pRB7'],['pRC0','pRC7']][rep]
    
    color = next(ax2._get_lines.prop_cycler)['color']
    
    gen=af_cols[0]
    nextgen=af_cols[1]       

    C_candidates=boot_var_dict_kelly[str(gen)+'_'+str(nextgen)]/(np.mean(pvals,1)*(1-np.mean(pvals,1)))
    Cdiff_bootstrap=(C_candidates[:,0]-C_candidates[:,-4])/(t**2*(np.mean(pvals[0])*(1-np.mean(pvals[0]))))
    
    Cdiff_mean=np.mean(Cdiff_bootstrap)
    
    
    ax2.errorbar(rep+1,Cdiff_mean,yerr=[[Cdiff_mean-np.percentile(Cdiff_bootstrap,q=2.5)],
                                               [np.percentile(Cdiff_bootstrap,q=97.5)-Cdiff_mean]],color=color,linewidth=1.)
                 
    ax2.plot(rep+1,Cdiff_mean,'.',markersize=8,color=color)
       
    ax2.set_xticks(range(1,4))
                 
    ax2.set_title(r'Kelly & Hughes',pad=-10,y=1,fontsize=9)
    
    ax2.annotate('B',[0.95,0.85],xycoords='axes fraction',fontsize=10)
    
    ax2.set_ylabel(r'Among-Locus Selection Coefficient Variance $\sigma^2(s|p)$',fontsize=10)
    
    ax2.set_xlabel(r'Replicate',fontsize=9)
    
    ax2.set_xlim([0.5,3.5])



delt=0.05
pvals=[[p, p+delt] for p in np.arange(0.5,1,delt)]

ax3=fig.add_subplot((gs[2]))


ax3.ticklabel_format(axis='y',style='sci',scilimits=(0,0))

af_cols=['PA_7_2009', 'PA_11_2009', 'PA_7_2010', 'PA_11_2010', 'PA_7_2011',
        'PA_10_2011', 'PA_11_2011']

season_labels=['S09','F09','S10','F10','S11','F11','LF11']

iteration_labels=[season_labels[ind]+'/'+season_labels[ind+1] for ind in range(len(season_labels)-1)]

for ind,gen in enumerate(af_cols[:-2]):
    color = next(ax2._get_lines.prop_cycler)['color']
        
    nextgen=af_cols[ind+1]


    if season_labels[ind][0]=='S':
        t=10
    else:
        #Number of winter generations is unknown
        #Bergland et al. suggest t=2
        #Creates a weird looking plot so leaving t same for ease of comparison
        t=10
        
    C_candidates=boot_var_dict_bergland[str(gen)+'_'+str(nextgen)]/(np.mean(pvals,1)*(1-np.mean(pvals,1)))
    Cdiff_bootstrap=(C_candidates[:,0]-C_candidates[:,-5])/(t**2*(np.mean(pvals[0])*(1-np.mean(pvals[0]))))
        
    Cdiff_mean=np.mean(Cdiff_bootstrap)
    
      
    ax3.errorbar(ind+1,Cdiff_mean,yerr=[[Cdiff_mean-np.percentile(Cdiff_bootstrap,q=2.5)],
                                               [np.percentile(Cdiff_bootstrap,q=97.5)-Cdiff_mean]],color=color,linewidth=1.)
                 
    ax3.plot(ind+1,Cdiff_mean,'.',markersize=8,color=color)



ax3.set_title(r'Bergland et al. (PA)',pad=-10,y=1,fontsize=9)

ax3.set_ylim([-1e-5,5e-4])

ax3.set_xticks(range(1,6))
ax3.set_xticklabels(iteration_labels[0:-1],rotation=45,fontsize=8)
ax3.set_xlabel("Iteration",fontsize=9)
#ax3.set_xlim([0,6])


ax3.annotate('C',[0.95,0.85],xycoords='axes fraction',fontsize=10)
    

#plt.savefig('C:\\Users\\jason\\Documents\\GitHub\\polygenic_variance\\sigma2s.eps',bbox_inches='tight')


#%%
#simulations

Cdiff_mean_all_reps={}

pvals=[[0.5,0.6],[0.875,0.925]]

for directory in os.listdir('D:\\slim\\'):
    path ='D:\\slim\\'+ directory+'\\boot\\'

    boot_var_dict={}
    boot_mean_dict={}
    for filename in os.listdir(path):
        if filename[:3]=='var':
            boot_var_dict[filename[3:-4]]=np.loadtxt(path+filename)
        elif filename[:4]=='mean':
            boot_mean_dict[filename[4:-4]]=np.loadtxt(path+filename)


    
    seeds=set([_[3:16] for _ in os.listdir(path) if os.path.isfile(os.path.join(path,_)) and _[:3]=='var'])
    
    # plt.figure()
    
    Cdiff_mean_all=[]
    for seed in seeds:
            
        #plt.figure()
        ax=plt.gca()
        
        gen='10000'
        nextgen='10010'
        
        Cdiff_bootstrap=[]        
       
        C_candidates=boot_var_dict[seed+'_'+gen+'_'+nextgen]/(np.mean(pvals,1)*(1-np.mean(pvals,1)))
        Cdiff_bootstrap.append(C_candidates[:,0]-C_candidates[:,-1])
        
        Cdiff_mean=np.mean(Cdiff_bootstrap,1)
        
        Cdiff_mean_all.append(Cdiff_mean[0])
        
        
        #plt.errorbar([1],Cdiff_mean,
        #            yerr=[Cdiff_mean-np.percentile(Cdiff_bootstrap,q=2.5,axis=1),np.percentile(Cdiff_bootstrap,q=97.5,axis=1)-Cdiff_mean]
        #            ,color=color)
        
        #plt.plot([1],Cdiff_mean,'.',markersize=15,color=color,label=str(gen))
        
        #plt.xlabel('Generation',fontsize=14)
        #plt.ylabel(r'Excess $\Delta p$ variance',fontsize=14)
        #plt.ylim([-0.004,0.004])
        #plt.legend(title="Reference generation",loc="upper left", ncol=2)    
        #plt.savefig("Cdiff_"+str(rep)+".pdf")
    
    Cdiff_mean_all_reps[directory]=Cdiff_mean_all
    
    # plt.hist(Cdiff_mean_all,label=directory)
    # plt.legend()
    # plt.xlim([-0.02,0.02])



#%%
fig = plt.figure(figsize=[3.,4],constrained_layout=True,dpi=300)
gs = gridspec.GridSpec(nrows=2, ncols=2, figure=fig)


ax2=fig.add_subplot((gs[0,:]))

ax2.ticklabel_format(axis='y',style='sci',scilimits=(0,0))

cols=['background_s5em2_u1em8',
      'background_s5em2_u1em9',
      'background_s1em2_u1em8',
      'background_s1em2_u1em9',
      'neutral',
      'positive_s1em2_u1em10',
      'positive_s1em2_u1em9',
      'positive_s2em2_u1em10',
      'positive_s2em2_u1em9']

Cdiff_df=pd.DataFrame(Cdiff_mean_all_reps)
Cdiff_df=Cdiff_df[cols]
Cdiff_array=np.array(Cdiff_df)

medianprops=dict(color='k')

bplot=ax2.boxplot(Cdiff_array,showfliers=False,patch_artist=True,notch=True,
            medianprops=medianprops)
ax2.set_ylim([-0.01,0.02])

strong_backroundsel_color=next(ax2._get_lines.prop_cycler)['color']
bplot['boxes'][0].set_facecolor(strong_backroundsel_color)

for patch in bplot['boxes'][1:-1]:
    patch.set_facecolor(next(ax2._get_lines.prop_cycler)['color'])
    
strong_positivesel_color=next(ax2._get_lines.prop_cycler)['color']
bplot['boxes'][-1].set_facecolor(strong_positivesel_color)

ax2.fill_between([0.5,4.5],[-0.01,-0.01],[0.02,0.02],alpha=0.3)
ax2.fill_between([4.5,5.5],[-0.01,-0.01],[0.02,0.02],alpha=0.3)
ax2.fill_between([5.5,9.5],[-0.01,-0.01],[0.02,0.02],alpha=0.3)


ax2.annotate(r'$U:\;\: 1\quad\: 0.1\quad\: 1\quad\: 0.1 \qquad\:\; 0.01 \:\: 0.1 \:\: 0.01 \:\: 0.1$',
             [-0.05,-0.12],xycoords='axes fraction',fontsize=8)
ax2.annotate(r'$s:\: -0.05 \quad\: -0.01 \qquad\qquad\:\: 0.01\qquad\quad 0.02$',
             [-0.05,-0.25],xycoords='axes fraction',fontsize=8)

ax2.annotate(r'Negative',[0.05,0.9],xycoords='axes fraction',fontsize=9)
ax2.annotate(r'Positive',[0.65,0.9],xycoords='axes fraction',fontsize=9)
ax2.annotate(r'Neutral',[0.47,0.6],xycoords='axes fraction',fontsize=9,rotation=90)
ax2.annotate('A',[0.94,0.05],xycoords='axes fraction',fontsize=10)
ax2.annotate(r"*",[0.8 , 1.2e-2], xycoords='data',fontsize=14)
ax2.annotate(r"*",[8.55 , 1.2e-2], xycoords='data',fontsize=14)
ax2.annotate(r"*",[9.05 , 1.2e-2], xycoords='data',fontsize=14)

ax2.set_xticklabels([' ',' ',' ',' ',' ',' ',' ',' ',' '],fontsize=20)

ax2.set_ylabel(r'Excess Variance',fontsize=10)

ax2.set_title(r'Simulation',fontsize=10)



def var_s_p(svals,pvals):
    return np.array([np.var([_[2] for _ in svals if ((_[0]<=p[1]) & (_[0]>=p[0])) | ((_[0]>=1-p[1]) & (_[0]<=1-p[0]))]) for p in pvals])



ax0=fig.add_subplot((gs[1,0]))

ax0.ticklabel_format(axis='y',style='sci',scilimits=(0,0))

path='D:\\slim\\background_s5em2_u1em8\\'

with open(path+'total_s.p', 'rb') as file:
    total_s_back = pickle.load(file)

delt=0.05
pvals=[[p, p+delt] for p in np.arange(0.5,1,delt)]
var_s_p_all_back=np.array([var_s_p(total_s_back[seed],pvals) for seed in total_s_back.keys()])

ps=np.array([p+delt/2 for p in np.arange(0.5,1,delt)])
bplot=ax0.boxplot(var_s_p_all_back*ps*(1-ps),widths=0.03,
                  showfliers=False,notch=True,medianprops=medianprops,
                  patch_artist=True,positions=ps)
ax0.set_xlim([0.5,0.95])
ax0.set_ylim([0.,0.045/100])
ax0.set_xticklabels(['0.5','','','','','','','','0.9',''])
ax0.set_xlabel(r'Major Allele Frequency',fontsize=10,position=[1.1,0.5])
ax0.set_ylabel('Selective divergence \n $p(1-p)\sigma^2(\overline{s}|p)$')
ax0.annotate('B',[0.85,0.9],xycoords='axes fraction',fontsize=10)
ax0.annotate(r"*",[.8 , 1.0], xycoords='axes fraction',fontsize=14)

for patch in bplot['boxes']:
    patch.set_facecolor(strong_backroundsel_color)
    

ax1=fig.add_subplot((gs[1,1]))

ax1.ticklabel_format(axis='y',style='sci',scilimits=(0,0))

path='D:\\slim\\positive_s2em2_u1em9\\'

with open(path+'total_s.p', 'rb') as file:
    total_s_pos = pickle.load(file)

delt=0.05
pvals=[[p, p+delt] for p in np.arange(0.5,1,delt)]
var_s_p_all_pos=np.array([var_s_p(total_s_pos[seed],pvals) for seed in total_s_pos.keys()])

bplot=ax1.boxplot(var_s_p_all_pos*ps*(1-ps),widths=0.03,
                  showfliers=False,notch=True,medianprops=medianprops,
                  patch_artist=True,positions=ps)
ax1.set_xlim([0.5,0.95])
ax1.set_ylim([0.,0.045/100])
ax1.set_xticklabels(['0.5','','','','','','','','0.9',''])
ax1.annotate('C',[0.85,0.9],xycoords='axes fraction',fontsize=10)
ax1.annotate(r"**",[.75 , 1.0], xycoords='axes fraction',fontsize=14)


for patch in bplot['boxes']:
    patch.set_facecolor(strong_positivesel_color)


plt.savefig('C:\\Users\\jason\\Documents\\GitHub\\polygenic_variance\\simulations.pdf',bbox_inches='tight')


#%%
#Supplemental figures

#%%
#Barghi et al. C difference with 95% block bootstrap ci
#All replicates


delt=0.025
pvals=[[p, p+delt] for p in np.arange(0.5,1,delt)]

fig = plt.figure(figsize=[7.,3.],constrained_layout=True,dpi=600)
gs = gridspec.GridSpec(nrows=2, ncols=5, figure=fig)


for rep in range(10):

    ax1=fig.add_subplot((gs[rep]))
    
    ax1.ticklabel_format(axis='y',style='sci',scilimits=(0,0)) 
    
    af_cols=['R'+str(rep+1)+'_F'+str(gen) for gen in range(0,61,10)]
    
    for ind,gen in enumerate(af_cols[:-1]):
        color = next(ax1._get_lines.prop_cycler)['color']
        Cdiff_bootstrap=[]        
        for nextgen in af_cols[ind+1:]:
    
            C_candidates=boot_var_dict_barghi[str(gen)+'_'+str(nextgen)]/(np.mean(pvals,1)*(1-np.mean(pvals,1)))
            Cdiff_bootstrap.append(C_candidates[:,0]-C_candidates[:,-4])
            
        Cdiff_mean=np.mean(Cdiff_bootstrap,1)    
        offs=0.5-ind/6
        
        
        ax1.errorbar(np.arange((ind+1)*10,61,10)+offs,Cdiff_mean,yerr=[Cdiff_mean-np.percentile(Cdiff_bootstrap,q=2.5,axis=1),
                                                   np.percentile(Cdiff_bootstrap,q=97.5,axis=1)-Cdiff_mean],color=color,linewidth=1.)
                     
        ax1.plot(np.arange((ind+1)*10,61,10)+offs,Cdiff_mean,'.',markersize=8,color=color,label=str(ind*10))
        
        ax1.set_xticks(range(10,61,10))
        
        ax1.set_xticklabels([10,'','','','',60])
        
        if rep==7:
            ax1.set_xlabel(r'Final Generation',fontsize=10)
        
        if rep==0:
            ax1.set_ylabel(r' ',loc='bottom',fontsize=10)
        
        ax1.text(0.05,0.85,'Rep. '+str(rep+1),transform=ax1.transAxes,fontsize=9)
        
        ax1.set_ylim([-1e-3,1.8e-2])

fig.text(0,0.35,'Excess Variance',rotation=90)

plt.savefig("Cdiff_barghi_all.eps",bbox_inches='tight')



#%%


fig = plt.figure(figsize=[4.,3],constrained_layout=True,dpi=300)
ax1=fig.add_subplot()

delt=0.05
pvals=[[p, p+delt] for p in np.arange(0.5,1,delt)]

#ax2=fig.add_subplot((gs[1]))

#ax2.ticklabel_format(axis='y',style='sci',scilimits=(0,0)) 

af_cols=['PA_7_2009', 'PA_11_2009', 'PA_7_2010', 'PA_11_2010', 'PA_7_2011',
       'PA_10_2011', 'PA_11_2011']

season_labels=['S09','F09','S10','F10','S11','F11','LF11']

for ind,gen in enumerate(af_cols[:-1]):
    color = next(ax2._get_lines.prop_cycler)['color']
    Cdiff_bootstrap=[]        
    #for nextgen in af_cols[ind+1:]:
    nextgen=af_cols[ind+1]

    C_candidates=boot_var_dict_bergland[str(gen)+'_'+str(nextgen)]/(np.mean(pvals,1)*(1-np.mean(pvals,1)))
    #Cdiff_bootstrap.append(C_candidates[:,0]-C_candidates[:,-5])
        
    C_mean=np.mean(C_candidates,axis=0)  
    offs=0.01-ind/500
    
    
    ax1.errorbar(np.mean(pvals,1)+offs,C_mean-np.min(C_mean),yerr=[C_mean-np.percentile(C_candidates,q=2.5,axis=0),
                                               np.percentile(C_candidates,q=97.5,axis=0)-C_mean],color=color,linewidth=1.)
                 
    ax1.plot(np.mean(pvals,1)+offs,C_mean-np.min(C_mean),'.',markersize=8,color=color,label=season_labels[ind]+r'$\rightarrow$'+season_labels[ind+1])
    
    ax1.set_ylabel(r'Variance Coefficient $C_t$',fontsize=10)
    
    ax1.set_xlim([0.5,0.9])
    ax1.set_ylim([-0.0015,0.012])
    
    
    
ax1.legend(loc="upper left", columnspacing=0, labelspacing=0, 
           handletextpad=0, ncol=3,title_fontsize=9, fontsize=8, borderaxespad=0, frameon=False)


ax1.set_xlabel(r'Major Allele Frequency',fontsize=10)

plt.savefig("Cvsp_bergland.pdf",bbox_inches='tight')

# #%%
# #Barghi mean delta p

# delt=0.025
# pvals=[[p, p+delt] for p in np.arange(0.5,1,delt)]

# fig = plt.figure(figsize=[3.,3.],constrained_layout=True,dpi=300)
# ax1=fig.add_subplot()

# ax1.ticklabel_format(axis='y',style='sci',scilimits=(0,0))

# for rep in range(10):
#     af_cols=['R'+str(rep+1)+'_F'+str(gen) for gen in range(0,61,10)]
       
    
#     gen=af_cols[0]
#     nextgen=af_cols[1]       

#     means=boot_mean_dict_barghi[str(gen)+'_'+str(nextgen)]#/(np.mean(pvals,1)*(1-np.mean(pvals,1)))
    
#     ax1.fill_between(np.mean(pvals,1),np.percentile(means,q=5,axis=0),
#                      np.percentile(means,q=95,axis=0),color='C2',alpha=0.5)
                
#     ax1.plot(np.mean(pvals,1),np.mean(means,0),label=str(rep+1))
    
#     ax1.plot([0.5,1],[0,0],'k')
    
#     ax1.set_xlabel(r'Major Allele Frequency',fontsize=10)
#     ax1.set_ylabel(r'Mean AF change $E[\Delta_t p | p]$',fontsize=10)
    
#     ax1.set_title(r'Barghi et al.',fontsize=10)
    
#     ax1.set_ylim([-0.02,0.02])

# ax1.legend(title="Replicate",loc="upper center", columnspacing=0.5, labelspacing=0, 
#            handletextpad=0.1, ncol=5,title_fontsize=9, fontsize=9, borderaxespad=0, frameon=False)
            
# ax1.savefig("mean_deltp_barghi.pdf",bbox_inches='tight')

# #%%
# #Kelly mean delta p

# delt=0.025
# pvals=[[p, p+delt] for p in np.arange(0.5,1,delt)]

# fig = plt.figure(figsize=[3.,3.],constrained_layout=True,dpi=300)
# ax1=fig.add_subplot()

# ax1.ticklabel_format(axis='y',style='sci',scilimits=(0,0))

# for rep in range(3):
    
#     af_cols=[['pRA0','pRA7'],['pRB0','pRB7'],['pRC0','pRC7']][rep]
    
#     for ind,gen in enumerate(af_cols[:-1]):    
#         nextgen=af_cols[ind+1]
    
#         meandeltp=boot_mean_dict_kelly[str(gen)+'_'+str(nextgen)]
        
        
#         ax1.fill_between(np.mean(pvals,1),np.percentile(meandeltp,q=5,axis=0),
#                          np.percentile(meandeltp,q=95,axis=0),color='C2',alpha=0.5)
                    
#         ax1.plot(np.mean(pvals,1),np.mean(meandeltp,0),label=['A','B','C'][rep])
        
#         ax1.plot([0.5,1],[0,0],'k')
        
#         ax1.set_xlabel(r'Major Allele Frequency',fontsize=10)
#         ax1.set_ylabel(r'Mean AF change $E[\Delta_t p | p]$',fontsize=10)
      
#         ax1.set_title(r'Kelly & Hughes',fontsize=10)
        
#         ax1.set_ylim([-0.02,0.02])
        
# ax1.legend(title="Replicate",loc="upper right", title_fontsize=9, fontsize=9, borderaxespad=0, frameon=False)

# plt.savefig("mean_deltp_kelly.pdf",bbox_inches='tight')

# #%%
# #Bergland mean delta p

# delt=0.05
# pvals=[[p, p+delt] for p in np.arange(0.5,1,delt)]

# af_cols=['PA_7_2009', 'PA_11_2009', 'PA_7_2010', 'PA_11_2010', 'PA_7_2011',
#        'PA_10_2011', 'PA_11_2011']

# season_labels=['S09','F09','S10','F10','S11','F11','LF11']

# fig = plt.figure(figsize=[3.,3.],constrained_layout=True,dpi=300)
# ax1=fig.add_subplot()

# ax1.ticklabel_format(axis='y',style='sci',scilimits=(0,0))

# for ind,gen in enumerate(af_cols[:-1]):
#     nextgen=af_cols[ind+1]
    
#     meandeltp=boot_mean_dict_bergland[str(gen)+'_'+str(nextgen)]
    
#     ax1.fill_between(np.mean(pvals,1),np.percentile(meandeltp,q=5,axis=0),
#                          np.percentile(meandeltp,q=95,axis=0),color='C2',alpha=0.5)
                    
#     ax1.plot(np.mean(pvals,1),np.mean(meandeltp,0),label=season_labels[ind]+r'$\rightarrow$'+season_labels[ind+1])
    
#     ax1.plot([0.5,1],[0,0],'k')
    
#     ax1.set_xlabel(r'Major Allele Frequency',fontsize=10)
#     ax1.set_ylabel(r'Mean AF change $E[\Delta_t p | p]$',fontsize=10)
  
#     ax1.set_title(r'Bergland et al. (PA)',fontsize=10)
    
#     ax1.set_ylim([-0.1,0.05])
        
# ax1.legend(loc="upper right", title_fontsize=9, fontsize=9, borderaxespad=0, frameon=False, ncol=2,
#            columnspacing=1, labelspacing=0)

# plt.savefig("mean_deltp_bergland.pdf",bbox_inches='tight')
