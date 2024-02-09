#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 11:07:49 2023

@author: joannafernandez
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#import pysam
#from pysam import FastaFile


def Collect_data(file1,fact):
    full = pd.read_csv(file1, sep="	")
    #reverse = pd.read_csv(file2)
    full['score'] = full['Watson'] + full['Crick']
  #  print(full)
    print (full['score'].sum())
    full['Wnorm']= (full['Watson'] / full['score'].sum())*fact
    full['Cnorm']= (full['Crick'] / full['score'].sum())*fact
    full['Cnorm']= full['Cnorm']*-1
    full['smoo_wnorm'] = full['Wnorm'].rolling(window = 10, center=True).mean()
    full['smoo_cnorm'] = full['Cnorm'].rolling(window = 10, center=True).mean()
    full['normscore']= (full['score']/ full['score'].sum())*fact
  #  print(full)
    chr1= full[(full["Chr"] == 1)]
    chr2= full[(full["Chr"] == 2)]
    chr3= full[(full["Chr"] == 3)]
    return chr1,chr2,chr3,full

wtchr1, wtchr2, wtchr3, wt = Collect_data('FullMap.pombase220208-23-cc-316-w-thiamine.txt', 10000)
sen1chr1, sen1chr2, sen1chr3, sen1 = Collect_data('FullMap.pombase220208-8-cc-317-w-thiamine.txt', 10000)
dbl8chr1, dbl8chr2, dbl8chr3, dbl8 = Collect_data('FullMap.pombase220208-9-cc-318-w-thiamine.txt', 10000)
dschr1, dschr2, dschr3, ds = Collect_data('FullMap.pombase220208-35-cc-319-w-thiamine.txt', 10000)

wtchr1, wtchr2, wtchr3, wt = Collect_data('FullMap.pombase220208-23-cc-316-w-thiamine_NOT3bpOffset_with_strandminTH_0.txt', 10000)
sen1chr1, sen1chr2, sen1chr3, sen1 = Collect_data('FullMap.pombase220208-8-cc-317-w-thiamine_NOT3bpOffset_with_strandminTH_0.txt', 10000)
dbl8chr1, dbl8chr2, dbl8chr3, dbl8 = Collect_data('FullMap.pombase220208-9-cc-318-w-thiamine_NOT3bpOffset_with_strandminTH_0.txt', 10000)
dschr1, dschr2, dschr3, ds = Collect_data('FullMap.pombase220208-35-cc-319-w-thiamine_NOT3bpOffset_with_strandminTH_0.txt', 10000)

wtchr1, wtchr2, wtchr3, wt = Collect_data('FullMap.pombase220208-23-cc-316-w-thiamine_3bpOffset_with_strandminTH_0.txt', 10000)
sen1chr1, sen1chr2, sen1chr3, sen1 = Collect_data('FullMap.pombase220208-8-cc-317-w-thiamine_3bpOffset_with_strandminTH_0.txt', 10000)
dbl8chr1, dbl8chr2, dbl8chr3, dbl8 = Collect_data('FullMap.pombase220208-9-cc-318-w-thiamine_3bpOffset_with_strandminTH_0.txt', 10000)
dschr1, dschr2, dschr3, ds = Collect_data('FullMap.pombase220208-35-cc-319-w-thiamine_3bpOffset_with_strandminTH_0.txt', 10000)

#new ethanol prefix

t2rchr1, t2rchr2, t2rchr3, t2r = Collect_data('FullMap.pombase220208-T2R.txt', 10000)
t2echr1, t2echr2, t2echr3, t2e = Collect_data('FullMap.pombase220208-T2E.txt', 10000)


#dup free
wtchr1, wtchr2, wtchr3, wt = Collect_data('FullMap.pombase220208-23-cc-316-w-thiamine_dupfree.txt', 10000)
sen1chr1, sen1chr2, sen1chr3, sen1 = Collect_data('FullMap.pombase220208-8-cc-317-w-thiamine_dupfree.txt', 10000)
dbl8chr1, dbl8chr2, dbl8chr3, dbl8 = Collect_data('FullMap.pombase220208-9-cc-318-w-thiamine_dupfree.txt', 10000)
dschr1, dschr2, dschr3, ds = Collect_data('FullMap.pombase220208-35-cc-319-w-thiamine_dupfree.txt', 10000)





wtchr1b, wtchr2b, wtchr3b, wtb = Collect_data('FullMap.pombase220208-WG141_WT6mins_0.1Kbp_binned_OS3_full_strand.txt', 10000)
t2tchr1b, t2tchr2b, t2tchr3b, t2tb = Collect_data('FullMap.pombase220208-WG142_top2-191t_0.1Kbp_binned_OS3_full_strand.txt', 10000)
t2mchr1b, t2mchr2b, t2mchr3b, t2mb = Collect_data('FullMap.pombase220208-WG143_top2-191m_0.1Kbp_binned_OS3_full_strand.txt', 10000)

wtchr1, wtchr2, wtchr3, wt = Collect_data('FullMap.pombase220208-WG141_WT6mins.txt', 10000)
t2tchr1, t2tchr2, t2tchr3, t2t = Collect_data('FullMap.pombase220208-WG142_top2-191t.txt',10000)
t2mchr1, t2mchr2, t2mchr3, t2m = Collect_data('FullMap.pombase220208-WG143_top2-191m.txt',10000)


#new senataxin cc-ethanol
wtchr1, wtchr2, wtchr3, wt = Collect_data('FullMap.pombase220208-316_cc_23.txt', 10000)
sen1chr1, sen1chr2, sen1chr3, sen1 = Collect_data('FullMap.pombase220208-317_cc_23.txt', 10000)
dbl8chr1, dbl8chr2, dbl8chr3, dbl8 = Collect_data('FullMap.pombase220208-318_cc_23.txt', 10000)
dschr1, dschr2, dschr3, ds = Collect_data('FullMap.pombase220208-319_cc_23.txt', 10000)

#dup free senataxin cc- ethanol 
wtchr1, wtchr2, wtchr3, wt = Collect_data('FullMap.pombase220208-316_cc_23_dupfree.txt', 10000)
sen1chr1, sen1chr2, sen1chr3, sen1 = Collect_data('FullMap.pombase220208-317_cc_23_dupfree.txt', 10000)
dbl8chr1, dbl8chr2, dbl8chr3, dbl8 = Collect_data('FullMap.pombase220208-318_cc_23_dupfree.txt', 10000)
dschr1, dschr2, dschr3, ds = Collect_data('FullMap.pombase220208-319_cc_23_dupfree.txt', 10000)


#binned 
wtchr1, wtchr2, wtchr3, wt = Collect_data('FullMap.pombase220208-316_cc_23_0.1Kbp_binned_OS0_full_strand.txt', 10000)
sen1chr1, sen1chr2, sen1chr3, sen1 = Collect_data('FullMap.pombase220208-317_cc_23_0.1Kbp_binned_OS0_full_strand.txt', 10000)
dbl8chr1, dbl8chr2, dbl8chr3, dbl8 = Collect_data('FullMap.pombase220208-318_cc_23_0.1Kbp_binned_OS0_full_strand.txt', 10000)
dschr1, dschr2, dschr3, ds = Collect_data('FullMap.pombase220208-319_cc_23_0.1Kbp_binned_OS0_full_strand.txt', 10000)

FullMap.pombase220208-316_cc_23_0.1Kbp_binned_OS0_full_strand




def Collect_data(file1,fact):
    full = pd.read_csv(file1, sep="	")
    #reverse = pd.read_csv(file2)
   # full['score'] = full['Watson'] + full['Crick']
  #  print(full)
   # print (full['score'].sum())

    return full

wt = Collect_data('FullMap.pombase220208-316_cc_23.txt', 10000)
sen1 = Collect_data('FullMap.pombase220208-317_cc_23.txt', 10000)
dbl8 = Collect_data('FullMap.pombase220208-318_cc_23.txt', 10000)
ds = Collect_data('FullMap.pombase220208-319_cc_23.txt', 10000)


#%%%

custom_palette = ["grey", "deepskyblue", "orange", "tomato" ]
sns.set_palette(custom_palette)


wt['source'] = 'WT'
sen1['source'] = 'sen1D'
dbl8['source'] = 'dbl8D'
ds['source'] = 'ds'

# Concatenate the dataframes
combined_averaged = pd.concat([wt, sen1, dbl8, ds], ignore_index=True)
combined_averaged['log2_score'] = np.log2(combined_averaged['normscore'])
import seaborn as sns
plt.figure(figsize=(10, 8))
sns.violinplot(data=combined_averaged, x="source", y="log2_score", hue='source')


normscore_df1 = wt['score']
normscore_df2 = sen1['score']

normscore_df3 = dbl8['score']
normscore_df4 = ds['score']

# Create a new DataFrame with the selected 'normscore' columns
new_df = pd.DataFrame({'wt': normscore_df1, 'sen1D': normscore_df2, 'dbl8D': normscore_df3, 'ds': normscore_df4})

iris_corr_matrix = new_df.corr()
print(iris_corr_matrix)
import seaborn as sns
# Create the heatmap using the `heatmap` function of Seaborn
sns.heatmap(iris_corr_matrix, cmap='coolwarm', annot=True, vmin= 0.60, vmax = 1)


#%%

def Find(file):
    genes = pd.read_csv(file, delimiter="\t")
 #   print(genes)
   # genes['length'] = genes['end'] - genes['start']
  #  print(genes['length'])
    genesfor = genes[(genes['coding_strand'] == 'forward')]
    genesrev = genes[(genes['coding_strand'] == 'reverse')]
    
    genes.loc[genes['chro'] == "chr1", 'chro'] = 1
    genes.loc[genes['chro'] == "chr2", 'chro'] = 2
    genes.loc[genes['chro'] == "chr3", 'chro'] = 3


    genesfor['Sbinpos'] = genesfor['start']/300
    genesfor['Sbinpos'] = genesfor['Sbinpos'].astype(int)
    genesfor['Sbinpos'] = genesfor['Sbinpos']*300 +150
    genesfor['Ebinpos'] = genesfor['end']/300
    genesfor['Ebinpos'] = genesfor['Ebinpos'].astype(int)
    genesfor['Ebinpos'] = genesfor['Ebinpos']*300 +150


    genesrev['Sbinposr'] = genesrev['end']/300
    genesrev['Sbinposr'] = genesrev['Sbinposr'].astype(int)
    genesrev['Sbinposr'] = genesrev['Sbinposr']*300 +150
    genesrev['Ebinposr'] = genesrev['start']/300
    genesrev['Ebinposr'] = genesrev['Ebinposr'].astype(int)
    genesrev['Ebinposr'] = genesrev['Ebinposr']*300 +150

    return genesfor, genesrev, genes

ggenesfor, ggenesrev, gggenes = Find("dbl8_stall_sites_direction.txt")


def Find(file):
    genes = pd.read_csv(file, delimiter=",")
 #   print(genes)
   # genes['length'] = genes['end'] - genes['start']
  #  print(genes['length'])
    genesfor = genes[(genes['coding_strand'] == 'forward')]
    genesrev = genes[(genes['coding_strand'] == 'reverse')]
    
    genes.loc[genes['chro'] == "chr1", 'chro'] = 1
    genes.loc[genes['chro'] == "chr2", 'chro'] = 2
    genes.loc[genes['chro'] == "chr3", 'chro'] = 3


    genesfor['Sbinpos'] = genesfor['start']/300
    genesfor['Sbinpos'] = genesfor['Sbinpos'].astype(int)
    genesfor['Sbinpos'] = genesfor['Sbinpos']*300 +150
    genesfor['Ebinpos'] = genesfor['end']/300
    genesfor['Ebinpos'] = genesfor['Ebinpos'].astype(int)
    genesfor['Ebinpos'] = genesfor['Ebinpos']*300 +150


    genesrev['Sbinposr'] = genesrev['end']/300
    genesrev['Sbinposr'] = genesrev['Sbinposr'].astype(int)
    genesrev['Sbinposr'] = genesrev['Sbinposr']*300 +150
    genesrev['Ebinposr'] = genesrev['start']/300
    genesrev['Ebinposr'] = genesrev['Ebinposr'].astype(int)
    genesrev['Ebinposr'] = genesrev['Ebinposr']*300 +150

    return genesfor, genesrev, genes
controlfor, controlrev, ccontrol = Find('new_control.csv')

sen1stall = gggenes[(gggenes['genotype'] == 'sen1D')]
dbl8stall = gggenes[(gggenes['genotype'] == 'dbl8D')]
doublestall = gggenes[(gggenes['genotype'] == 'sen1dbl8DD_unique')]


sen1stalls = sen1stall['ID'].to_list()
dbl8stalls = dbl8stall['ID'].to_list()
dsstalls = doublestall['ID'].to_list()

def Findfeat(file):
    genes = pd.read_csv(file, delimiter="\t")
    print(genes)
    genes.loc[genes['Chromosome'] == "I", 'chro'] = 'chr1'
    genes.loc[genes['Chromosome'] == "II", 'chro'] = 'chr2'
    genes.loc[genes['Chromosome'] == "III", 'chro'] = 'chr3'
    
    featfor = genes[(genes['Strand'] == 'forward')]
    featrev = genes[(genes['Strand'] == 'reverse')]
    
    featrev.loc[featrev['Start position'] < 25, 'Start position'] = 25
     
    featfor['Sbinpos'] = featfor['Start position']/50
    featfor['Sbinpos'] = featfor['Sbinpos'].astype(int)
    featfor['Sbinpos'] = featfor['Sbinpos']*50 + 25
    featfor['Ebinpos'] = featfor['End position']/50
    featfor['Ebinpos'] = featfor['Ebinpos'].astype(int)
    featfor['Ebinpos'] = featfor['Ebinpos']*50 +25
    
    featrev['Sbinpos'] = featrev['End position']/50 
    featrev['Sbinpos'] = featrev['Sbinpos'].astype(int)
    featrev['Sbinpos'] = featrev['Sbinpos']*50 +25
    featrev['Ebinpos'] = featrev['Start position']/50
    featrev['Ebinpos'] = featrev['Ebinpos'].astype(int)
    featrev['Ebinpos'] = featrev['Ebinpos']*50 +25
    
    return featfor, featrev, genes

featfor, featrev, ffeat = Findfeat('protein_coding_gene_list.tsv')

con1 = ccontrol[(ccontrol['chro'] == 1)]
con2 = ccontrol[(ccontrol['chro'] == 2)]
con3 = ccontrol[(ccontrol['chro'] == 3)]

feat1 = ffeat[(ffeat['chro'] == 'chr1')]
feat2 = ffeat[(ffeat['chro'] == 'chr2')]
feat3 = ffeat[(ffeat['chro'] == 'chr3')]




gene1 = gggenes[(gggenes['chro'] == 1)]
gene2 = gggenes[(gggenes['chro'] == 2)]
gene3 = gggenes[(gggenes['chro'] == 3)]


aaaah = {'start':[3753687,1602264,1070904], 'end': [3789421,1644747,1137003], 'chro': [1,2,3]}
comcentro = pd.DataFrame(data=aaaah)
comcentro['Sbinpos'] = comcentro['start']/50
comcentro['Sbinpos'] = comcentro['Sbinpos'].astype(int)
comcentro['Sbinpos'] = comcentro['Sbinpos']*50 +25
comcentro['Ebinpos'] = comcentro['end']/50
comcentro['Ebinpos'] = comcentro['Ebinpos'].astype(int)
comcentro['Ebinpos'] = comcentro['Ebinpos']*50 +25


d = {'start': [3753687], 'end': [3789421]}
centro1 = pd.DataFrame(data=d)

dd = {'start': [1602264], 'end': [1644747]}
centro2 = pd.DataFrame(data=dd)

ddd = {'start': [1070904], 'end': [1137003]}
centro3 = pd.DataFrame(data=ddd)


te1 = {'start': [1, 5554844], 'end': [29663,5579133]}
telo1 = pd.DataFrame(data=te1)
te2 ={'start': [1, 4500619], 'end': [39186,4539800]}
telo2 = pd.DataFrame(data=te2)
te3 ={'start': [], 'end': []}
telo3 = pd.DataFrame(data=te3)


def Chromosome_plot (cc1, cc2, cc3, cc4, featurex, centro, genee, telo, ccon):
    ff, (ax2, ax3, ax4, ax5, ax6) = plt.subplots(5,1)
   
    #ax1.set_title('WT')
    ax2.plot(cc1['Pos'], cc1['Wnorm'], color ='firebrick', alpha=0.8, linewidth = 1)
    ax2.plot(cc1['Pos'], cc1['Cnorm'], color ='steelblue', alpha=0.8, linewidth = 1)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.set_ylim(-0.2,0.2)
    ax3.set_ylim(-0.2,0.2)
    ax4.set_ylim(-0.2,0.2)
    ax5.set_ylim(-0.2,0.2)

    ax2.set_ylabel('WT (HpM)')
    #ax2.set_ylim(-5,5)
    ax2.set_ylabel('WT')

    ax3.plot(cc2['Pos'], cc2['Wnorm'], color ='firebrick', alpha=0.8, linewidth = 1)
    ax3.plot(cc2['Pos'], cc2['Cnorm'], color ='steelblue', alpha=0.8, linewidth = 1)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.spines['bottom'].set_visible(False)
  #  ax3.set_ylim(-0.5,0.5)
   # ax3.set_ylabel('SEN1d (HpM)')
    ax3.set_ylabel('sen1D (HpM)')

    ax6.set_xlabel('Chromosome position')
    ax4.plot(cc3['Pos'], cc3['Wnorm'], color ='firebrick', alpha=0.8,linewidth = 1)
    ax4.plot(cc3['Pos'], cc3['Cnorm'], color ='steelblue', alpha=0.8, linewidth = 1)
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)
    ax4.spines['bottom'].set_visible(False)
  #  ax4.set_ylim(-0.5,0.5)
    ax4.set_ylabel('dbl8D (HpM)')

    ax5.plot(cc4['Pos'], cc4['Wnorm'], color ='firebrick', alpha=0.8, linewidth = 1)
    ax5.plot(cc4['Pos'], cc4['Cnorm'], color ='steelblue', alpha=0.8, linewidth = 1)
  #  ax5.set_ylim(-0.5,0.5)
    ax5.spines['top'].set_visible(False)
    ax5.spines['right'].set_visible(False)
    ax5.spines['bottom'].set_visible(False)
    ax5.set_ylabel('sen1Ddbl8D (HpM)')
 #   ax4.set_xlabel('Chromosome position')

    ax2.get_xaxis().set_ticks([])
    ax3.get_xaxis().set_ticks([])
    ax4.get_xaxis().set_ticks([])
    ax5.get_xaxis().set_ticks([])
    
    ax2.set_xlim(3800000,3950000)
    ax6.set_xlim(3800000,3950000)
    ax3.set_xlim(3800000,3950000)
    ax4.set_xlim(3800000,3950000)
    ax5.set_xlim(3800000,3950000)


            
                  
    for fe in featurex.itertuples(index=False, name=None):
        if fe[5] == 'reverse':
            if fe[7] == 'sen1D':
                ax6.axvspan(fe[2],fe[3],0.3,0.5,color="deepskyblue",alpha=0.5)
               # ax6.annotate(fe[0],xy=(fe[2],0.3))
            if fe[7] == 'dbl8D':
                ax6.axvspan(fe[2],fe[3],0.3,0.5,color="orange",alpha=0.5)
                #ax6.annotate(fe[0],xy=(fe[2],0.3))
            if fe[7] == 'sen1dbl8DD_unique':
                ax6.axvspan(fe[2],fe[3],0.3,0.5,color="tomato",alpha=0.5)
                #ax6.annotate(fe[0],xy=(fe[2],0.3))
            
        elif fe[5] == 'forward':
            if fe[7] == 'sen1D':
                ax6.axvspan(fe[2],fe[3],0.5,0.8,color="deepskyblue",alpha=0.5)
               # ax6.annotate(fe[0],xy=(fe[3],0.6))
            if fe[7] == 'dbl8D':
                ax6.axvspan(fe[2],fe[3],0.5,0.8,color="orange",alpha=0.5)
               # ax6.annotate(fe[0],xy=(fe[3],0.6))
            if fe[7] == 'sen1dbl8DD_unique':
                ax6.axvspan(fe[2],fe[3],0.5,0.8,color="tomato",alpha=0.5)
               # ax6.annotate(fe[0],xy=(fe[3],0.6))
                ax6.set_ylabel('Gene annotations')
                
    for c in centro.itertuples(index=False, name=None):
            ax6.axvspan(c[0],c[1],0.1,0.15,color="black",alpha=0.8)
            

            
    for t in telo.itertuples(index=False, name=None):
            ax6.axvspan(t[0],t[1],0.1,0.15,color="grey",alpha=0.8)
            
    for ge in genee.itertuples(index=False, name=None):
        if ge[3] == 'protein coding':
            if ge[7] == 'reverse':
                ax6.axvspan(ge[4],ge[5],0.2,0.3,color="rosybrown",alpha=0.3)
            elif ge[7] == 'forward':
                ax6.axvspan(ge[4],ge[5],0.8,0.9,color="rosybrown",alpha=0.3)
                
    for cc in ccon.itertuples(index=False, name=None):
            if cc[4] == 'reverse':
                ax6.axvspan(cc[2],cc[3],0.3,0.5,color="royalblue",alpha=0.5)
            elif cc[4] == 'forward':
                ax6.axvspan(cc[2],cc[3],0.5,0.8,color="royalblue",alpha=0.5)

    return ff

chromosome1 = Chromosome_plot(wtchr1, t2tchr1, t2rchr1, t2echr1, gene1, centro1, feat1, telo1)
chromosome1 = Chromosome_plot(wtchr1, sen1chr1, dbl8chr1, dschr1, gene1, centro1, feat1, telo1,con1)
chromosome2 = Chromosome_plot(wtchr2, sen1chr2, dbl8chr2, dschr2, gene2, centro2, feat2, telo2)
chromosome3 = Chromosome_plot(wtchr3, sen1chr3, dbl8chr3, dschr3, gene3, centro3, feat3, telo3)
chromosome1 = Chromosome_plot(wtchr3, t2tchr3, t2rchr3, t2echr3, gene3, centro3, feat3, telo3)

2+2
t2rchr1, t2rchr2, t2rchr3, t2r = Collect_data('FullMap.pombase220208-T2R.txt', 10000)
t2echr1, t2echr2, t2echr3, t2e = Collect_data('FullMap.pombase220208-T2E.txt', 10000)


def Chromosome_plot (cc1, cc2, cc3, cc4, centro, genee, telo):
    ff, (ax2, ax3, ax4, ax5, ax6) = plt.subplots(5,1, sharex=True)
    #ax1.set_title('WT')
    ax2.plot(cc1['Pos'], cc1['Wnorm'], color ='firebrick', alpha=0.8)
    ax2.plot(cc1['Pos'], cc1['Cnorm'], color ='steelblue', alpha=0.8)
   # ax2.set_ylim(-0.2,0.2)
    ax2.set_ylabel('WT (HpM)')
    #ax2.set_ylim(-5,5)
    ax2.set_ylabel('WT')

    ax3.plot(cc2['Pos'], cc2['Wnorm'], color ='firebrick', alpha=0.8)
    ax3.plot(cc2['Pos'], cc2['Cnorm'], color ='steelblue', alpha=0.8)
  #  ax3.set_ylim(-0.5,0.5)
   # ax3.set_ylabel('SEN1d (HpM)')
    ax3.set_ylabel('t2t (HpM)')

    ax6.set_xlabel('Chromosome position')
    ax4.plot(cc3['Pos'], cc3['Wnorm'], color ='firebrick', alpha=0.8)
    ax4.plot(cc3['Pos'], cc3['Cnorm'], color ='steelblue', alpha=0.8)
  #  ax4.set_ylim(-0.5,0.5)
    ax4.set_ylabel('t2r (HpM)')

    ax5.plot(cc4['Pos'], cc4['Wnorm'], color ='firebrick', alpha=0.8)
    ax5.plot(cc4['Pos'], cc4['Cnorm'], color ='steelblue', alpha=0.8)
  #  ax5.set_ylim(-0.5,0.5)
    ax5.set_ylabel('ETHANOL (HpM)')
 #   ax4.set_xlabel('Chromosome position')

 
    for c in centro.itertuples(index=False, name=None):
            ax6.axvspan(c[0],c[1],0.1,0.15,color="black",alpha=0.8)
            
    for t in telo.itertuples(index=False, name=None):
            ax6.axvspan(t[0],t[1],0.1,0.15,color="grey",alpha=0.8)
            
    for ge in genee.itertuples(index=False, name=None):
        if ge[3] == 'protein coding':
            if ge[7] == 'reverse':
                ax6.axvspan(ge[4],ge[5],0.2,0.3,color="rosybrown",alpha=0.3)
            elif ge[7] == 'forward':
                ax6.axvspan(ge[4],ge[5],0.8,0.9,color="rosybrown",alpha=0.3)

    return ff

chromosome1 = Chromosome_plot(wtchr1, t2tchr1, t2rchr1, t2echr1, centro1, feat1, telo1)






def Chromosome_plot (cc3, cc4, featurex, centro, genee, telo):
    ff, (ax4, ax5) = plt.subplots(2,1, sharex=True, sharey=True)
    #ax1.set_title('WT')

    ax5.set_xlabel('Chromosome position')
    ax4.plot(cc3['Pos'], cc3['Wnorm'], color ='firebrick', alpha=0.8)
    ax4.plot(cc3['Pos'], cc3['Cnorm'], color ='steelblue', alpha=0.8)
    ax4.set_ylim(-0.5,0.5)
    ax4.set_ylabel('t2r (HpM)')

    ax5.plot(cc4['Pos'], cc4['Wnorm'], color ='firebrick', alpha=0.8)
    ax5.plot(cc4['Pos'], cc4['Cnorm'], color ='steelblue', alpha=0.8)
    ax5.set_ylim(-0.5,0.5)
    ax5.set_ylabel('T2E (HpM)')
 #   ax4.set_xlabel('Chromosome position')


                
    for c in centro.itertuples(index=False, name=None):
            ax4.axvspan(c[0],c[1],0.1,0.15,color="black",alpha=0.8)
            ax5.axvspan(c[0],c[1],0.1,0.15,color="black",alpha=0.8)
            
    for t in telo.itertuples(index=False, name=None):
            ax4.axvspan(t[0],t[1],0.1,0.15,color="grey",alpha=0.8)
            ax5.axvspan(t[0],t[1],0.1,0.15,color="grey",alpha=0.8)


    return ff

chromosome1 = Chromosome_plot(t2rchr2, t2echr2, gene2, centro2, feat2, telo2)

    
#%%%


def filter_df1_by_ggenes(df1, ggenes):
    # Create an empty DataFrame to store the filtered rows
    filtered_df1 = pd.DataFrame()

    # Iterate through each row in ggenes
    for index, row in ggenes.iterrows():
        start = row['start']
        stop = row['end']
        chro = row['chro']

        # Filter rows in df1 based on the conditions 'Pos' > start, 'Pos' < stop, and 'Chr' == 'chro'
        filtered_rows = df1[(df1['Pos'] > start) & (df1['Pos'] < stop) & (df1['Chr'] == chro)]

        # Append the filtered rows to the result DataFrame
        filtered_df1 = pd.concat([filtered_df1, filtered_rows])

    # Reset the index of the resulting DataFrame (if needed)
    filtered_df1.reset_index(drop=True, inplace=True)

    return filtered_df1

sen1filt = filter_df1_by_ggenes(sen1, gggenes)
wtfilt = filter_df1_by_ggenes(wt, gggenes)
dbl8filt = filter_df1_by_ggenes(dbl8, gggenes)
dsfilt = filter_df1_by_ggenes(ds, gggenes)

controlsen1filt = filter_df1_by_ggenes(sen1, ccontrol)
controlwtfilt = filter_df1_by_ggenes(wt, ccontrol)
controldbl8filt = filter_df1_by_ggenes(dbl8, ccontrol)
conttroldsfilt = filter_df1_by_ggenes(ds, ccontrol)

sen1filt.to_csv("sen1filt.csv", index=False)
wtfilt.to_csv("wtfilt.csv", index=False)
dbl8filt.to_csv("dbl8filt.csv", index=False)
dsfilt.to_csv("dsfilt.csv", index=False)

controlsen1filt.to_csv("controlsen1filt.csv", index=False)
controlwtfilt.to_csv("controlwtfilt.csv", index=False)
controldbl8filt.to_csv("controldbl8filt.csv", index=False)
conttroldsfilt.to_csv("conttroldsfilt.csv", index=False)


#%%

output_path = "//Users/patricfernandez/Documents/911filter/wtstallfilt"
output_path = "/Users/patricfernandez/Documents/911filter/sen1stallfilter"


def filter_df1_by_ggenes_and_export(df1, ggenes, output_path):
    # Create an empty dictionary to store the filtered DataFrames
    filtered_dfs = {}

    # Iterate through each row in ggenes
    for index, row in ggenes.iterrows():
        start = row['start']
        stop = row['end']
        chro = row['chro']

        # Filter rows in df1 based on the conditions 'Pos' > start, 'Pos' < stop, and 'Chr' == 'chro'
        filtered_rows = df1[(df1['Pos'] > start) & (df1['Pos'] < stop) & (df1['Chr'] == chro)]

        # Store the filtered rows in the dictionary with a key representing the gene
        gene_name = row['ID']
        filtered_dfs[gene_name] = filtered_rows

        # Export the filtered DataFrame to a CSV file
        output_file = f"{output_path}/{gene_name}.tsv"
        filtered_rows.to_csv(output_file,sep='\t', index=False)

    return filtered_dfs


sen1filt = filter_df1_by_ggenes_and_export(sen1, gggenes, output_path)
dbl8filt = filter_df1_by_ggenes_and_export(dbl8, gggenes, output_path)
dsfitl = filter_df1_by_ggenes_and_export(ds, gggenes, output_path)
wtfilt = filter_df1_by_ggenes_and_export(wt, gggenes, output_path)



output_path = "/Users/patricfernandez/Documents/wt control filt"
output_path = '/Users/patricfernandez/Documents/911filter/dbl8controlfilter'

def filter_df1_by_ggenes_and_export(df1, ggenes, output_path):
    # Create an empty dictionary to store the filtered DataFrames
    filtered_dfs = {}

    # Iterate through each row in ggenes
    for index, row in ggenes.iterrows():
        start = row['start']
        stop = row['end']
        chro = row['chro']

        # Filter rows in df1 based on the conditions 'Pos' > start, 'Pos' < stop, and 'Chr' == 'chro'
        filtered_rows = df1[(df1['Pos'] > start) & (df1['Pos'] < stop) & (df1['Chr'] == chro)]

        # Store the filtered rows in the dictionary with a key representing the gene
        gene_name = row['ID']
        filtered_dfs[gene_name] = filtered_rows

        # Export the filtered DataFrame to a CSV file
        output_file = f"{output_path}/{gene_name}.csv"
        filtered_rows.to_csv(output_file, sep='\t', index=False)

    return filtered_dfs


sen1filt = filter_df1_by_ggenes_and_export(sen1, ccontrol, output_path)
dbl8filt = filter_df1_by_ggenes_and_export(dbl8, ccontrol, output_path)
dsfitl = filter_df1_by_ggenes_and_export(ds, ccontrol, output_path)
wtfilt = filter_df1_by_ggenes_and_export(wt, ccontrol, output_path)


sen1filt = filter_df1_by_ggenes_and_export(sen1, ccontrol, output_path)
dbl8filt = filter_df1_by_ggenes_and_export(dbl8, ccontrol, output_path)
dsfitl = filter_df1_by_ggenes_and_export(ds, ccontrol, output_path)
wtfilt = filter_df1_by_ggenes_and_export(wt, ccontrol, output_path)



#%%

def Findfeat(file):
    genes = pd.read_csv(file, delimiter="\t")
    print(genes)
    genes.loc[genes['Chromosome'] == "I", 'chro'] = 1
    genes.loc[genes['Chromosome'] == "II", 'chro'] = 2
    genes.loc[genes['Chromosome'] == "III", 'chro'] = 3
    genes.loc[genes['Strand'] == "forward", 'coding_strand'] = 'forward'
    genes.loc[genes['Strand'] == "reverse", 'coding_strand'] = 'reverse'
    
    genes = genes.rename(columns={'Start position': 'start'})
    genes = genes.rename(columns={'End position': 'end'})
    genes = genes.rename(columns={'Systematic ID': 'ID'})
    featfor = genes[(genes['coding_strand'] == 'forward')]
    featrev = genes[(genes['coding_strand'] == 'reverse')]
    
    
    
    return featfor, featrev, genes

featfor, featrev, ffeat = Findfeat('protein_coding_gene_list.tsv')


output_path = "/Users/patricfernandez/Documents/all_genes_filter_ds"

output_path = '/Users/patricfernandez/Documents/all_genes_filter_wt/'


def filter_df1_by_ggenes_and_export(df1, ggenes, output_path):
    # Create an empty dictionary to store the filtered DataFrames
    filtered_dfs = {}

    # Iterate through each row in ggenes
    for index, row in ggenes.iterrows():
        start = row['start']
        stop = row['end']
        chro = row['chro']

        # Filter rows in df1 based on the conditions 'Pos' > start, 'Pos' < stop, and 'Chr' == 'chro'
        filtered_rows = df1[(df1['Pos'] > start) & (df1['Pos'] < stop) & (df1['Chr'] == chro)]

        # Store the filtered rows in the dictionary with a key representing the gene
        gene_name = row['ID']
        filtered_dfs[gene_name] = filtered_rows

        # Export the filtered DataFrame to a CSV file
        output_file = f"{output_path}/{gene_name}.tsv"
        filtered_rows.to_csv(output_file,sep='\t', index=False)

    return filtered_dfs


dsfiltall = filter_df1_by_ggenes_and_export(ds, ffeat, output_path)
wtfiltall = filter_df1_by_ggenes_and_export(wt, ffeat, output_path)

#first step is replace "." with "a"
#then add FullMap.
#then replace .tsv with .txt

#%%

select10 = ffeat.sample(10)
select200 = ffeat.sample(200)


def rename_column(df, old_name, new_name):
    df.rename(columns={old_name: new_name}, inplace=True)
    return df

gggenes = rename_column(gggenes, "ID", "Systematic ID")

def discard_duplicates(ggenes, select200):
    merged = pd.merge(ggenes, select200, on="Systematic ID", how="inner")
    filtered = select200[~select200["Systematic ID"].isin(merged["Systematic ID"])]
    return filtered

select200 = discard_duplicates(gggenes, select200)
select10 = discard_duplicates(gggenes, select10)

select200.to_csv('select200.csv', index=False)
select10.to_csv('select10.csv', index=False)


select200 = rename_column(select200, "Strand", "coding_strand")
select10 = rename_column(select10, "Strand", "coding_strand")

select200 = rename_column(select200, "Start position", "start")
select10 = rename_column(select10, "Start position", "start")

select10.loc[select10['chro'] == "chr1", 'chro'] = 1
select10.loc[select10['chro'] == "chr2", 'chro'] = 2
select10.loc[select10['chro'] == "chr3", 'chro'] = 3

select200.loc[select200['chro'] == "chr1", 'chro'] = 1
select200.loc[select200['chro'] == "chr2", 'chro'] = 2
select200.loc[select200['chro'] == "chr3", 'chro'] = 3



#write a subset for let's say 1kb either side of a gene
def flank(file,frame1):
    data_output = pd.DataFrame()
    data_output_cnorm = pd.DataFrame()
    for i in range(len(file)):
        tempstart = file.iloc[i]["start"]
        tempchro = file.iloc[i]["chro"]
        tempstrand = file.iloc[i]['coding_strand']
        
        startminus = (tempstart) - 1000
        tempstarty = (tempstart) + 100

        
        if tempstrand == 'forward':
            tempsubsetfor = frame1.loc[(frame1['Pos']>= startminus) & (frame1['Pos'] <= tempstarty) & (frame1['Chr'] == tempchro)]
           # print(tempsubsetfor)
            tempsubsetfor = tempsubsetfor.reset_index(drop = True)
            line = pd.DataFrame({"Pos": tempstarty}, index=[1100])
            line2 = pd.DataFrame({"Pos": startminus}, index=[0])
            tempsubsetfor = tempsubsetfor.append(line, ignore_index=False)
            tempsubsetfor = tempsubsetfor.append(line2, ignore_index=False)
            
            tempsubsetfor = tempsubsetfor.sort_values('Pos', ascending=True)
        #?   # tempsubsetfor = tempsubsetfor.reset_index(drop = True)
          #  print(tempsubsetfor)
            #pos = tempsubsetfor.iloc[i]['Pos']
            #tpos = tempsubsetfor.iloc[0]['Pos']
            for x in range(len(tempsubsetfor)):
                pos = tempsubsetfor.iloc[x]['Pos']
               # print(pos)
                tpos = tempsubsetfor.iloc[0]['Pos']
               # print(tpos)
                
                #pos = tempsubsetfor.iloc[x]['Pos']
                tempsubsetfor.loc[tempsubsetfor.index[x], 'aPos'] = (pos - tpos)
                #print(tempsubsetfor) 
                #tempsubsetfor = tempsubsetfor.sort_values('aPos', ascending=False)
                #tempsubsetfor = tempsubsetfor.reset_index(drop = True)
           # print(tempsubsetfor)    
            tempsubsetfor = tempsubsetfor.set_index('aPos')
            
            

            data_blank = pd.DataFrame({
    
                "aPos" : np.arange(0.0, 1100, 1)}).merge(
        
                    tempsubsetfor, on = "aPos", how = "left"

                    )

            data_blank = data_blank[["aPos", "Wnorm", "Cnorm"]]
            
            data_output = pd.concat([data_output, data_blank['Wnorm']], axis =1, ignore_index=True)
            data_output_cnorm = pd.concat([data_output_cnorm, data_blank['Cnorm']], axis =1, ignore_index=True)
            
            data_output = data_output.fillna(0)

            data_output_cnorm = data_output_cnorm.fillna(0)
 #   data_output = data_output.T
  #  data_output_cnorm = data_output_cnorm.T
            
   # data_output['Wmean'] = data_output.mean(axis=1)
    #data_output_cnorm['Cmean'] = data_output_cnorm.mean(axis=1)
            
           #print(data_output)

    return data_output,data_output_cnorm

#hello = flank(gene1, csen1chr1)

wtforflanksenW, wtforflanksenC = flank(sen1stall, wt)
senforflanksenW, senforflanksenC = flank(sen1stall, sen1)
dblforflanksenW, dblforflanksenC = flank(sen1stall, dbl8)
dsforflanksenW, dsforflanksenC = flank(sen1stall, ds)

wtforflankdblW, wtforflankdblC = flank(dbl8stall, wt)
senforflankdblW, senforflankdblC = flank(dbl8stall, sen1)
dblforflankdblW, dblforflankdblC = flank(dbl8stall, dbl8)
dsforflankdblW, dsforflankdblC = flank(dbl8stall, ds)

wtforflankconW, wtforflankconC = flank(ccontrol, wt)
senforflankconW, senforflankconC = flank(ccontrol, sen1)
dblforflankconW, dblforflankconC = flank(ccontrol, dbl8)
dsforflankconW, dsforflankconC = flank(ccontrol, ds)

wtforflankdsW, wtforflankdsC = flank(doublestall, wt)
senforflankdsW, senforflankdsC = flank(doublestall, sen1)
dblforflankdsW, dblforflankdsC = flank(doublestall, dbl8)
dsforflankdsW, dsforflankdsC = flank(doublestall, ds)





#VERSION 2 INCLUDE REVERSE STRAND GENES
def flank(file,frame1):
    data_output = pd.DataFrame()
    data_output_cnorm = pd.DataFrame()
    
    for i in range(len(file)):
        tempstrand = file.iloc[i]['coding_strand']
        tempchro = file.iloc[i]["chro"]
        if tempstrand == 'forward':
            tempstart = file.iloc[i]["start"]
            tempchro = file.iloc[i]["chro"]

            startminus = (tempstart) - 1000
            tempstarty = (tempstart) + 100

        
            tempsubsetfor = frame1.loc[(frame1['Pos']>= startminus) & (frame1['Pos'] <= tempstarty) & (frame1['Chr'] == tempchro)]

            tempsubsetfor = tempsubsetfor.reset_index(drop = True)
            line = pd.DataFrame({"Pos": tempstarty}, index=[1100])
            line2 = pd.DataFrame({"Pos": startminus}, index=[0])
            tempsubsetfor = tempsubsetfor.append(line, ignore_index=False)
            tempsubsetfor = tempsubsetfor.append(line2, ignore_index=False)
            
            tempsubsetfor = tempsubsetfor.sort_values('Pos', ascending=True)

            for x in range(len(tempsubsetfor)):
                pos = tempsubsetfor.iloc[x]['Pos']

                tpos = tempsubsetfor.iloc[0]['Pos']

                tempsubsetfor.loc[tempsubsetfor.index[x], 'aPos'] = (pos - tpos)

            tempsubsetfor = tempsubsetfor.set_index('aPos')
                        

            data_blank = pd.DataFrame({
    
                "aPos" : np.arange(0.0, 1100, 1)}).merge(
        
                    tempsubsetfor, on = "aPos", how = "left"

                    )

            data_blank = data_blank[["aPos", "Wnorm", "Cnorm"]]
            
            data_output = pd.concat([data_output, data_blank['Wnorm']], axis =1, ignore_index=True)
            data_output_cnorm = pd.concat([data_output_cnorm, data_blank['Cnorm']], axis =1, ignore_index=True)
            
            data_output = data_output.fillna(0)

            data_output_cnorm = data_output_cnorm.fillna(0)
        
        
        if tempstrand == 'reverse':
            tempstartr = file.iloc[i]["start"]            

            startminusr = (tempstartr) + 1000
            tempstartyr = (tempstartr) - 100

        
            tempsubsetrev = frame1.loc[(frame1['Pos']>= tempstartyr) & (frame1['Pos'] <= startminusr) & (frame1['Chr'] == tempchro)]

            tempsubsetrev = tempsubsetrev.reset_index(drop = True)
            line = pd.DataFrame({"Pos": tempstartyr}, index=[1100])
            line2 = pd.DataFrame({"Pos": startminusr}, index=[0])
            tempsubsetrev = tempsubsetrev.append(line, ignore_index=False)
            tempsubsetrev = tempsubsetrev.append(line2, ignore_index=False)
            
            tempsubsetrev = tempsubsetrev.sort_values('Pos', ascending=True)

            for x in range(len(tempsubsetrev)):
                pos = tempsubsetrev.iloc[x]['Pos']

                tpos = tempsubsetrev.iloc[0]['Pos']

                tempsubsetrev.loc[tempsubsetrev.index[x], 'aPos'] = (pos - tpos)

            tempsubsetrev = tempsubsetrev.set_index('aPos')
                        

            data_blankr = pd.DataFrame({
    
                "aPos" : np.arange(0.0, 1100, 1)}).merge(
        
                    tempsubsetrev, on = "aPos", how = "left"

                    )

            data_blankr = data_blankr[["aPos", "Wnorm", "Cnorm"]]
            
            data_output = pd.concat([data_output, data_blankr['Wnorm']], axis =1, ignore_index=True)
            data_output_cnorm = pd.concat([data_output_cnorm, data_blankr['Cnorm']], axis =1, ignore_index=True)
            
            data_output = data_output.fillna(0)

            data_output_cnorm = data_output_cnorm.fillna(0)

    return data_output,data_output_cnorm

wtforflanksenW, wtforflanksenC = flank(sen1stall, wt)
senforflanksenW, senforflanksenC = flank(sen1stall, sen1)
dblforflanksenW, dblforflanksenC = flank(sen1stall, dbl8)
dsforflanksenW, dsforflanksenC = flank(sen1stall, ds)

wtforflankdblW, wtforflankdblC = flank(dbl8stall, wt)
senforflankdblW, senforflankdblC = flank(dbl8stall, sen1)
dblforflankdblW, dblforflankdblC = flank(dbl8stall, dbl8)
dsforflankdblW, dsforflankdblC = flank(dbl8stall, ds)

wtforflankconW, wtforflankconC = flank(ccontrol, wt)
senforflankconW, senforflankconC = flank(ccontrol, sen1)
dblforflankconW, dblforflankconC = flank(ccontrol, dbl8)
dsforflankconW, dsforflankconC = flank(ccontrol, ds)

wtforflankdsW, wtforflankdsC = flank(doublestall, wt)
senforflankdsW, senforflankdsC = flank(doublestall, sen1)
dblforflankdsW, dblforflankdsC = flank(doublestall, dbl8)
dsforflankdsW, dsforflankdsC = flank(doublestall, ds)


wtforflanksel10W, wtforflanksel10C = flank(select10, wt)
senforflanksel10W, senforflanksel10C = flank(select10, sen1)
dblforflanksel10W, dblforflanksel10C = flank(select10, dbl8)
dsforflanksel10W, dsforflanksel10C = flank(select10, ds)


wtforflanksel200W, wtforflanksel200C = flank(select200, wt)
senforflanksel200W, senforflanksel200C = flank(select200, sen1)
dblforflanksel200W, dblforflanksel200C = flank(select200, dbl8)
dsforflanksel200W, dsforflanksel200C = flank(select200, ds)






#VERSION 3 INCLUDE REVERSE STRAND GENES and ttake normscore
def flank(file,frame1):
    data_output = pd.DataFrame()
   # data_output_cnorm = pd.DataFrame()
    
    for i in range(len(file)):
        tempstrand = file.iloc[i]['coding_strand']
        tempchro = file.iloc[i]["chro"]
        if tempstrand == 'forward':
            tempstart = file.iloc[i]["start"]
            tempchro = file.iloc[i]["chro"]

            startminus = (tempstart) - 3000
            tempstarty = (tempstart) + 3000

        
            tempsubsetfor = frame1.loc[(frame1['Pos']>= startminus) & (frame1['Pos'] <= tempstarty) & (frame1['Chr'] == tempchro)]

            tempsubsetfor = tempsubsetfor.reset_index(drop = True)
            line = pd.DataFrame({"Pos": tempstarty}, index=[6000])
            line2 = pd.DataFrame({"Pos": startminus}, index=[0])
            tempsubsetfor = tempsubsetfor.append(line, ignore_index=False)
            tempsubsetfor = tempsubsetfor.append(line2, ignore_index=False)
            
            tempsubsetfor = tempsubsetfor.sort_values('Pos', ascending=True)

            for x in range(len(tempsubsetfor)):
                pos = tempsubsetfor.iloc[x]['Pos']

                tpos = tempsubsetfor.iloc[0]['Pos']

                tempsubsetfor.loc[tempsubsetfor.index[x], 'aPos'] = (pos - tpos)

            tempsubsetfor = tempsubsetfor.set_index('aPos')
                        

            data_blank = pd.DataFrame({
    
                "aPos" : np.arange(0.0, 6000, 1)}).merge(
        
                    tempsubsetfor, on = "aPos", how = "left"

                    )

            data_blank = data_blank[["aPos", "normscore"]]
            
            data_output = pd.concat([data_output, data_blank['normscore']], axis =1, ignore_index=True)
         #   data_output_cnorm = pd.concat([data_output_cnorm, data_blank['Cnorm']], axis =1, ignore_index=True)
            
            data_output = data_output.fillna(0)

         #   data_output_cnorm = data_output_cnorm.fillna(0)
        
        
        if tempstrand == 'reverse':
            tempstartr = file.iloc[i]["end"]            

            startminusr = (tempstartr) + 3000
            tempstartyr = (tempstartr) - 3000

        
            tempsubsetrev = frame1.loc[(frame1['Pos']>= tempstartyr) & (frame1['Pos'] <= startminusr) & (frame1['Chr'] == tempchro)]

            tempsubsetrev = tempsubsetrev.reset_index(drop = True)
            line = pd.DataFrame({"Pos": tempstartyr}, index=[6000])
            line2 = pd.DataFrame({"Pos": startminusr}, index=[0])
            tempsubsetrev = tempsubsetrev.append(line, ignore_index=False)
            tempsubsetrev = tempsubsetrev.append(line2, ignore_index=False)
            
            tempsubsetrev = tempsubsetrev.sort_values('Pos', ascending=True)

            for x in range(len(tempsubsetrev)):
                pos = tempsubsetrev.iloc[0]['Pos']

                tpos = tempsubsetrev.iloc[x]['Pos']

                tempsubsetrev.loc[tempsubsetrev.index[x], 'aPos'] = (pos - tpos)

            tempsubsetrev = tempsubsetrev.set_index('aPos')
                        

            data_blankr = pd.DataFrame({
    
                "aPos" : np.arange(0.0, 6000, 1)}).merge(
        
                    tempsubsetrev, on = "aPos", how = "left"

                    )

            data_blankr = data_blankr[["aPos", "normscore"]]
            
            data_output = pd.concat([data_output, data_blankr['normscore']], axis =1, ignore_index=True)
           # data_output_cnorm = pd.concat([data_output_cnorm, data_blankr['Cnorm']], axis =1, ignore_index=True)
            
            data_output = data_output.fillna(0)

           # data_output_cnorm = data_output_cnorm.fillna(0)

    return data_output



def flank(file,frame1):
    data_output = pd.DataFrame()
   # data_output_cnorm = pd.DataFrame()
    
    for i in range(len(file)):
        tempstrand = file.iloc[i]['coding_strand']
        tempchro = file.iloc[i]["chro"]
        if tempstrand == 'forward':
            tempstart = file.iloc[i]["start"]
            tempchro = file.iloc[i]["chro"]

            startminus = (tempstart) - 3000
            tempstarty = (tempstart) + 3000

        
            tempsubsetfor = frame1.loc[(frame1['Pos']>= startminus) & (frame1['Pos'] <= tempstarty) & (frame1['Chr'] == tempchro)]

            tempsubsetfor = tempsubsetfor.reset_index(drop = True)
            line = pd.DataFrame({"Pos": tempstarty}, index=[6000])
            line2 = pd.DataFrame({"Pos": startminus}, index=[0])
            tempsubsetfor = pd.concat([tempsubsetfor, line, line2], axis=0, ignore_index=False)
            
            tempsubsetfor = tempsubsetfor.sort_values('Pos', ascending=True)

            for x in range(len(tempsubsetfor)):
                pos = tempsubsetfor.iloc[x]['Pos']

                tpos = tempsubsetfor.iloc[0]['Pos']

                tempsubsetfor.loc[tempsubsetfor.index[x], 'aPos'] = (pos - tpos)

            tempsubsetfor = tempsubsetfor.set_index('aPos')
                        

            data_blank = pd.DataFrame({
    
                "aPos" : np.arange(0.0, 6000, 1)}).merge(
        
                    tempsubsetfor, on = "aPos", how = "left"

                    )

            data_blank = data_blank[["aPos", "normscore"]]
            
            data_output = pd.concat([data_output, data_blank['normscore']], axis =1, ignore_index=True)
         #   data_output_cnorm = pd.concat([data_output_cnorm, data_blank['Cnorm']], axis =1, ignore_index=True)
            
            data_output = data_output.fillna(0)

         #   data_output_cnorm = data_output_cnorm.fillna(0)
        
        
        if tempstrand == 'reverse':
            tempstartr = file.iloc[i]["end"]            

            startminusr = (tempstartr) + 3000
            tempstartyr = (tempstartr) - 3000

        
            tempsubsetrev = frame1.loc[(frame1['Pos']>= tempstartyr) & (frame1['Pos'] <= startminusr) & (frame1['Chr'] == tempchro)]

            tempsubsetrev = tempsubsetrev.reset_index(drop = True)
            line = pd.DataFrame({"Pos": tempstartyr}, index=[6000])
            line2 = pd.DataFrame({"Pos": startminusr}, index=[0])
            tempsubsetrev = pd.concat([tempsubsetrev, line, line2], axis=0, ignore_index=False)
            
            tempsubsetrev = tempsubsetrev.sort_values('Pos', ascending=True)

            for x in range(len(tempsubsetrev)):
                pos = tempsubsetrev.iloc[0]['Pos']

                tpos = tempsubsetrev.iloc[x]['Pos']

                tempsubsetrev.loc[tempsubsetrev.index[x], 'aPos'] = (pos - tpos)

            tempsubsetrev = tempsubsetrev.set_index('aPos')
                        

            data_blankr = pd.DataFrame({
    
                "aPos" : np.arange(0.0, 6000, 1)}).merge(
        
                    tempsubsetrev, on = "aPos", how = "left"

                    )

            data_blankr = data_blankr[["aPos", "normscore"]]
            
            data_output = pd.concat([data_output, data_blankr['normscore']], axis =1, ignore_index=True)
           # data_output_cnorm = pd.concat([data_output_cnorm, data_blankr['Cnorm']], axis =1, ignore_index=True)
            
            data_output = data_output.fillna(0)

           # data_output_cnorm = data_output_cnorm.fillna(0)

    return data_output

wtforflanksenW = flank(sen1stall, wt)
senforflanksenW = flank(sen1stall, sen1)
dblforflanksenW = flank(sen1stall, dbl8)
dsforflanksenW = flank(sen1stall, ds)

wtforflankdblW = flank(dbl8stall, wt)
senforflankdblW = flank(dbl8stall, sen1)
dblforflankdblW = flank(dbl8stall, dbl8)
dsforflankdblW = flank(dbl8stall, ds)

wtforflankdsW = flank(doublestall, wt)
senforflankdsW = flank(doublestall, sen1)
dblforflankdsW = flank(doublestall, dbl8)
dsforflankdsW = flank(doublestall, ds)

wtforflankconW = flank(ccontrol, wt)
senforflankconW = flank(ccontrol, sen1)
dblforflankconW = flank(ccontrol, dbl8)
dsforflankconW = flank(ccontrol, ds)




wtforflanksel10W = flank(select10, wt)
senforflanksel10W = flank(select10, sen1)
dblforflanksel10W = flank(select10, dbl8)
dsforflanksel10W = flank(select10, ds)


wtforflanksel200W = flank(select200, wt)
senforflanksel200W = flank(select200, sen1)
dblforflanksel200W = flank(select200, dbl8)
dsforflanksel200W = flank(select200, ds)

#%%%
#using pol2 plot 

def Find(file):
    genes = pd.read_csv(file, delimiter=",")
    print(genes)
    genesfor = genes[(genes['coding_strand'] == 'forward')]
    genesrev = genes[(genes['coding_strand'] == 'reverse')]


    genesfor['Sbinpos'] = genesfor['start']/100
    genesfor['Sbinpos'] = genesfor['Sbinpos'].astype(int)
    genesfor['Sbinpos'] = genesfor['Sbinpos']*100
    genesfor['Ebinpos'] = genesfor['end']/100
    genesfor['Ebinpos'] = genesfor['Ebinpos'].astype(int)
    genesfor['Ebinpos'] = genesfor['Ebinpos']*100


    genesrev['Sbinposr'] = genesrev['end']/100
    genesrev['Sbinposr'] = genesrev['Sbinposr'].astype(int)
    genesrev['Sbinposr'] = genesrev['Sbinposr']*100
    genesrev['Ebinposr'] = genesrev['start']/100
    genesrev['Ebinposr'] = genesrev['Ebinposr'].astype(int)
    genesrev['Ebinposr'] = genesrev['Ebinposr']*100

    return genesfor, genesrev, genes

genesfor, genesrev, ggenes = Find("dbl8_stall_sites_direction.txt")
controlfor, controlrev, ccontrol = Find('new_control.csv')


sen1stall = ggenes[(ggenes['genotype'] == 'sen1D')]
dbl8stall = ggenes[(ggenes['genotype'] == 'dbl8D')]
doublestall = ggenes[(ggenes['genotype'] == 'sen1dbl8DD_unique')]



def RStart(genesfor, chr1, chr2, chr3, p, k):
    xx =[]
    
    for g in genesfor.itertuples():
        if g[7] == k:
            if g[8] == p:
                if g[5] == 'chr1':
                    x = chr1.index[chr1['Pos'] == g[9]].tolist()
                    xx.append(x) 

                if g[5] == 'chr2':
                    x2 = chr2.index[chr2['Pos'] == g[9]].tolist()
                    xx.append(x2)

                if g[5] == 'chr3':
                    x3 = chr3.index[chr3['Pos'] == g[9]].tolist()
                    xx.append(x3)
                
    return xx
dsforxxR = RStart(genesfor, sen1chr1, sen1chr2, sen1chr3, 'sen1dbl8DD_unique', 'rightward')
dsrevxxR = RStart(genesrev, sen1chr1, sen1chr2, sen1chr3, 'sen1dbl8DD_unique', 'rightward')
dsforxxL = RStart(genesfor, sen1chr1, sen1chr2, sen1chr3, 'sen1dbl8DD_unique', 'leftward')
dsrevxxL = RStart(genesrev, sen1chr1, sen1chr2, sen1chr3, 'sen1dbl8DD_unique', 'leftward')


def REnd(genesfor, chr1, chr2, chr3, p, k):
    xxe =[]
    for ge in genesfor.itertuples():
        if ge[7] == k:
            if ge[8] == p:
                if ge[5] == 'chr1':
                    xe = chr1.index[chr1['Pos'] == ge[10]].tolist()
                    xxe.append(xe) 

                if ge[5] == 'chr2':
                    xe2 = chr2.index[chr2['Pos'] == ge[10]].tolist()
                    xxe.append(xe2)

                if ge[5] == 'chr3':
                    xe3 = chr3.index[chr3['Pos'] == ge[10]].tolist()
                    xxe.append(xe3)
    return xxe

dsforxxeR = REnd(genesfor, sen1chr1, sen1chr2, sen1chr3, 'sen1dbl8DD_unique', 'rightward')
dsrevxxeR = REnd(genesrev, sen1chr1, sen1chr2, sen1chr3, 'sen1dbl8DD_unique', 'rightward')
dsforxxeL = REnd(genesfor, sen1chr1, sen1chr2, sen1chr3, 'sen1dbl8DD_unique', 'leftward')
dsrevxxeL = REnd(genesrev, sen1chr1, sen1chr2, sen1chr3, 'sen1dbl8DD_unique', 'leftward')

#lets make two list, one is ds hh start 
# rev right and 
# for left 




#for wt , forward coding strand genes only 
#TSS
def Start(genesfor, chr1, chr2, chr3, p):
    xx =[]
    
    for g in genesfor.itertuples():
        if g[8] == p:
            if g[5] == 'chr1':
                x = chr1.index[chr1['Pos'] == g[9]].tolist()
                xx.append(x) 

            if g[5] == 'chr2':
                x2 = chr2.index[chr2['Pos'] == g[9]].tolist()
                xx.append(x2)

            if g[5] == 'chr3':
                x3 = chr3.index[chr3['Pos'] == g[9]].tolist()
                xx.append(x3)
                
    return xx

senforxx = Start(genesfor, sen1chr1, sen1chr2, sen1chr3, 'sen1D')
senrevxx = Start(genesrev, sen1chr1, sen1chr2, sen1chr3, 'sen1D')
dbl8forxx = Start(genesfor, dbl8chr1, dbl8chr2, dbl8chr3, 'dbl8D')
dbl8revxx = Start(genesrev, dbl8chr1, dbl8chr2, dbl8chr3, 'dbl8D')
#dsforxx = Start(genesfor, dschr1, dschr2, dschr3, 'sen1dbl8DD_unique')
#dsrevxx = Start(genesrev, dschr1, dschr2, dschr3, 'sen1dbl8DD_unique')


def ControlStart(genesfor, chr1, chr2, chr3):
    xx =[]
    
    for g in genesfor.itertuples():
            if g[2] == 'chr1':
                x = chr1.index[chr1['Pos'] == g[8]].tolist()
                xx.append(x) 

            if g[2] == 'chr2':
                x2 = chr2.index[chr2['Pos'] == g[8]].tolist()
                xx.append(x2)

            if g[2] == 'chr3':
                x3 = chr3.index[chr3['Pos'] == g[8]].tolist()
                xx.append(x3)
                
    return xx
controlforxx = ControlStart(controlfor, wtchr1, wtchr2, wtchr3)
controlrevxx = ControlStart(controlrev, wtchr1, wtchr2, wtchr3)


def ControlEnd(genesfor, chr1, chr2, chr3):
    xx =[]
    
    for g in genesfor.itertuples():
            if g[2] == 'chr1':
                x = chr1.index[chr1['Pos'] == g[9]].tolist()
                xx.append(x) 

            if g[2] == 'chr2':
                x2 = chr2.index[chr2['Pos'] == g[9]].tolist()
                xx.append(x2)

            if g[2] == 'chr3':
                x3 = chr3.index[chr3['Pos'] == g[9]].tolist()
                xx.append(x3)
                
    return xx
controlforxxe = ControlEnd(controlfor, wtchr1, wtchr2, wtchr3)
controlrevxxe = ControlEnd(controlrev, wtchr1, wtchr2, wtchr3)

#TES
def End(genesfor, chr1, chr2, chr3, p):
    xxe =[]
    for ge in genesfor.itertuples():
        if ge[8] == p:
            if ge[5] == 'chr1':
                xe = chr1.index[chr1['Pos'] == ge[10]].tolist()
                xxe.append(xe) 

            if ge[5] == 'chr2':
                xe2 = chr2.index[chr2['Pos'] == ge[10]].tolist()
                xxe.append(xe2)

            if ge[5] == 'chr3':
                xe3 = chr3.index[chr3['Pos'] == ge[10]].tolist()
                xxe.append(xe3)
    return xxe

senforxxe = End(genesfor, sen1chr1, sen1chr2, sen1chr3, 'sen1D')
senrevxxe = End(genesrev, sen1chr1, sen1chr2, sen1chr3, 'sen1D')
dbl8forxxe = End(genesfor, dbl8chr1, dbl8chr2, dbl8chr3, 'dbl8D')
dbl8revxxe = End(genesrev, dbl8chr1, dbl8chr2, dbl8chr3, 'dbl8D')
#dsforxxe = End(genesfor, dschr1, dschr2, dschr3, 'sen1dbl8DD_unique')
#dsrevxxe = End(genesrev, dschr1, dschr2, dschr3, 'sen1dbl8DD_unique')


#selection of random 200 genes 

#####PLEASE NOTE RIGHT NOW ITS JUST ALL GENES
def selectionlStart(genesfor, chr1, chr2, chr3):
    xx =[]
    
    for g in genesfor.itertuples():
            if g[9] == 'chr1':
                x = chr1.index[chr1['pos'] == g[10]].tolist()
                xx.append(x) 

            if g[9] == 'chr2':
                x2 = chr2.index[chr2['pos'] == g[10]].tolist()
                xx.append(x2)

            if g[9] == 'chr3':
                x3 = chr3.index[chr3['pos'] == g[10]].tolist()
                xx.append(x3)
                
    return xx
selforxx = selectionlStart(featfor, wtchr1, wtchr2, wtchr3)
selrevxx = selectionlStart(featrev, wtchr1, wtchr2, wtchr3)


def selectionEnd(genesfor, chr1, chr2, chr3):
    xx =[]
    
    for g in genesfor.itertuples():
            if g[9] == 'chr1':
                x = chr1.index[chr1['pos'] == g[11]].tolist()
                xx.append(x) 

            if g[9] == 'chr2':
                x2 = chr2.index[chr2['pos'] == g[11]].tolist()
                xx.append(x2)

            if g[9] == 'chr3':
                x3 = chr3.index[chr3['pos'] == g[11]].tolist()
                xx.append(x3)
                
    return xx
selforxxe = selectionEnd(featfor, wtchr1, wtchr2, wtchr3)
selrevxxe = selectionEnd(featrev, wtchr1, wtchr2, wtchr3)


def Gene_bits(list1,list2,frame):

    yucky=[]
    for z, ze in zip(list1, list2):
        yy = frame.loc[z[0]:ze[0],'normscore'].tolist()
        yucky.append(yy)
    stalls = pd.DataFrame(yucky)   

    return stalls

sen1stallswt = Gene_bits(senforxx, senforxxe, wt)
sen1stallssen1 = Gene_bits(senforxx, senforxxe, sen1)
sen1stallsdbl8 = Gene_bits(senforxx, senforxxe, dbl8)
sen1stallsds = Gene_bits(senforxx, senforxxe, ds)

dbl8stallswt = Gene_bits(dbl8forxx, dbl8forxxe, wt)
dbl8stallssen1 = Gene_bits(dbl8forxx, dbl8forxxe, sen1)
dbl8stallsdbl8 = Gene_bits(dbl8forxx, dbl8forxxe, dbl8)
dbl8stallsds = Gene_bits(dbl8forxx, dbl8forxxe, ds)

#for left double
dsstallswtL = Gene_bits(dsforxxL, dsforxxeL, wt)
dsstallssen1L = Gene_bits(dsforxxL, dsforxxeL, sen1)
dsstallsdbl8L = Gene_bits(dsforxxL, dsforxxeL, dbl8)
dsstallsdsL = Gene_bits(dsforxxL, dsforxxeL, ds)

#for right double 
dsstallswtR = Gene_bits(dsforxxR, dsforxxeR, wt)
dsstallssen1R = Gene_bits(dsforxxR, dsforxxeR, sen1)
dsstallsdbl8R = Gene_bits(dsforxxR, dsforxxeR, dbl8)
dsstallsdsR = Gene_bits(dsforxxR, dsforxxeR, ds)

#for control 
constallswt = Gene_bits(controlforxx, controlforxxe, wt)
constallssen1 = Gene_bits(controlforxx, controlforxxe, sen1)
constallsdbl8 = Gene_bits(controlforxx, controlforxxe, dbl8)
constallsds = Gene_bits(controlforxx, controlforxxe, ds)

#for random selection
selbodywt = Gene_bits(selforxx, selforxxe, wt)
selbodysen1 = Gene_bits(selforxx, selforxxe, sen1)
selbodydbl8 = Gene_bits(selforxx, selforxxe, dbl8)
selbodyds = Gene_bits(selforxx, selforxxe, ds)



def Rev_gene_bits (list1, list2, frame):
    alldatar =[]
    for y, ye  in zip(list1, list2):
        yyr = frame.loc[ye[0]:y[0],'normscore'].tolist()
        alldatar.append(yyr)
    stalls = pd.DataFrame(alldatar)
    return stalls

sen1stallsrwt = Rev_gene_bits(senrevxx, senrevxxe, wt)
sen1stallsrsen1 = Rev_gene_bits(senrevxx, senrevxxe, sen1)
sen1stallsrdbl8 = Rev_gene_bits(senrevxx, senrevxxe, dbl8)
sen1stallsrds = Rev_gene_bits(senrevxx, senrevxxe, ds)

dbl8stallsrwt = Rev_gene_bits(dbl8revxx, dbl8revxxe, wt)
dbl8stallsrsen1 = Rev_gene_bits(dbl8revxx, dbl8revxxe, sen1)
dbl8stallsrdbl8 = Rev_gene_bits(dbl8revxx, dbl8revxxe, dbl8)
dbl8stallsrds = Rev_gene_bits(dbl8revxx, dbl8revxxe, ds)

#rev right double
dsstallsrwtR = Rev_gene_bits(dsrevxxR, dsrevxxeR, wt)
dsstallsrsen1R = Rev_gene_bits(dsrevxxR, dsrevxxeR, sen1)
dsstallsrdbl8R = Rev_gene_bits(dsrevxxR, dsrevxxeR, dbl8)
dsstallsrdsR = Rev_gene_bits(dsrevxxR, dsrevxxeR, ds)

#rev left double 
dsstallsrwtL = Rev_gene_bits(dsrevxxL, dsrevxxeL, wt)
dsstallsrsen1L = Rev_gene_bits(dsrevxxL, dsrevxxeL, sen1)
dsstallsrdbl8L = Rev_gene_bits(dsrevxxL, dsrevxxeL, dbl8)
dsstallsrdsL = Rev_gene_bits(dsrevxxL, dsrevxxeL, ds)


#rev control 
constallsrwt = Rev_gene_bits(controlrevxx, controlrevxxe, wt)
constallsrsen1 = Rev_gene_bits(controlrevxx, controlrevxxe, sen1)
constallsrdbl8 = Rev_gene_bits(controlrevxx, controlrevxxe, dbl8)
constallsrds = Rev_gene_bits(controlrevxx, controlrevxxe, ds)

#rev random


selbodyrwt = Rev_gene_bits(selrevxx, selrevxxe, wt)
selbodyrsen1 = Rev_gene_bits(selrevxx, selrevxxe, sen1)
selbodyrdbl8 = Rev_gene_bits(selrevxx, selrevxxe, dbl8)
selbodyrds = Rev_gene_bits(selrevxx, selrevxxe, ds)

#%%
#seperate W and C norm now 

def Gene_bits(list1,list2,frame):

    yucky=[]
    for z, ze in zip(list1, list2):
        yy = frame.loc[z[0]:ze[0],'Wnorm'].tolist()
        yucky.append(yy)
    stalls = pd.DataFrame(yucky)   

    return stalls

wsen1stallswt = Gene_bits(senforxx, senforxxe, wt)
wsen1stallssen1 = Gene_bits(senforxx, senforxxe, sen1)
wsen1stallsdbl8 = Gene_bits(senforxx, senforxxe, dbl8)
wsen1stallsds = Gene_bits(senforxx, senforxxe, ds)

wdbl8stallswt = Gene_bits(dbl8forxx, dbl8forxxe, wt)
wdbl8stallssen1 = Gene_bits(dbl8forxx, dbl8forxxe, sen1)
wdbl8stallsdbl8 = Gene_bits(dbl8forxx, dbl8forxxe, dbl8)
wdbl8stallsds = Gene_bits(dbl8forxx, dbl8forxxe, ds)

#for left double
wdsstallswtL = Gene_bits(dsforxxL, dsforxxeL, wt)
wdsstallssen1L = Gene_bits(dsforxxL, dsforxxeL, sen1)
wdsstallsdbl8L = Gene_bits(dsforxxL, dsforxxeL, dbl8)
wdsstallsdsL = Gene_bits(dsforxxL, dsforxxeL, ds)

#for right double 
wdsstallswtR = Gene_bits(dsforxxR, dsforxxeR, wt)
wdsstallssen1R = Gene_bits(dsforxxR, dsforxxeR, sen1)
wdsstallsdbl8R = Gene_bits(dsforxxR, dsforxxeR, dbl8)
wdsstallsdsR = Gene_bits(dsforxxR, dsforxxeR, ds)

#for control 
wconstallswt = Gene_bits(controlforxx, controlforxxe, wt)
wconstallssen1 = Gene_bits(controlforxx, controlforxxe, sen1)
wconstallsdbl8 = Gene_bits(controlforxx, controlforxxe, dbl8)
wconstallsds = Gene_bits(controlforxx, controlforxxe, ds)




def Rev_gene_bits (list1, list2, frame):
    alldatar =[]
    for y, ye  in zip(list1, list2):
        yyr = frame.loc[ye[0]:y[0],'Wnorm'].tolist()
        alldatar.append(yyr)
    stalls = pd.DataFrame(alldatar)
    return stalls

wsen1stallsrwt = Rev_gene_bits(senrevxx, senrevxxe, wt)
wsen1stallsrsen1 = Rev_gene_bits(senrevxx, senrevxxe, sen1)
wsen1stallsrdbl8 = Rev_gene_bits(senrevxx, senrevxxe, dbl8)
wsen1stallsrds = Rev_gene_bits(senrevxx, senrevxxe, ds)

wdbl8stallsrwt = Rev_gene_bits(dbl8revxx, dbl8revxxe, wt)
wdbl8stallsrsen1 = Rev_gene_bits(dbl8revxx, dbl8revxxe, sen1)
wdbl8stallsrdbl8 = Rev_gene_bits(dbl8revxx, dbl8revxxe, dbl8)
wdbl8stallsrds = Rev_gene_bits(dbl8revxx, dbl8revxxe, ds)

#rev right double
wdsstallsrwtR = Rev_gene_bits(dsrevxxR, dsrevxxeR, wt)
wdsstallsrsen1R = Rev_gene_bits(dsrevxxR, dsrevxxeR, sen1)
wdsstallsrdbl8R = Rev_gene_bits(dsrevxxR, dsrevxxeR, dbl8)
wdsstallsrdsR = Rev_gene_bits(dsrevxxR, dsrevxxeR, ds)

#rev left double 
wdsstallsrwtL = Rev_gene_bits(dsrevxxL, dsrevxxeL, wt)
wdsstallsrsen1L = Rev_gene_bits(dsrevxxL, dsrevxxeL, sen1)
wdsstallsrdbl8L = Rev_gene_bits(dsrevxxL, dsrevxxeL, dbl8)
wdsstallsrdsL = Rev_gene_bits(dsrevxxL, dsrevxxeL, ds)


#rev control 
wconstallsrwt = Rev_gene_bits(controlrevxx, controlrevxxe, wt)
wconstallsrsen1 = Rev_gene_bits(controlrevxx, controlrevxxe, sen1)
wconstallsrdbl8 = Rev_gene_bits(controlrevxx, controlrevxxe, dbl8)
wconstallsrds = Rev_gene_bits(controlrevxx, controlrevxxe, ds)





def Gene_bits(list1,list2,frame):

    yucky=[]
    for z, ze in zip(list1, list2):
        yy = frame.loc[z[0]:ze[0],'Cnorm'].tolist()
        yucky.append(yy)
    stalls = pd.DataFrame(yucky)   

    return stalls

sen1stallswt = Gene_bits(senforxx, senforxxe, wt)
sen1stallssen1 = Gene_bits(senforxx, senforxxe, sen1)
sen1stallsdbl8 = Gene_bits(senforxx, senforxxe, dbl8)
sen1stallsds = Gene_bits(senforxx, senforxxe, ds)

dbl8stallswt = Gene_bits(dbl8forxx, dbl8forxxe, wt)
dbl8stallssen1 = Gene_bits(dbl8forxx, dbl8forxxe, sen1)
dbl8stallsdbl8 = Gene_bits(dbl8forxx, dbl8forxxe, dbl8)
dbl8stallsds = Gene_bits(dbl8forxx, dbl8forxxe, ds)

#for left double
dsstallswtL = Gene_bits(dsforxxL, dsforxxeL, wt)
dsstallssen1L = Gene_bits(dsforxxL, dsforxxeL, sen1)
dsstallsdbl8L = Gene_bits(dsforxxL, dsforxxeL, dbl8)
dsstallsdsL = Gene_bits(dsforxxL, dsforxxeL, ds)

#for right double 
dsstallswtR = Gene_bits(dsforxxR, dsforxxeR, wt)
dsstallssen1R = Gene_bits(dsforxxR, dsforxxeR, sen1)
dsstallsdbl8R = Gene_bits(dsforxxR, dsforxxeR, dbl8)
dsstallsdsR = Gene_bits(dsforxxR, dsforxxeR, ds)

#for control 
constallswt = Gene_bits(controlforxx, controlforxxe, wt)
constallssen1 = Gene_bits(controlforxx, controlforxxe, sen1)
constallsdbl8 = Gene_bits(controlforxx, controlforxxe, dbl8)
constallsds = Gene_bits(controlforxx, controlforxxe, ds)



def Rev_gene_bits (list1, list2, frame):
    alldatar =[]
    for y, ye  in zip(list1, list2):
        yyr = frame.loc[ye[0]:y[0],'Cnorm'].tolist()
        alldatar.append(yyr)
    stalls = pd.DataFrame(alldatar)
    return stalls

sen1stallsrwt = Rev_gene_bits(senrevxx, senrevxxe, wt)
sen1stallsrsen1 = Rev_gene_bits(senrevxx, senrevxxe, sen1)
sen1stallsrdbl8 = Rev_gene_bits(senrevxx, senrevxxe, dbl8)
sen1stallsrds = Rev_gene_bits(senrevxx, senrevxxe, ds)

dbl8stallsrwt = Rev_gene_bits(dbl8revxx, dbl8revxxe, wt)
dbl8stallsrsen1 = Rev_gene_bits(dbl8revxx, dbl8revxxe, sen1)
dbl8stallsrdbl8 = Rev_gene_bits(dbl8revxx, dbl8revxxe, dbl8)
dbl8stallsrds = Rev_gene_bits(dbl8revxx, dbl8revxxe, ds)

#rev right double
dsstallsrwtR = Rev_gene_bits(dsrevxxR, dsrevxxeR, wt)
dsstallsrsen1R = Rev_gene_bits(dsrevxxR, dsrevxxeR, sen1)
dsstallsrdbl8R = Rev_gene_bits(dsrevxxR, dsrevxxeR, dbl8)
dsstallsrdsR = Rev_gene_bits(dsrevxxR, dsrevxxeR, ds)

#rev left double 
dsstallsrwtL = Rev_gene_bits(dsrevxxL, dsrevxxeL, wt)
dsstallsrsen1L = Rev_gene_bits(dsrevxxL, dsrevxxeL, sen1)
dsstallsrdbl8L = Rev_gene_bits(dsrevxxL, dsrevxxeL, dbl8)
dsstallsrdsL = Rev_gene_bits(dsrevxxL, dsrevxxeL, ds)


#rev control 
constallsrwt = Rev_gene_bits(controlrevxx, controlrevxxe, wt)
constallsrsen1 = Rev_gene_bits(controlrevxx, controlrevxxe, sen1)
constallsrdbl8 = Rev_gene_bits(controlrevxx, controlrevxxe, dbl8)
constallsrds = Rev_gene_bits(controlrevxx, controlrevxxe, ds)






def Expand_genes(stall_df):
    alldata = []
    for row in stall_df.iterrows():
        xx = row[1]
        xx = xx.dropna()
        length = len(xx)
        num = 1000/length
        print (num*length)

        te = np.arange(0,1000,num, dtype = int)

        expand = np.zeros(1000)

        for itter, val in enumerate(te):
            if itter<=length-2:
                expand[val:te[itter+1]] = xx[itter]
            else: continue
        expand[te[length-1]:] = xx[length-1]
        alldata.append(expand)
    return alldata

#control
cw = Expand_genes(constallswt)
cs = Expand_genes(constallssen1)
cb = Expand_genes(constallsdbl8)
cds = Expand_genes(constallsds)

cwr = Expand_genes(constallsrwt)
csr = Expand_genes(constallsrsen1)
cbr = Expand_genes(constallsrdbl8)
cdsr = Expand_genes(constallsrds)



#sen1
aw = Expand_genes(sen1stallswt)
ass = Expand_genes(sen1stallssen1)
ad = Expand_genes(sen1stallsdbl8)
ads = Expand_genes(sen1stallsds)

raw= Expand_genes(sen1stallsrwt)
ras = Expand_genes(sen1stallsrsen1)
rad = Expand_genes(sen1stallsrdbl8)
rads = Expand_genes(sen1stallsrds)



#dbl8
dbw = Expand_genes(dbl8stallswt)
dbss = Expand_genes(dbl8stallssen1)
dbd = Expand_genes(dbl8stallsdbl8)
dbds = Expand_genes(dbl8stallsds)

rdbw= Expand_genes(dbl8stallsrwt)
rdbs = Expand_genes(dbl8stallsrsen1)
rdbd = Expand_genes(dbl8stallsrdbl8)
rdbds = Expand_genes(dbl8stallsrds)


#ds right 
dswR = Expand_genes(dsstallswtR)
dssR = Expand_genes(dsstallssen1R)
dsbR = Expand_genes(dsstallsdbl8R)
dsdsR = Expand_genes(dsstallsdsR)

rdswR = Expand_genes(dsstallsrwtR)
rdssR = Expand_genes(dsstallsrsen1R)
rdsbR = Expand_genes(dsstallsrdbl8R)
rdsdsR = Expand_genes(dsstallsrdsR)

#ds left
dswL = Expand_genes(dsstallswtL)
dssL = Expand_genes(dsstallssen1L)
dsbL = Expand_genes(dsstallsdbl8L)
dsdsL = Expand_genes(dsstallsdsL)

rdswL = Expand_genes(dsstallsrwtL)
rdssL = Expand_genes(dsstallsrsen1L)
rdsbL = Expand_genes(dsstallsrdbl8L)
rdsdsL = Expand_genes(dsstallsrdsL)



def Expand_genes(stall_df):
    alldata = []
    for row in stall_df.iterrows():
        xx = row[1]
        xx = xx.dropna()
        length = len(xx)
        num = 1000/length
        print (num*length)

        te = np.arange(0,1000,num, dtype = int)

        expand = np.zeros(1000)

        for itter, val in enumerate(te):
            if itter<=length-2:
                expand[val:te[itter+1]] = xx[itter]
            else: continue
        expand[te[length-1]:] = xx[length-1]
        alldata.append(expand)
    return alldata

#control
wcw = Expand_genes(wconstallswt)
wcs = Expand_genes(wconstallssen1)
wcb = Expand_genes(wconstallsdbl8)
wcds = Expand_genes(wconstallsds)

wcwr = Expand_genes(wconstallsrwt)
wcsr = Expand_genes(wconstallsrsen1)
wcbr = Expand_genes(wconstallsrdbl8)
wcdsr = Expand_genes(wconstallsrds)



#sen1
waw = Expand_genes(wsen1stallswt)
wass = Expand_genes(wsen1stallssen1)
wad = Expand_genes(wsen1stallsdbl8)
wads = Expand_genes(wsen1stallsds)

wraw= Expand_genes(wsen1stallsrwt)
wras = Expand_genes(wsen1stallsrsen1)
wrad = Expand_genes(wsen1stallsrdbl8)
wrads = Expand_genes(wsen1stallsrds)



#dbl8
wdbw = Expand_genes(wdbl8stallswt)
wdbss = Expand_genes(wdbl8stallssen1)
wdbd = Expand_genes(wdbl8stallsdbl8)
wdbds = Expand_genes(wdbl8stallsds)

wrdbw= Expand_genes(wdbl8stallsrwt)
wrdbs = Expand_genes(wdbl8stallsrsen1)
wrdbd = Expand_genes(wdbl8stallsrdbl8)
wrdbds = Expand_genes(wdbl8stallsrds)


#ds right 
wdswR = Expand_genes(wdsstallswtR)
wdssR = Expand_genes(wdsstallssen1R)
wdsbR = Expand_genes(wdsstallsdbl8R)
wdsdsR = Expand_genes(wdsstallsdsR)

wrdswR = Expand_genes(wdsstallsrwtR)
wrdssR = Expand_genes(wdsstallsrsen1R)
wrdsbR = Expand_genes(wdsstallsrdbl8R)
wrdsdsR = Expand_genes(wdsstallsrdsR)

#ds left
wdswL = Expand_genes(wdsstallswtL)
wdssL = Expand_genes(wdsstallssen1L)
wdsbL = Expand_genes(wdsstallsdbl8L)
wdsdsL = Expand_genes(wdsstallsdsL)

wrdswL = Expand_genes(wdsstallsrwtL)
wrdssL = Expand_genes(wdsstallsrsen1L)
wrdsbL = Expand_genes(wdsstallsrdbl8L)
wrdsdsL = Expand_genes(wdsstallsrdsL)


def concati(rev, forward):
    forwards = np.stack(forward, axis=0 )
    revs = np.stack(rev, axis=0 )
    rev1 = revs[:, ::-1]
    
    new = np.concatenate([forwards,rev1])
    return new

#reverse LEFT, forward right 
HTw = concati(rdswL, dswR)
HTs = concati(rdssL, dssR)
HTb = concati(rdsbL, dsbR)
HTds = concati (rdsdsL,dsdsR)

#reverse right + forward leeft 
HHw = concati(rdswR, dswL)
HHs = concati(rdssR, dssL)
HHb = concati(rdsbR, dsbL)
HHds = concati(rdsdsR, dsdsL)

sa = concati(raw, aw)
sas = concati(ras,ass)
sab = concati(rad,ad)
sads = concati(rads,ads)

dbl8w = concati(rdbw, dbw)
dbl8s = concati(rdbs, dbss)
dbl8b = concati(rdbd, dbd)
dbl8ds = concati(rdbds, dbds)

caw = concati(cwr, cw)
cas = concati(csr, cs)
cab = concati(cbr, cb)
cads = concati(cdsr, cds)


def concati(rev, forward):
    forwards = np.stack(forward, axis=0 )
    revs = np.stack(rev, axis=0 )
    rev1 = revs[:, ::-1]
    
    new = np.concatenate([forwards,rev1])
    return new

#reverse LEFT, forward right 
wHTw = concati(wrdswL, wdswR)
wHTs = concati(wrdssL, wdssR)
wHTb = concati(wrdsbL, wdsbR)
wHTds = concati (wrdsdsL, wdsdsR)

#reverse right + forward leeft 
wHHw = concati(wrdswR, wdswL)
wHHs = concati(wrdssR, wdssL)
wHHb = concati(wrdsbR, wdsbL)
wHHds = concati(wrdsdsR, wdsdsL)

wsa = concati(wraw, waw)
wsas = concati(wras,wass)
wsab = concati(wrad,wad)
wsads = concati(wrads,wads)

wdbl8w = concati(wrdbw, wdbw)
wdbl8s = concati(wrdbs,wdbss)
wdbl8b = concati(wrdbd, wdbd)
wdbl8ds = concati(wrdbds, wdbds)

wcaw = concati(wcwr, wcw)
wcas = concati(wcsr, wcs)
wcab = concati(wcbr, wcb)
wcads = concati(wcdsr, wcds)


import seaborn as sns
fxdsl,((ax1, ax3, ax5,ax7)) = plt.subplots(4,1)  
fxdsl.tight_layout()
cbar_ax = fxdsl.add_axes([.91, .3, .03, .4])
sns.heatmap(HHw, cmap = 'coolwarm', ax=ax1, cbar=True, cbar_ax=cbar_ax)
ax1.axes.yaxis.set_ticks([])
ax1.axes.xaxis.set_ticks([])
ax1.set_title('S∆D∆ Head-Head collisions')
ax1.set_ylabel('WT')
sns.heatmap(HHs, cmap = 'coolwarm', ax=ax3, cbar=True, cbar_ax=cbar_ax)
ax3.axes.yaxis.set_ticks([])
ax3.set_ylabel('Sen1∆')
ax3.axes.xaxis.set_ticks([])
sns.heatmap(HHb, cmap = 'coolwarm', ax=ax5,cbar=True, cbar_ax=cbar_ax)
ax5.axes.yaxis.set_ticks([])
ax5.set_ylabel('Dbl8∆')
ax5.axes.xaxis.set_ticks([])
sns.heatmap(HHds, cmap = 'coolwarm', ax=ax7, cbar=True, cbar_ax=cbar_ax)
ax7.axes.yaxis.set_ticks([])
ax7.set_ylabel('Sen1∆Dbl8∆')
ax7.set_xticks([0,1000])
ax7.set_xticklabels(['TSS','TES'])

fuckkkkkk, (ax2, ax4, ax6, ax8) = plt.subplots(4,1)
fuckkkkkk.tight_layout()
cbar_ax = fuckkkkkk.add_axes([.91, .3, .03, .4])
sns.heatmap(HTw, cmap = 'coolwarm', ax=ax2,  cbar=True, cbar_ax=cbar_ax)
ax2.set_ylabel('WT')
ax2.axes.yaxis.set_ticks([])
ax2.set_title('S∆D∆ co-directional collisions')
ax2.axes.xaxis.set_ticks([])
sns.heatmap(HTs, cmap = 'coolwarm', ax=ax4, cbar=True, cbar_ax=cbar_ax)
ax4.axes.yaxis.set_ticks([])
ax4.set_ylabel('Sen1∆')
ax4.axes.xaxis.set_ticks([])
sns.heatmap(HTb, cmap = 'coolwarm', ax=ax6, cbar=True, cbar_ax=cbar_ax)
ax6.axes.yaxis.set_ticks([])
ax6.set_ylabel('Dbl8∆')  
ax6.axes.xaxis.set_ticks([])
sns.heatmap(HTds, cmap = 'coolwarm', ax=ax8,cbar=True, cbar_ax=cbar_ax)
ax8.axes.yaxis.set_ticks([])
ax8.set_ylabel('Sen1∆Dbl8∆')

ax8.set_xticks([0,1000])
ax8.set_xticklabels(['TSS','TES'])


#CONTROL HEAT MAP
fix, ((ax1, ax3, ax5,ax7)) = plt.subplots(4,1)  
fix.tight_layout()
cbar_ax = fix.add_axes([.91, .3, .03, .4])
sns.heatmap(wcaw, cmap = 'coolwarm_r', ax=ax1,cbar=True, cbar_ax=cbar_ax)
ax1.axes.yaxis.set_ticks([])
ax1.axes.xaxis.set_ticks([])
ax1.set_title('Forward Strand')
ax1.set_ylabel('WT')
sns.heatmap(wcas, cmap = 'coolwarm_r', ax=ax3,cbar=True, cbar_ax=cbar_ax)

ax3.axes.yaxis.set_ticks([])
ax3.set_ylabel('Sen1∆')
ax3.axes.xaxis.set_ticks([])
sns.heatmap(wcab, cmap = 'coolwarm_r', ax=ax5, cbar=True, cbar_ax= cbar_ax)
ax5.axes.yaxis.set_ticks([])
ax5.set_ylabel('Dbl8∆')
ax5.axes.xaxis.set_ticks([])
    
sns.heatmap(wcads, cmap = 'coolwarm_r', ax=ax7,cbar=True, cbar_ax=cbar_ax)
ax7.axes.yaxis.set_ticks([])
ax7.set_ylabel('Sen1∆Dbl8∆')

ax7.set_xticks([0,1000])
ax7.set_xticklabels(['TSS','TES'])


#DBL8 heatmap !!!!!
fxb, ((ax1, ax3, ax5,ax7)) = plt.subplots(4,1)  
fxb.tight_layout()
cbar_ax = fxb.add_axes([.91, .3, .03, .4])
sns.heatmap(wdbl8w, cmap = 'coolwarm_r', ax=ax1,  cbar=True, cbar_ax=cbar_ax)
ax1.axes.yaxis.set_ticks([])
ax1.axes.xaxis.set_ticks([])
ax1.set_title('Forward Strand')
ax1.set_ylabel('WT')
sns.heatmap(wdbl8s, cmap = 'coolwarm_r', ax=ax3,  cbar=True, cbar_ax=cbar_ax)
ax3.axes.yaxis.set_ticks([])
ax3.set_ylabel('Sen1∆')
ax3.axes.xaxis.set_ticks([])
sns.heatmap(wdbl8b, cmap = 'coolwarm_r', ax=ax5,  cbar=True, cbar_ax=cbar_ax)
ax5.axes.yaxis.set_ticks([])
ax5.set_ylabel('Dbl8∆')
ax5.axes.xaxis.set_ticks([])
sns.heatmap(wdbl8ds, cmap = 'coolwarm_r', ax=ax7,  cbar=True, cbar_ax=cbar_ax)
ax7.axes.yaxis.set_ticks([])
ax7.set_ylabel('Sen1∆Dbl8∆')

ax7.set_xticks([0,1000])
ax7.set_xticklabels(['TSS','TES'])

#SEN1 heatmap bitches 
fx, ((ax1, ax3, ax5,ax7)) = plt.subplots(4,1)  
fx.tight_layout()
cbar_ax = fx.add_axes([.91, .3, .03, .4])
sns.heatmap(sa, cmap = 'coolwarm', ax=ax1, cbar=True, cbar_ax=cbar_ax)
ax1.axes.yaxis.set_ticks([])
ax1.axes.xaxis.set_ticks([])
ax1.set_title('Forward Strand')
ax1.set_ylabel('WT')
sns.heatmap(sas, cmap = 'coolwarm', ax=ax3, cbar=True, cbar_ax=cbar_ax)
ax3.axes.yaxis.set_ticks([])
ax3.set_ylabel('Sen1∆')
ax3.axes.xaxis.set_ticks([])
sns.heatmap(sab, cmap = 'coolwarm', ax=ax5,  cbar=True, cbar_ax=cbar_ax)
ax5.axes.yaxis.set_ticks([])
ax5.set_ylabel('Dbl8∆')
ax5.axes.xaxis.set_ticks([])
sns.heatmap(sads, cmap = 'coolwarm', ax=ax7, cbar=True, cbar_ax=cbar_ax)
ax7.axes.yaxis.set_ticks([])
ax7.set_ylabel('Sen1∆Dbl8∆')


ax7.set_xticks([0,1000])
ax7.set_xticklabels(['TSS','TES'])
  
  
def pileup (thing):
    #thing = -thing
    #print(thing)
    #print(-thing)
    whack = np.stack( thing, axis=0 )
    stuff = whack.mean(axis=0)

    wow = pd.DataFrame(stuff, columns=['OriginalData'])

    # Apply smoothing using a rolling window
    wow['SmoothedData'] = wow['OriginalData'].rolling(window=50, min_periods=1).mean()
    wowe = wow[['SmoothedData']]

    return wowe
#sen compiled


#hh ht 
HHwline = pileup(HHw)
HHsline = pileup(HHs)
HHbline = pileup(HHb)
HHdsline = pileup(HHds)
HTwline = pileup(HTw)
HTsline = pileup(HTs)
HTbline = pileup(HTb)
HTdsline = pileup(HTds)


#control gene body
cwline = pileup(caw)
csline = pileup(cas)
cbline = pileup(cab)
cdsline = pileup(cads)

rcwline = pileup(cwr)
rcsline = pileup(csr)
rcbline = pileup(cbr)
rcdsline = pileup(cdsr)


#sen
awline = pileup(sa)
assline = pileup(sas)
adline = pileup(sab)
adsline= pileup(sads)


#dbl8
dbwline = pileup(dbl8w)
dbssline = pileup(dbl8s)
dbdline = pileup(dbl8b)
dbdsline= pileup(dbl8ds)

#WATSON

#hh ht 
wHHwline = pileup(wHHw)
wHHsline = pileup(wHHs)
wHHbline = pileup(wHHb)
wHHdsline = pileup(wHHds)
wHTwline = pileup(wHTw)
wHTsline = pileup(wHTs)
wHTbline = pileup(wHTb)
wHTdsline = pileup(wHTds)


#control gene body
wcwline = pileup(wcaw)
wcsline = pileup(wcas)
wcbline = pileup(wcab)
wcdsline = pileup(wcads)

wrcwline = pileup(wcwr)
wrcsline = pileup(wcsr)
wrcbline = pileup(wcbr)
wrcdsline = pileup(wcdsr)


#sen
wawline = pileup(wsa)
wassline = pileup(wsas)
wadline = pileup(wsab)
wadsline= pileup(wsads)


#dbl8
wdbwline = pileup(wdbl8w)
wdbssline = pileup(wdbl8s)
wdbdline = pileup(wdbl8b)
wdbdsline= pileup(wdbl8ds)




figgy, (ax1) =plt.subplots(1, sharey=True)

x = np.arange(0, 1000, 1)
ax1.plot(cwline.index, cwline['SmoothedData'], color = 'black', alpha =0.5, linewidth=1)
ax1.plot(wcwline.index, wcwline['SmoothedData'], color = 'black', alpha =0.5, linewidth=1)
ax1.plot(csline.index, csline['SmoothedData'], color = 'deepskyblue', alpha=0.5, linewidth=1)
ax1.plot(wcsline.index, wcsline['SmoothedData'], color = 'deepskyblue', alpha=0.5, linewidth=1)
ax1.plot(cbline.index, cbline['SmoothedData'], color = 'orange', alpha=0.5, linewidth=1)
ax1.plot(wcbline.index, wcbline['SmoothedData'], color = 'orange', alpha=0.5, linewidth=1)
ax1.plot(cdsline.index, cdsline['SmoothedData'], color = 'tomato', alpha=0.5, linewidth=1)
ax1.plot(wcdsline.index, wcdsline['SmoothedData'], color = 'tomato', alpha=0.5, linewidth=1)

ax1.legend(loc='best')



figgy, (ax1) =plt.subplots(1, sharey=True)
ax1.plot(awline.index, awline['SmoothedData'], color = 'black', alpha =0.8, linewidth=1)
ax1.plot(wawline.index, wawline['SmoothedData'], color = 'black', alpha =0.8, linewidth=1)


ax1.plot(assline.index, assline['SmoothedData'], color = 'deepskyblue', alpha=0.8, linewidth=1)
ax1.plot(wassline.index, wassline['SmoothedData'], color = 'deepskyblue', alpha=0.8, linewidth=1)

ax1.plot(adline.index, adline['SmoothedData'], color = 'orange', alpha=0.8, linewidth=1)
ax1.plot(wadline.index, wadline['SmoothedData'], color = 'orange', alpha=0.8, linewidth=1)

ax1.plot(adsline.index, adsline['SmoothedData'], color = 'tomato', alpha=0.8, linewidth=1)
ax1.plot(wadsline.index, wadsline['SmoothedData'], color = 'tomato', alpha=0.8, linewidth=1)
x = np.arange(0, 1000, 1)

ax1.set_xticks([0,1000])
ax1.set_xticklabels(['TSS','TES'])
#ax1.set_ylim([-4, 124])
ax1.set_ylabel('normalised reads')
ax1.legend(loc='best')

ax1.set_xticks([0,1000])
ax1.set_xticklabels(['TSS','TES'])
ax1.legend(loc='best')
ax1.set_ylabel('normalised reads')
#ax1.fill_between(x, y1a, y2a, facecolor='mediumaquamarine', alpha =0.2)
#y2 = linerev
#ax2.fill_between(x, y2, facecolor='green', alpha =0.3)
#filkggy, (ax1) =plt.subplots(1,1, sharey=True)

figgy, (ax1) =plt.subplots(1,1, sharey=True)
ax1.plot(dbwline.index, dbwline['SmoothedData'], color = 'black', alpha =0.8, linewidth=1)
ax1.plot(wdbwline.index, wdbwline['SmoothedData'], color = 'black', alpha =0.8, linewidth=1)
ax1.plot(dbssline.index, dbssline['SmoothedData'], color = 'deepskyblue', alpha=0.8, linewidth=1)
ax1.plot(wdbssline.index, wdbssline['SmoothedData'], color = 'deepskyblue', alpha=0.8, linewidth=1)
ax1.plot(dbdline.index, dbdline['SmoothedData'], color = 'orange', alpha=0.8, linewidth=1)
ax1.plot(wdbdline.index, wdbdline['SmoothedData'], color = 'orange', alpha=0.8, linewidth=1)
ax1.plot(dbdsline.index, dbdsline['SmoothedData'], color = 'tomato', alpha=0.8, linewidth=1)
ax1.plot(wdbdsline.index, wdbdsline['SmoothedData'], color = 'tomato', alpha=0.8, linewidth=1)
y3 = np.array(bwCON['lower']).flatten()
y4 = np.array(bwCON['upper']).flatten()
ax1.fill_between(x, y3, y4, facecolor='grey', alpha =0.2)
y3s = np.array(bsCON['lower']).flatten()
y4s = np.array(bsCON['upper']).flatten()
ax1.fill_between(x, y3s, y4s, facecolor='deepskyblue', alpha =0.1,)
y3b = np.array(bbCON['lower']).flatten()
y4b = np.array(bbCON['upper']).flatten()
ax1.fill_between(x, y3b, y4b, facecolor='gold', alpha =0.2)
y3ds = np.array(bdsCON['lower']).flatten()
y4ds = np.array(bdsCON['upper']).flatten()
ax1.fill_between(x, y3ds, y4ds, facecolor='indianred', alpha =0.1)
ax1.legend(loc='best')
ax1.set_xticks([0,1000])
ax1.set_xticklabels(['TSS','TES'])
#ax1.set_ylim([-4, 124])
ax1.legend(loc='best')
ax1.set_ylabel('normalised reads')


figgy, (ax1) =plt.subplots(1, sharey=True)
#ds left
ax1.plot(HHwline.index, HHwline['SmoothedData'], color = 'black', alpha =0.8, linewidth=1)
ax1.plot(wHHwline.index, wHHwline['SmoothedData'], color = 'black', alpha =0.8, linewidth=1)
ax1.plot(HHsline.index, HHsline['SmoothedData'], color = 'deepskyblue', alpha=0.8, linewidth=1)
ax1.plot(wHHsline.index, wHHsline['SmoothedData'], color = 'deepskyblue', alpha=0.8, linewidth=1)
ax1.plot(HHbline.index, HHbline['SmoothedData'], color = 'orange', alpha=0.8, linewidth=1)
ax1.plot(wHHbline.index, wHHbline['SmoothedData'], color = 'orange', alpha=0.8, linewidth=1)
ax1.plot(HHdsline.index, HHdsline['SmoothedData'], color = 'tomato', alpha=0.8, linewidth=1)
ax1.plot(wHHdsline.index, wHHdsline['SmoothedData'], color = 'tomato', alpha=0.8, linewidth=1)
y5 = np.array(HHwCON['lower']).flatten()
y6= np.array(HHwCON['upper']).flatten()
ax1.fill_between(x, y5, y6, facecolor='grey', alpha =0.2)
y5s = np.array(HHsCON['lower']).flatten()
y6s= np.array(HHsCON['upper']).flatten()
#ax1.fill_between(x, y5s, y6s, facecolor='teal', alpha =0.2)
y5b = np.array(HHbCON['lower']).flatten()
y6b= np.array(HHbCON['upper']).flatten()
#ax1.fill_between(x, y5b, y6b, facecolor='teal', alpha =0.2)
y5ds = np.array(HHdsCON['lower']).flatten()
y6ds= np.array(HHdsCON['upper']).flatten()
ax1.fill_between(x, y5ds, y6ds, facecolor='indianred', alpha =0.2)
#ax1.legend(loc='best')


ax1.set_xticks([0,1000])
ax1.set_xticklabels(['TSS','TES'])
ax1.legend(loc='best')
ax1.set_ylabel('normalised reads')

figgy, (ax1) =plt.subplots(1, sharey=True)
#ds right
ax1.plot(HTwline.index, HTwline['SmoothedData'], color = 'black', alpha =0.8, linewidth=1)
ax1.plot(wHTwline.index, wHTwline['SmoothedData'], color = 'black', alpha =0.8, linewidth=1)
ax1.plot(HTsline.index, HTsline['SmoothedData'], color = 'deepskyblue', alpha=0.8, linewidth=1)
ax1.plot(wHTsline.index, wHTsline['SmoothedData'], color = 'deepskyblue', alpha=0.8, linewidth=1)
ax1.plot(HTbline.index, HTbline['SmoothedData'], color = 'orange', alpha=0.8, linewidth=1)
ax1.plot(wHTbline.index, wHTbline['SmoothedData'], color = 'orange', alpha=0.8, linewidth=1)
ax1.plot(HTdsline.index, HTdsline['SmoothedData'], color = 'tomato', alpha=0.8, linewidth=1)
ax1.plot(wHTdsline.index, wHTdsline['SmoothedData'], color = 'tomato', alpha=0.8, linewidth=1)
y7 = np.array(HTwCON['lower']).flatten()
y8= np.array(HTwCON['upper']).flatten()
ax1.fill_between(x, y7, y8, facecolor='grey', alpha =0.2)
y7s = np.array(HTsCON['lower']).flatten()
y8s= np.array(HTsCON['upper']).flatten()
#ax1.fill_between(x, y7s, y8s, facecolor='salmon', alpha =0.2)
y7b = np.array(HTbCON['lower']).flatten()
y8b= np.array(HTbCON['upper']).flatten()
#ax1.fill_between(x, y7b, y8b, facecolor='salmon', alpha =0.2)
y7ds = np.array(HTdsCON['lower']).flatten()
y8ds= np.array(HTdsCON['upper']).flatten()
ax1.fill_between(x, y7ds, y8ds, facecolor='salmon', alpha =0.2, )
#ax1.legend(loc='best')

ax1.set_xticks([0,1000])
ax1.set_xticklabels(['TSS','TES'])
ax1.legend(loc='best')
ax1.set_ylabel('normalised reads')
#ax1.set_ylim([-4, 124])


figgy, (ax1) =plt.subplots(1, sharey=True)
ax1.plot(cwline.index, cwline['SmoothedData'], color = 'black', alpha =0.5, linewidth=1)
ax1.plot(csline.index, csline['SmoothedData'], color = 'deepskyblue', alpha=0.5, linewidth=1)
ax1.plot(cbline.index, cbline['SmoothedData'], color = 'orange', alpha=0.5, linewidth=1)
ax1.plot(cdsline.index, cdsline['SmoothedData'], color = 'tomato', alpha=0.5, linewidth=1)
y11 = np.array(cawCON['lower']).flatten()
y10= np.array(cawCON['upper']).flatten()
ax1.fill_between(x, y10, y11, facecolor='grey', alpha =0.3,)
y11s = np.array(casCON['lower']).flatten()
y10s= np.array(casCON['upper']).flatten()
#ax1.fill_between(x, y10s, y11s, facecolor='mediumaquamarine', alpha =0.3)
y11b = np.array(cabCON['lower']).flatten()
y10b = np.array(cabCON['upper']).flatten()
#ax1.fill_between(x, y10b, y11b, facecolor='mediumaquamarine', alpha =0.3)
y11ds = np.array(cadsCON['lower']).flatten()
y10ds= np.array(cadsCON['upper']).flatten()
ax1.fill_between(x, y10ds, y11ds, facecolor='mediumaquamarine', alpha =0.3)
ax1.legend(loc='best')



ax1.set_xticks([0,1000])
ax1.set_xticklabels(['TSS','TES'])
ax1.legend(loc='best')
ax1.set_ylabel('normalised reads')
#ax1.set_title('forward strand')


figgy, (ax1) =plt.subplots(1, sharey=True)
ax1.plot(swline.index, swline[0], color = 'black', alpha =0.5, linewidth=1)
ax1.plot(ssline.index, ssline[0], color = 'deepskyblue', alpha=0.5, linewidth=1)
ax1.plot(sbline.index, sbline[0], color = 'orange', alpha=0.5, linewidth=1)
ax1.plot(sdsline.index, sdsline[0], color = 'tomato', alpha=0.5, linewidth=1)
y9 = np.array(sewCON['lower']).flatten()
y10= np.array(sewCON['upper']).flatten()
ax1.fill_between(x, y9, y10, facecolor='grey', alpha =0.2)



y9s = np.array(sesCON['lower']).flatten()
y10s= np.array(sesCON['upper']).flatten()
ax1.fill_between(x, y9s, y10s, facecolor='deepskyblue', alpha =0.3)
y9b = np.array(sebCON['lower']).flatten()
y10b= np.array(sebCON['upper']).flatten()
ax1.fill_between(x, y9b, y10b, facecolor='orange', alpha =0.3)
y9ds = np.array(sedsCON['lower']).flatten()
y10ds= np.array(sedsCON['upper']).flatten()
ax1.fill_between(x, y9ds, y10ds, facecolor='tomato', alpha =0.3)
ax1.legend(loc='best')

ax1.set_xticks([0,1000])
ax1.set_xticklabels(['TSS','TES'])
ax1.legend(loc='best')
ax1.set_ylabel('normalised reads')


ax1.set_ylim([-4, 124])

#%%

def Expand_genes(stall_df):
    alldata = []
    for row in stall_df.iterrows():
        xx = row[1]
        xx = xx.dropna()
        length = len(xx)
        num = 1000/length
        print (num*length)

        te = np.arange(0,1000,num, dtype = int)

        expand = np.zeros(1000)

        for itter, val in enumerate(te):
            if itter<=length-2:
                expand[val:te[itter+1]] = xx[itter]
            else: continue
        expand[te[length-1]:] = xx[length-1]
        alldata.append(expand)
    return alldata

#control
cw = Expand_genes(constallswt)
cs = Expand_genes(constallssen1)
cb = Expand_genes(constallsdbl8)
cds = Expand_genes(constallsds)

cwr = Expand_genes(constallsrwt)
csr = Expand_genes(constallsrsen1)
cbr = Expand_genes(constallsrdbl8)
cdsr = Expand_genes(constallsrds)



#sen1
aw = Expand_genes(sen1stallswt)
ass = Expand_genes(sen1stallssen1)
ad = Expand_genes(sen1stallsdbl8)
ads = Expand_genes(sen1stallsds)

raw= Expand_genes(sen1stallsrwt)
ras = Expand_genes(sen1stallsrsen1)
rad = Expand_genes(sen1stallsrdbl8)
rads = Expand_genes(sen1stallsrds)



#dbl8
dbw = Expand_genes(dbl8stallswt)
dbss = Expand_genes(dbl8stallssen1)
dbd = Expand_genes(dbl8stallsdbl8)
dbds = Expand_genes(dbl8stallsds)

rdbw= Expand_genes(dbl8stallsrwt)
rdbs = Expand_genes(dbl8stallsrsen1)
rdbd = Expand_genes(dbl8stallsrdbl8)
rdbds = Expand_genes(dbl8stallsrds)


#ds right 
dswR = Expand_genes(dsstallswtR)
dssR = Expand_genes(dsstallssen1R)
dsbR = Expand_genes(dsstallsdbl8R)
dsdsR = Expand_genes(dsstallsdsR)

rdswR = Expand_genes(dsstallsrwtR)
rdssR = Expand_genes(dsstallsrsen1R)
rdsbR = Expand_genes(dsstallsrdbl8R)
rdsdsR = Expand_genes(dsstallsrdsR)

#ds left
dswL = Expand_genes(dsstallswtL)
dssL = Expand_genes(dsstallssen1L)
dsbL = Expand_genes(dsstallsdbl8L)
dsdsL = Expand_genes(dsstallsdsL)

rdswL = Expand_genes(dsstallsrwtL)
rdssL = Expand_genes(dsstallsrsen1L)
rdsbL = Expand_genes(dsstallsrdbl8L)
rdsdsL = Expand_genes(dsstallsrdsL)

#selection 
sw = Expand_genes(selbodywt)
ss = Expand_genes(selbodysen1)
sb = Expand_genes(selbodydbl8)
sds = Expand_genes(selbodyds)

rsw = Expand_genes(selbodyrwt)
rss = Expand_genes(selbodyrsen1)
rsb = Expand_genes(selbodyrdbl8)
rsds = Expand_genes(selbodyrds)


def concati(rev, forward):
    forwards = np.stack(forward, axis=0 )
    revs = np.stack(rev, axis=0 )
    rev1 = revs[:, ::-1]
    
    new = np.concatenate([forwards,rev1])
    return new

#reverse LEFT, forward right 
HTw = concati(rdswL, dswR)
HTs = concati(rdssL, dssR)
HTb = concati(rdsbL, dsbR)
HTds = concati (rdsdsL,dsdsR)

#reverse right + forward leeft 
HHw = concati(rdswR, dswL)
HHs = concati(rdssR, dssL)
HHb = concati(rdsbR, dsbL)
HHds = concati(rdsdsR, dsdsL)

sa = concati(raw, aw)
sas = concati(ras,ass)
sab = concati(rad,ad)
sads = concati(rads,ads)

dbl8w = concati(rdbw, dbw)
dbl8s = concati(rdbs, dbss)
dbl8b = concati(rdbd, dbd)
dbl8ds = concati(rdbds, dbds)

caw = concati(cwr, cw)
cas = concati(csr, cs)
cab = concati(cbr, cb)
cads = concati(cdsr, cds)

sew = concati(rsw, sw)
ses = concati(rss, ss)
seb = concati(rsb, sb)
seds = concati(rsds, sds)




#%%

import scipy.stats as stats
from scipy.stats import t

x = pd.DataFrame(sa)
m = x[0].mean()

s = x[0].std() 
dof = len(x[0])-1 
confidence = 0.95
t_crit = np.abs(t.ppf((1-confidence)/2,dof))
print(t_crit)
(m-s*t_crit/np.sqrt(len(x)), m+s*t_crit/np.sqrt(len(x))) 


def conint(array):
    confidence = 0.95
    trythis = []
    x = pd.DataFrame(array)
    #print(x)
    for column in x:

        m = (x[column]).mean()
        s = x[column].std()

        dof = len(x[column])-1 
        t_crit = np.abs(t.ppf((1-confidence)/2,dof))
        interval = (m-s*t_crit/np.sqrt(len(x)), m+s*t_crit/np.sqrt(len(x)))
        
        trythis.append(interval)
        saint = pd.DataFrame(trythis, columns=['lower', 'upper'])
        #oft = saint.T
    return saint

saCON = conint(sa)
sasCON = conint(sas)
sabCON = conint(sab)
sadsCON = conint(sads)

bwCON = conint(dbl8w)
bsCON = conint(dbl8s)
bbCON = conint(dbl8b)
bdsCON = conint(dbl8ds)

HHwCON = conint(HHw)
HHsCON = conint(HHs)
HHbCON = conint(HHb)
HHdsCON = conint(HHds)


HTwCON = conint(HTw)
HTsCON = conint(HTs)
HTbCON = conint(HTb)
HTdsCON = conint(HTds)

cawCON = conint(caw)
casCON = conint(cas)
cabCON = conint(cab)
cadsCON = conint(cads)

sewCON = conint(sew)
sesCON = conint(ses)
sebCON = conint(seb)
sedsCON = conint(seds)


        
    
    
def pileup (thing):
    #thing = -thing
    #print(thing)
    #print(-thing)
    whack = np.stack( thing, axis=0 )
    stuff = whack.mean(axis=0)
    wow = pd.DataFrame(stuff)
    
    return wow
#sen compiled


#hh ht 
HHwline = pileup(HHw)
HHsline = pileup(HHs)
HHbline = pileup(HHb)
HHdsline = pileup(HHds)
HTwline = pileup(HTw)
HTsline = pileup(HTs)
HTbline = pileup(HTb)
HTdsline = pileup(HTds)

#selection
swline = pileup(sew)
ssline = pileup(ses)
sbline = pileup(seb)
sdsline = pileup(seds)


#control gene body
cwline = pileup(caw)
csline = pileup(cas)
cbline = pileup(cab)
cdsline = pileup(cads)

rcwline = pileup(cwr)
rcsline = pileup(csr)
rcbline = pileup(cbr)
rcdsline = pileup(cdsr)


#sen
awline = pileup(sa)
assline = pileup(sas)
adline = pileup(sab)
adsline= pileup(sads)


#dbl8
dbwline = pileup(dbl8w)
dbssline = pileup(dbl8s)
dbdline = pileup(dbl8b)
dbdsline= pileup(dbl8ds)


import seaborn as sns
fxdsl,((ax1, ax3, ax5,ax7)) = plt.subplots(4,1)  
fxdsl.tight_layout()
cbar_ax = fxdsl.add_axes([.91, .3, .03, .4])
sns.heatmap(HHw, cmap = 'coolwarm', ax=ax1,vmin= 0, vmax=0.5, cbar=True, cbar_ax=cbar_ax)
ax1.axes.yaxis.set_ticks([])
ax1.axes.xaxis.set_ticks([])
ax1.set_title('S∆D∆ Head-Head collisions')
ax1.set_ylabel('WT')
sns.heatmap(HHs, cmap = 'coolwarm', ax=ax3,vmin= 0, vmax=0.5, cbar=True, cbar_ax=cbar_ax)
ax3.axes.yaxis.set_ticks([])
ax3.set_ylabel('Sen1∆')
ax3.axes.xaxis.set_ticks([])
sns.heatmap(HHb, cmap = 'coolwarm', ax=ax5, vmin= 0, vmax=0.5,cbar=True, cbar_ax=cbar_ax)
ax5.axes.yaxis.set_ticks([])
ax5.set_ylabel('Dbl8∆')
ax5.axes.xaxis.set_ticks([])
sns.heatmap(HHds, cmap = 'coolwarm', ax=ax7, vmin= 0, vmax=0.5, cbar=True, cbar_ax=cbar_ax)
ax7.axes.yaxis.set_ticks([])
ax7.set_ylabel('Sen1∆Dbl8∆')
ax7.set_xticks([0,1000])
ax7.set_xticklabels(['TSS','TES'])

fuckkkkkk, (ax2, ax4, ax6, ax8) = plt.subplots(4,1)
fuckkkkkk.tight_layout()
cbar_ax = fuckkkkkk.add_axes([.91, .3, .03, .4])
sns.heatmap(HTw, cmap = 'coolwarm', ax=ax2,vmin= 0, vmax=0.5,  cbar=True, cbar_ax=cbar_ax)
ax2.set_ylabel('WT')
ax2.axes.yaxis.set_ticks([])
ax2.set_title('S∆D∆ co-directional collisions')
ax2.axes.xaxis.set_ticks([])
sns.heatmap(HTs, cmap = 'coolwarm', ax=ax4, vmin= 0, vmax=0.5, cbar=True, cbar_ax=cbar_ax)
ax4.axes.yaxis.set_ticks([])
ax4.set_ylabel('Sen1∆')
ax4.axes.xaxis.set_ticks([])
sns.heatmap(HTb, cmap = 'coolwarm', ax=ax6, vmin= 0, vmax=0.5,cbar=True, cbar_ax=cbar_ax)
ax6.axes.yaxis.set_ticks([])
ax6.set_ylabel('Dbl8∆')  
ax6.axes.xaxis.set_ticks([])
sns.heatmap(HTds, cmap = 'coolwarm', ax=ax8,  vmin= 0, vmax=0.5,cbar=True, cbar_ax=cbar_ax)
ax8.axes.yaxis.set_ticks([])
ax8.set_ylabel('Sen1∆Dbl8∆')

ax8.set_xticks([0,1000])
ax8.set_xticklabels(['TSS','TES'])


#CONTROL HEAT MAP
fix, ((ax1, ax3, ax5,ax7)) = plt.subplots(4,1)  
fix.tight_layout()
cbar_ax = fix.add_axes([.91, .3, .03, .4])
sns.heatmap(caw, cmap = 'coolwarm', ax=ax1,cbar=True, cbar_ax=cbar_ax)
ax1.axes.yaxis.set_ticks([])
ax1.axes.xaxis.set_ticks([])
ax1.set_title('Forward Strand')
ax1.set_ylabel('WT')
sns.heatmap(cas, cmap = 'coolwarm', ax=ax3,cbar=True, cbar_ax=cbar_ax)

ax3.axes.yaxis.set_ticks([])
ax3.set_ylabel('Sen1∆')
ax3.axes.xaxis.set_ticks([])
sns.heatmap(cab, cmap = 'coolwarm', ax=ax5, cbar=True, cbar_ax= cbar_ax)
ax5.axes.yaxis.set_ticks([])
ax5.set_ylabel('Dbl8∆')
ax5.axes.xaxis.set_ticks([])
    
sns.heatmap(cads, cmap = 'coolwarm', ax=ax7,cbar=True, cbar_ax=cbar_ax)
ax7.axes.yaxis.set_ticks([])
ax7.set_ylabel('Sen1∆Dbl8∆')

ax7.set_xticks([0,1000])
ax7.set_xticklabels(['TSS','TES'])


#DBL8 heatmap !!!!!
fxb, ((ax1, ax3, ax5,ax7)) = plt.subplots(4,1)  
fxb.tight_layout()
cbar_ax = fxb.add_axes([.91, .3, .03, .4])
sns.heatmap(dbl8w, cmap = 'coolwarm', ax=ax1,  cbar=True, cbar_ax=cbar_ax)
ax1.axes.yaxis.set_ticks([])
ax1.axes.xaxis.set_ticks([])
ax1.set_title('Forward Strand')
ax1.set_ylabel('WT')
sns.heatmap(dbl8s, cmap = 'coolwarm', ax=ax3,  cbar=True, cbar_ax=cbar_ax)
ax3.axes.yaxis.set_ticks([])
ax3.set_ylabel('Sen1∆')
ax3.axes.xaxis.set_ticks([])
sns.heatmap(dbl8b, cmap = 'coolwarm', ax=ax5,  cbar=True, cbar_ax=cbar_ax)
ax5.axes.yaxis.set_ticks([])
ax5.set_ylabel('Dbl8∆')
ax5.axes.xaxis.set_ticks([])
sns.heatmap(dbl8ds, cmap = 'coolwarm', ax=ax7,  cbar=True, cbar_ax=cbar_ax)
ax7.axes.yaxis.set_ticks([])
ax7.set_ylabel('Sen1∆Dbl8∆')

ax7.set_xticks([0,1000])
ax7.set_xticklabels(['TSS','TES'])

#SEN1 heatmap bitches 
fx, ((ax1, ax3, ax5,ax7)) = plt.subplots(4,1)  
fx.tight_layout()
cbar_ax = fx.add_axes([.91, .3, .03, .4])
sns.heatmap(sa, cmap = 'coolwarm', ax=ax1, cbar=True, cbar_ax=cbar_ax)
ax1.axes.yaxis.set_ticks([])
ax1.axes.xaxis.set_ticks([])
ax1.set_title('Forward Strand')
ax1.set_ylabel('WT')
sns.heatmap(sas, cmap = 'coolwarm', ax=ax3, cbar=True, cbar_ax=cbar_ax)
ax3.axes.yaxis.set_ticks([])
ax3.set_ylabel('Sen1∆')
ax3.axes.xaxis.set_ticks([])
sns.heatmap(sab, cmap = 'coolwarm', ax=ax5,  cbar=True, cbar_ax=cbar_ax)
ax5.axes.yaxis.set_ticks([])
ax5.set_ylabel('Dbl8∆')
ax5.axes.xaxis.set_ticks([])
sns.heatmap(sads, cmap = 'coolwarm', ax=ax7, cbar=True, cbar_ax=cbar_ax)
ax7.axes.yaxis.set_ticks([])
ax7.set_ylabel('Sen1∆Dbl8∆')


ax7.set_xticks([0,1000])
ax7.set_xticklabels(['TSS','TES'])



  
    
def pileup (thing):
    #thing = -thing
    #print(thing)
    #print(-thing)
    whack = np.stack( thing, axis=0 )
    stuff = whack.mean(axis=0)

    wow = pd.DataFrame(stuff, columns=['OriginalData'])

    # Apply smoothing using a rolling window
    wow['SmoothedData'] = wow['OriginalData'].rolling(window=50, min_periods=1).mean()
    wowe = wow[['SmoothedData']]

    return wowe
#sen compiled


#hh ht 
HHwline = pileup(HHw)
HHsline = pileup(HHs)
HHbline = pileup(HHb)
HHdsline = pileup(HHds)
HTwline = pileup(HTw)
HTsline = pileup(HTs)
HTbline = pileup(HTb)
HTdsline = pileup(HTds)

#selection
swline = pileup(sew)
ssline = pileup(ses)
sbline = pileup(seb)
sdsline = pileup(seds)


#control gene body
cwline = pileup(caw)
csline = pileup(cas)
cbline = pileup(cab)
cdsline = pileup(cads)

rcwline = pileup(cwr)
rcsline = pileup(csr)
rcbline = pileup(cbr)
rcdsline = pileup(cdsr)


#sen
awline = pileup(sa)
assline = pileup(sas)
adline = pileup(sab)
adsline= pileup(sads)


#dbl8
dbwline = pileup(dbl8w)
dbssline = pileup(dbl8s)
dbdline = pileup(dbl8b)
dbdsline= pileup(dbl8ds)

#%%
figgy, (ax1) =plt.subplots(1, sharey=True)

x = np.arange(0, 1000, 1)
y11 = np.array(cawCON['lower']).flatten()
y10= np.array(cawCON['upper']).flatten()
ax1.fill_between(x, y10, y11, facecolor='grey', alpha =0.4,)
y1s =np.array(sasCON['lower']).flatten()
y2s=np.array(sasCON['upper']).flatten()
ax1.fill_between(x, y1s, y2s, facecolor='deepskyblue', alpha =0.4, )
y3b = np.array(bbCON['lower']).flatten()
y4b = np.array(bbCON['upper']).flatten()
ax1.fill_between(x, y3b, y4b, facecolor='gold', alpha =0.4)
y5ds = np.array(HHdsCON['lower']).flatten()
y6ds= np.array(HHdsCON['upper']).flatten()
ax1.fill_between(x, y5ds, y6ds, facecolor='indianred', alpha =0.4)
y7ds = np.array(HTdsCON['lower']).flatten()
y8ds= np.array(HTdsCON['upper']).flatten()
ax1.fill_between(x, y7ds, y8ds, facecolor='salmon', alpha =0.4, )
ax1.legend(loc='best')


ax1.set_xticks([0,1000])
ax1.set_xticklabels(['TSS','TES'])
ax1.legend(loc='best')
ax1.set_ylabel('normalised reads')






figgy, (ax1) =plt.subplots(1, sharey=True)

x = np.arange(0, 1000, 1)
ax1.plot(cwline.index, cwline['SmoothedData'], color = 'black', alpha =0.5, linewidth=1)
ax1.plot(csline.index, csline['SmoothedData'], color = 'deepskyblue', alpha=0.5, linewidth=1)
ax1.plot(cbline.index, cbline['SmoothedData'], color = 'orange', alpha=0.5, linewidth=1)
ax1.plot(cdsline.index, cdsline['SmoothedData'], color = 'tomato', alpha=0.5, linewidth=1)
y11 = np.array(cawCON['lower']).flatten()
y10= np.array(cawCON['upper']).flatten()
ax1.fill_between(x, y10, y11, facecolor='grey', alpha =0.3,)
y11s = np.array(casCON['lower']).flatten()
y10s= np.array(casCON['upper']).flatten()
#ax1.fill_between(x, y10s, y11s, facecolor='mediumaquamarine', alpha =0.3)
y11b = np.array(cabCON['lower']).flatten()
y10b = np.array(cabCON['upper']).flatten()
#ax1.fill_between(x, y10b, y11b, facecolor='mediumaquamarine', alpha =0.3)
y11ds = np.array(cadsCON['lower']).flatten()
y10ds= np.array(cadsCON['upper']).flatten()
#ax1.fill_between(x, y10ds, y11ds, facecolor='mediumaquamarine', alpha =0.3)
ax1.legend(loc='best')


ax1.set_xticks([0,1000])
ax1.set_xticklabels(['TSS','TES'])
ax1.legend(loc='best')
ax1.set_ylabel('normalised reads')



figgy, (ax1) =plt.subplots(1, sharey=True)
ax1.plot(awline.index, awline['SmoothedData'], color = 'black', alpha =0.8, linewidth=1)
ax1.plot(assline.index, assline['SmoothedData'], color = 'deepskyblue', alpha=0.8, linewidth=1)
ax1.plot(adline.index, adline['SmoothedData'], color = 'orange', alpha=0.8, linewidth=1)
ax1.plot(adsline.index, adsline['SmoothedData'], color = 'tomato', alpha=0.8, linewidth=1)
x = np.arange(0, 1000, 1)
y1 = np.array(saCON['lower']).flatten()
y2= np.array(saCON['upper']).flatten()
ax1.fill_between(x, y1, y2, facecolor='grey', alpha =0.2)
y1s =np.array(sasCON['lower']).flatten()
y2s=np.array(sasCON['upper']).flatten()
ax1.fill_between(x, y1s, y2s, facecolor='deepskyblue', alpha =0.2, )
y1b =np.array(sabCON['lower']).flatten()
y2b=np.array(sabCON['upper']).flatten()
#ax1.fill_between(x, y1b, y2b, facecolor='gold', alpha =0.2)
y1ds =np.array(sadsCON['lower']).flatten()
y2ds=np.array(sadsCON['upper']).flatten()
ax1.fill_between(x, y1ds, y2ds, facecolor='indianred', alpha =0.1)
ax1.set_xticks([0,1000])
ax1.set_xticklabels(['TSS','TES'])
#ax1.set_ylim([-4, 124])
ax1.set_ylabel('normalised reads')
ax1.legend(loc='best')

ax1.set_xticks([0,1000])
ax1.set_xticklabels(['TSS','TES'])
ax1.legend(loc='best')
ax1.set_ylabel('normalised reads')
#ax1.fill_between(x, y1a, y2a, facecolor='mediumaquamarine', alpha =0.2)
#y2 = linerev
#ax2.fill_between(x, y2, facecolor='green', alpha =0.3)
#filkggy, (ax1) =plt.subplots(1,1, sharey=True)

figgy, (ax1) =plt.subplots(1,1, sharey=True)
ax1.plot(dbwline.index, dbwline['SmoothedData'], color = 'black', alpha =0.8, linewidth=1)
ax1.plot(dbssline.index, dbssline['SmoothedData'], color = 'deepskyblue', alpha=0.8, linewidth=1)
ax1.plot(dbdline.index, dbdline['SmoothedData'], color = 'orange', alpha=0.8, linewidth=1)
ax1.plot(dbdsline.index, dbdsline['SmoothedData'], color = 'tomato', alpha=0.8, linewidth=1)
y3 = np.array(bwCON['lower']).flatten()
y4 = np.array(bwCON['upper']).flatten()
ax1.fill_between(x, y3, y4, facecolor='grey', alpha =0.2)
y3s = np.array(bsCON['lower']).flatten()
y4s = np.array(bsCON['upper']).flatten()
#ax1.fill_between(x, y3s, y4s, facecolor='deepskyblue', alpha =0.1,)
y3b = np.array(bbCON['lower']).flatten()
y4b = np.array(bbCON['upper']).flatten()
ax1.fill_between(x, y3b, y4b, facecolor='gold', alpha =0.2)
y3ds = np.array(bdsCON['lower']).flatten()
y4ds = np.array(bdsCON['upper']).flatten()
ax1.fill_between(x, y3ds, y4ds, facecolor='indianred', alpha =0.1)
ax1.legend(loc='best')
ax1.set_xticks([0,1000])
ax1.set_xticklabels(['TSS','TES'])
#ax1.set_ylim([-4, 124])
ax1.legend(loc='best')
ax1.set_ylabel('normalised reads')
ax1.set_xticks([0,1000])
ax1.set_xticklabels(['TSS','TES'])
ax1.legend(loc='best')
ax1.set_ylabel('normalised reads')


figgy, (ax1) =plt.subplots(1, sharey=True)
#ds left
ax1.plot(HHwline.index, HHwline['SmoothedData'], color = 'black', alpha =0.8, linewidth=1)
ax1.plot(HHsline.index, HHsline['SmoothedData'], color = 'deepskyblue', alpha=0.8, linewidth=1)
ax1.plot(HHbline.index, HHbline['SmoothedData'], color = 'orange', alpha=0.8, linewidth=1)
ax1.plot(HHdsline.index, HHdsline['SmoothedData'], color = 'tomato', alpha=0.8, linewidth=1)
y5 = np.array(HHwCON['lower']).flatten()
y6= np.array(HHwCON['upper']).flatten()
ax1.fill_between(x, y5, y6, facecolor='grey', alpha =0.2)
y5s = np.array(HHsCON['lower']).flatten()
y6s= np.array(HHsCON['upper']).flatten()
ax1.fill_between(x, y5s, y6s, facecolor='deepskyblue', alpha =0.2)
y5b = np.array(HHbCON['lower']).flatten()
y6b= np.array(HHbCON['upper']).flatten()
ax1.fill_between(x, y5b, y6b, facecolor='gold', alpha =0.2)
y5ds = np.array(HHdsCON['lower']).flatten()
y6ds= np.array(HHdsCON['upper']).flatten()
ax1.fill_between(x, y5ds, y6ds, facecolor='indianred', alpha =0.2)
#ax1.legend(loc='best')


ax1.set_xticks([0,1000])
ax1.set_xticklabels(['TSS','TES'])
ax1.legend(loc='best')
ax1.set_ylabel('normalised reads')

figgy, (ax1) =plt.subplots(1, sharey=True)
#ds right
ax1.plot(HTwline.index, HTwline['SmoothedData'], color = 'black', alpha =0.8, linewidth=1)
ax1.plot(HTsline.index, HTsline['SmoothedData'], color = 'deepskyblue', alpha=0.8, linewidth=1)
ax1.plot(HTbline.index, HTbline['SmoothedData'], color = 'orange', alpha=0.8, linewidth=1)
ax1.plot(HTdsline.index, HTdsline['SmoothedData'], color = 'tomato', alpha=0.8, linewidth=1)
y7 = np.array(HTwCON['lower']).flatten()
y8= np.array(HTwCON['upper']).flatten()
ax1.fill_between(x, y7, y8, facecolor='grey', alpha =0.2)
y7s = np.array(HTsCON['lower']).flatten()
y8s= np.array(HTsCON['upper']).flatten()
ax1.fill_between(x, y7s, y8s, facecolor='deepskyblue', alpha =0.2)
y7b = np.array(HTbCON['lower']).flatten()
y8b= np.array(HTbCON['upper']).flatten()
ax1.fill_between(x, y7b, y8b, facecolor='gold', alpha =0.2)
y7ds = np.array(HTdsCON['lower']).flatten()
y8ds= np.array(HTdsCON['upper']).flatten()
ax1.fill_between(x, y7ds, y8ds, facecolor='salmon', alpha =0.2, )
#ax1.legend(loc='best')

ax1.set_xticks([0,1000])
ax1.set_xticklabels(['TSS','TES'])
ax1.legend(loc='best')
ax1.set_ylabel('normalised reads')
#ax1.set_ylim([-4, 124])


figgy, (ax1) =plt.subplots(1, sharey=True)
ax1.plot(cwline.index, cwline['SmoothedData'], color = 'black', alpha =0.5, linewidth=1)
ax1.plot(csline.index, csline['SmoothedData'], color = 'deepskyblue', alpha=0.5, linewidth=1)
ax1.plot(cbline.index, cbline['SmoothedData'], color = 'orange', alpha=0.5, linewidth=1)
ax1.plot(cdsline.index, cdsline['SmoothedData'], color = 'tomato', alpha=0.5, linewidth=1)
y11 = np.array(cawCON['lower']).flatten()
y10= np.array(cawCON['upper']).flatten()
ax1.fill_between(x, y10, y11, facecolor='grey', alpha =0.3,)
y11s = np.array(casCON['lower']).flatten()
y10s= np.array(casCON['upper']).flatten()
#ax1.fill_between(x, y10s, y11s, facecolor='mediumaquamarine', alpha =0.3)
y11b = np.array(cabCON['lower']).flatten()
y10b = np.array(cabCON['upper']).flatten()
#ax1.fill_between(x, y10b, y11b, facecolor='mediumaquamarine', alpha =0.3)
y11ds = np.array(cadsCON['lower']).flatten()
y10ds= np.array(cadsCON['upper']).flatten()
ax1.fill_between(x, y10ds, y11ds, facecolor='mediumaquamarine', alpha =0.3)
ax1.legend(loc='best')



ax1.set_xticks([0,1000])
ax1.set_xticklabels(['TSS','TES'])
ax1.legend(loc='best')
ax1.set_ylabel('normalised reads')
#ax1.set_title('forward strand')


figgy, (ax1) =plt.subplots(1, sharey=True)
ax1.plot(swline.index, swline[0], color = 'black', alpha =0.5, linewidth=1)
ax1.plot(ssline.index, ssline[0], color = 'deepskyblue', alpha=0.5, linewidth=1)
ax1.plot(sbline.index, sbline[0], color = 'orange', alpha=0.5, linewidth=1)
ax1.plot(sdsline.index, sdsline[0], color = 'tomato', alpha=0.5, linewidth=1)
y9 = np.array(sewCON['lower']).flatten()
y10= np.array(sewCON['upper']).flatten()
ax1.fill_between(x, y9, y10, facecolor='grey', alpha =0.2)



y9s = np.array(sesCON['lower']).flatten()
y10s= np.array(sesCON['upper']).flatten()
ax1.fill_between(x, y9s, y10s, facecolor='deepskyblue', alpha =0.3)
y9b = np.array(sebCON['lower']).flatten()
y10b= np.array(sebCON['upper']).flatten()
ax1.fill_between(x, y9b, y10b, facecolor='orange', alpha =0.3)
y9ds = np.array(sedsCON['lower']).flatten()
y10ds= np.array(sedsCON['upper']).flatten()
ax1.fill_between(x, y9ds, y10ds, facecolor='tomato', alpha =0.3)
ax1.legend(loc='best')

ax1.set_xticks([0,1000])
ax1.set_xticklabels(['TSS','TES'])
ax1.legend(loc='best')
ax1.set_ylabel('normalised reads')


ax1.set_ylim([-4, 124])
#%%%

# %% Importing genome sequence and calculating AT content
from pysam import FastaFile
data_fasta = FastaFile("/Users/joannafernandez/Desktop/python/Schizosaccharomyces_pombe_all_chromosomes (1).fa")
data_AT_content = pd.DataFrame()
binwidth = 30

for chro in ["chr1", "chr2", "chr3"] :

    tmp_data = pd.DataFrame({"chro" : chro,
                             "pos" : range(len(data_fasta.fetch(chro))),
                             "seq" : [*data_fasta.fetch(chro)]})
   
    tmp_data["pos"] = tmp_data["pos"] + 1
    tmp_data["AT_score"] = np.where(tmp_data["seq"].isin(["A", "a", "T", "t"]), 1, 0)
    tmp_bins = list(np.repeat(list(range((int(binwidth / 2)), len(tmp_data) + (int(binwidth / 2)), binwidth)), binwidth))[0:len(tmp_data)]
    tmp_data["bin"] = tmp_bins
    data_AT_content = pd.concat([data_AT_content, tmp_data])
    del([tmp_data, tmp_bins])
   
data_AT_content_binned = data_AT_content.groupby(
    ["chro", "bin"], as_index = False).agg(
        AT_bias = ("AT_score", lambda x : sum(x) / binwidth))
       
data_AT_content_binned = data_AT_content_binned.rename(columns = {"bin" : "pos"})

data_AT_content_binned["AT_bias_smooth"] = data_AT_content_binned.groupby(
    ["chro"], as_index = False)["AT_bias"].transform(
        lambda x : x.rolling(12, center = True).mean())

#%%


from Bio import SeqIO

genome_file = "/Users/joannafernandez/Desktop/python/Schizosaccharomyces_pombe_all_chromosomes (1).fa"

# Load the S. pombe genome from the FASTA file
genome_records = list(SeqIO.parse(genome_file, "fasta"))

print(len(genome_records))

# Access individual records
for record in genome_records:
    print(f"ID: {record.id}")
    print(f"Sequence length: {len(record)}")

# Check the number of records loaded (should be greater than 1 for multiple sequences)
print(len(genome_records))

first_genome_sequence = genome_records[0]

# Accessing specific locations in the genome sequence
position_1 = 1000
position_2 = 2000

# Retrieve nucleotides at the specified positions
nucleotide_at_position_1 = first_genome_sequence.seq[position_1 - 1]  # Biopython uses 0-based indexing
nucleotide_at_position_2 = first_genome_sequence.seq[position_2 - 1]

# Print the results
print(f"Nucleotide at position {position_1}: {nucleotide_at_position_1}")
print(f"Nucleotide at position {position_2}: {nucleotide_at_position_2}")



#try this

genome_file = "/Users/joannafernandez/Desktop/python/Schizosaccharomyces_pombe_all_chromosomes (1).fa"

# Load the S. pombe genome from the FASTA file
genome_records = list(SeqIO.parse(genome_file, "fasta"))

# Initialize lists to store data
chro_list = []
position_list = []
base_list = []

# Iterate through genome records
for record in genome_records:
    chro_id = record.id
    sequence = record.seq

    # Iterate through the sequence and store position and base
    for position, base in enumerate(sequence, start=1):
        chro_list.append(chro_id)
        position_list.append(position)
        base_list.append(base)

# Create a DataFrame
genome_df = pd.DataFrame({'Chr': chro_list, 'Pos': position_list, 'base': base_list})

# Display the DataFrame
print(genome_df)

unique_values = genome_df['base'].unique()

genome_df.loc[genome_df['Chr'] == 'I', 'Chr'] = 1
genome_df.loc[genome_df['Chr'] == 'II', 'Chr'] = 2
genome_df.loc[genome_df['Chr'] == 'III', 'Chr'] = 3




#VERSION 3 INCLUDE REVERSE STRAND GENES and ttake normscore
def flank(file,frame1,gd):
    data_output = pd.DataFrame()
   # data_output_cnorm = pd.DataFrame()
    
    for i in range(len(file)):
        tempstrand = file.iloc[i]['coding_strand']
        tempchro = file.iloc[i]["chro"]
        if tempstrand == 'forward':
            if 
            tempstart = file.iloc[i]["start"]
            tempchro = file.iloc[i]["chro"]

            tempend = file.iloc[i]["end"]


            tempsubsetfor = frame1.loc[(frame1['Pos']>= tempstart) & (frame1['Pos'] <= tempend) & (frame1['Chr'] == tempchro)]

            tempsubsetfor = tempsubsetfor.reset_index(drop = True)
            line = pd.DataFrame({"Pos": tempstarty}, index=[6000])
            line2 = pd.DataFrame({"Pos": startminus}, index=[0])
            tempsubsetfor = tempsubsetfor.append(line, ignore_index=False)
            tempsubsetfor = tempsubsetfor.append(line2, ignore_index=False)
            
            tempsubsetfor = tempsubsetfor.sort_values('Pos', ascending=True)


#this fills empty with 0
            for x in range(len(tempsubsetfor)):
                pos = tempsubsetfor.iloc[x]['Pos']

                tpos = tempsubsetfor.iloc[0]['Pos']

                tempsubsetfor.loc[tempsubsetfor.index[x], 'aPos'] = (pos - tpos)

            tempsubsetfor = tempsubsetfor.set_index('aPos')
                        

            data_blank = pd.DataFrame({
    
                "aPos" : np.arange(0.0, 6000, 1)}).merge(
        
                    tempsubsetfor, on = "aPos", how = "left"

                    )

            data_blank = data_blank[["aPos", "normscore"]]
            
            data_output = pd.concat([data_output, data_blank['normscore']], axis =1, ignore_index=True)
         #   data_output_cnorm = pd.concat([data_output_cnorm, data_blank['Cnorm']], axis =1, ignore_index=True)
            
            data_output = data_output.fillna(0)

         #   data_output_cnorm = data_output_cnorm.fillna(0)
        
        
        if tempstrand == 'reverse':
            tempstartr = file.iloc[i]["end"]            

            startminusr = (tempstartr) + 3000
            tempstartyr = (tempstartr) - 3000

        
            tempsubsetrev = frame1.loc[(frame1['Pos']>= tempstartyr) & (frame1['Pos'] <= startminusr) & (frame1['Chr'] == tempchro)]

            tempsubsetrev = tempsubsetrev.reset_index(drop = True)
            line = pd.DataFrame({"Pos": tempstartyr}, index=[6000])
            line2 = pd.DataFrame({"Pos": startminusr}, index=[0])
            tempsubsetrev = tempsubsetrev.append(line, ignore_index=False)
            tempsubsetrev = tempsubsetrev.append(line2, ignore_index=False)
            
            tempsubsetrev = tempsubsetrev.sort_values('Pos', ascending=True)

            for x in range(len(tempsubsetrev)):
                pos = tempsubsetrev.iloc[0]['Pos']

                tpos = tempsubsetrev.iloc[x]['Pos']

                tempsubsetrev.loc[tempsubsetrev.index[x], 'aPos'] = (pos - tpos)

            tempsubsetrev = tempsubsetrev.set_index('aPos')
                        

            data_blankr = pd.DataFrame({
    
                "aPos" : np.arange(0.0, 6000, 1)}).merge(
        
                    tempsubsetrev, on = "aPos", how = "left"

                    )

            data_blankr = data_blankr[["aPos", "normscore"]]
            
            data_output = pd.concat([data_output, data_blankr['normscore']], axis =1, ignore_index=True)
           # data_output_cnorm = pd.concat([data_output_cnorm, data_blankr['Cnorm']], axis =1, ignore_index=True)
            
            data_output = data_output.fillna(0)

           # data_output_cnorm = data_output_cnorm.fillna(0)

    return data_output

wtforflanksenW = flank(sen1stall, wt)
senforflanksenW = flank(sen1stall, sen1)
dblforflanksenW = flank(sen1stall, dbl8)
dsforflanksenW = flank(sen1stall, ds)




#%%%


import scipy.stats as stats
from scipy.stats import t

def conint(array):
    confidence = 0.95
    trythis = []
    reshaped_df = array.transpose()
    x = pd.DataFrame(reshaped_df)
    #print(x)
    for row in x:

        m = (x[row]).mean()
        s = x[row].std()

        dof = len(x[row])-1 
        t_crit = np.abs(t.ppf((1-confidence)/2,dof))
        interval = (m-s*t_crit/np.sqrt(len(x)), m+s*t_crit/np.sqrt(len(x)))
        
        trythis.append(interval)
        saint = pd.DataFrame(trythis, columns=['lower', 'upper'])
        saint['lower'] = saint['lower'].rolling(50).mean()
        saint['upper'] = saint['upper'].rolling(50).mean()
    
        #oft = saint.T
    return saint

senWw = conint(wtforflanksenW)
#senCw = conint(wtforflanksenC)
senWs = conint(senforflanksenW)
senWb = conint(dblforflanksenW)
senWd = conint(dsforflanksenW)

dblWw = conint(wtforflankdblW)
#dblCw = conint(wtforflankdblC)
dblWb = conint(dblforflankdblW)
dblWs = conint(senforflankdblW)
dblWd = conint(dsforflankdblW)



#dblC = conint(dblforflankdblC)

dsWw = conint(wtforflankdsW)
#dsCw = conint(wtforflankdsC)
dsWs = conint(dsforflankdsW)
dsWb = conint(dblforflankdsW)
dsWd = conint(dsforflankdsW)

#dsC = conint(dsforflankdsC)

conW = conint(wtforflankconW)
conWs = conint(senforflankconW)
conWb = conint(dblforflankconW)
conWd = conint(dsforflankconW)
#conC = conint(wtforflankconC)

sel10W = conint(wtforflanksel10W)
sel10C = conint(wtforflanksel10C)
sel200W = conint(wtforflanksel200W)
sel200C = conint(wtforflanksel200C)


####okay so i need to redo the plots for everything again because idiot forgot something llower upper here shoulld be flipped?
##fix it now, and redolater bitch
'''
figgy, (ax1) =plt.subplots(1, sharey=True)
x = np.arange(0, 1101, 1)
y9 = np.array(conC['upper']).flatten()
y10= np.array(conW['lower']).flatten()
ax1.fill_between(x, y9, y10, facecolor='grey', alpha =0.5)

y9s = np.array(senW['upper']).flatten()
y10s= np.array(senC['lower']).flatten()
ax1.fill_between(x, y9s, y10s, facecolor='deepskyblue', alpha =0.5)
y9b = np.array(dblW['upper']).flatten()
y10b= np.array(dblC['lower']).flatten()
ax1.fill_between(x, y9b, y10b, facecolor='orange', alpha =0.5)
y9ds = np.array(dsW['upper']).flatten()
y10ds= np.array(dsC['lower']).flatten()
ax1.fill_between(x, y9ds, y10ds, facecolor='tomato', alpha =0.5)
ax1.set_xticks([0,100,1100])
#ax1.set_ylim(-0.012, 0.012)
ax1.set_xticklabels(['+0.1Kb','TSS','-1kb'])
'''




#ax1.legend(loc='best')
##have not fixed any other code except above ^
#what i should do is pllot upper and lower for both W and C 
fixed, (ax1) =plt.subplots(1, sharey=True)
x = np.arange(0, 6001, 1)
ax1.set_xticks([0,1000,6000])
#ax1.set_ylim(-0.02, 0.04)
ax1.set_xticklabels(['+1Kb','TSS','-1kb'])
#ax1.legend(loc='best')


y9s = np.array(dsWw['upper']).flatten()
y10s= np.array(dsWw['lower']).flatten()
ax1.fill_between(x, y9s, y10s, facecolor='black', alpha =0.5)
#y9sc = np.array(senC['upper']).flatten()
#y10sc= np.array(senC['lower']).flatten()
#ax1.fill_between(x, y9sc, y10sc, facecolor='deepskyblue', alpha =0.5)

y9b = np.array(dsW['upper']).flatten()
y10b= np.array(dsW['lower']).flatten()
ax1.fill_between(x, y9b, y10b, facecolor='tomato', alpha =0.5)
#y9bc = np.array(dblC['upper']).flatten()
#y10bc= np.array(dblC['lower']).flatten()
#ax1.fill_between(x, y9bc, y10bc, facecolor='orange', alpha =0.5)

y9ds = np.array(dsW['upper']).flatten()
y10ds= np.array(dsW['lower']).flatten()
ax1.fill_between(x, y9ds, y10ds, facecolor='tomato', alpha =0.5)
#y9dsc = np.array(dsC['upper']).flatten()
#y10dcs= np.array(dsC['lower']).flatten()
#ax1.fill_between(x, y9dsc, y10dcs, facecolor='tomato', alpha =0.5)

y9 = np.array(conW['upper']).flatten()
y10= np.array(conW['lower']).flatten()
ax1.fill_between(x, y9, y10, facecolor='grey', alpha =0.5)
#y9c = np.array(conC['upper']).flatten()
#y1c0= np.array(conC['lower']).flatten()
#ax1.fill_between(x, y9c, y1c0, facecolor='grey', alpha =0.5)
ax1.set_xticks([0,1000,2000])
ax1.set_ylim(-0.02, 0.04)
ax1.set_xticklabels(['+1Kb','TSS','-1kb'])
#ax1.legend(loc='best')
##have not fixed any other code except above ^
#what i should do is pllot upper and lower for both W and C 





dblWw = conint(wtforflankdblW)
dblCw = conint(wtforflankdblC)
dblW = conint(dblforflankdblW)
dblC = conint(dblforflankdblC)

figgy, (ax1) =plt.subplots(1, sharey=True)
x = np.arange(0, 6001, 1)

y9sr = np.array(dsW['lower']).flatten()
y10sr= np.array(dsW['upper']).flatten()
ax1.fill_between(x, y9sr, y10sr, facecolor='orange', alpha =0.5)
y9sg = np.array(dsC['lower']).flatten()
y10sg= np.array(dsC['upper']).flatten()
ax1.fill_between(x, y9sg, y10sg, facecolor='orange', alpha =0.5)

y9bx = np.array(dblWw['lower']).flatten()
y10bx= np.array(dsWw['upper']).flatten()

ax1.fill_between(x, y9bx, y10bx, facecolor='black', alpha =0.5)

y9bxx = np.array(dblCw['lower']).flatten()
y10bxx= np.array(dsCw['upper']).flatten()

ax1.fill_between(x, y9bxx, y10bxx, facecolor='black', alpha =0.5)


ax1.set_xticks([0,100,1100])
#ax1.set_ylim(-0.012, 0.012)
ax1.set_xticklabels(['+0.1Kb','TSS','-1kb'])



#try upper only

figgy, (ax1) =plt.subplots(1, sharey=True)
#x = np.arange(0, 1101, 1)

y9sr = np.array(dsW['lower']).flatten()
y10sr= np.array(dsW['upper']).flatten()

y9sg = np.array(dsC['lower']).flatten()
y10sg= np.array(dsC['upper']).flatten()
ax1.fill_between(x, y10sr, y9sg, facecolor='orange', alpha =0.5)

y9bx = np.array(dblWw['lower']).flatten()
y10bx= np.array(dsWw['upper']).flatten()

y9bxx = np.array(dblCw['lower']).flatten()
y10bxx= np.array(dsCw['upper']).flatten()

ax1.fill_between(x, y10bx, y9bxx, facecolor='black', alpha =0.5)


ax1.set_xticks([0,100,1100])
#ax1.set_ylim(-0.012, 0.012)
ax1.set_xticklabels(['+0.1Kb','TSS','-1kb'])




figgy, (ax1) =plt.subplots(1, sharey=True)
ax1.set_ylim(-0.035, 0.035)
#ax1.set_ylim(-0.022, 0.022)
#ax1.set_ylim(-0.012, 0.012)
y9b1 = np.array(sel10W['lower']).flatten()
y10b1= np.array(sel10W['upper']).flatten()
ax1.fill_between(x, y9b1, y10b1, facecolor='forestgreen', alpha =0.5)
y9b1c = np.array(sel10C['lower']).flatten()
y10b1c= np.array(sel10C['upper']).flatten()
ax1.fill_between(x, y9b1c, y10b1c, facecolor='forestgreen', alpha =0.5)

y9ds2 = np.array(sel200W['lower']).flatten()
y10ds2= np.array(sel200W['upper']).flatten()
ax1.fill_between(x, y9ds2, y10ds2, facecolor='mediumpurple', alpha =0.5)
y9ds2c = np.array(sel200C['lower']).flatten()
y10ds2c= np.array(sel200C['upper']).flatten()
ax1.fill_between(x, y9ds2c, y10ds2c, facecolor='mediumpurple', alpha =0.5)
#ax1.legend(loc='best')

ax1.set_xticks([0,100,1100])
ax1.set_xticklabels(['+0.1Kb','TSS','-1kb'])
#ax1.legend(loc='best')
#ax1.set_ylabel('normalised reads')




def take_mean(all_data):
    all_data['mean'] = all_data.mean(axis=1)
    all_data['smoo_mean'] = all_data['mean'].rolling(window = 200, center=True).mean()
    
    return all_data

#wtforflanksenW = take_mean(wtforflanksenW, 100)
#wtforflanksenC = take_mean(wtforflanksenC)
wtforflanksenW = take_mean(wtforflanksenW)
senforflanksenW = take_mean(senforflanksenW)
#senforflanksenC = take_mean(senforflanksenC)
dblforflanksenW = take_mean(dblforflanksenW)
#dblforflanksenC = take_mean(dblforflanksenC)
dsforflanksenW = take_mean(dsforflanksenW)
#dsforflanksenC = take_mean(dsforflanksenC)

#wtforflankdblC = take_mean(wtforflankdblC)
wtforflankdblW = take_mean(wtforflankdblW)
senforflankdblW = take_mean(senforflankdblW)
#senforflankdblC = take_mean(senforflankdblC)
dblforflankdblW = take_mean(dblforflankdblW)
#dblforflankdblC = take_mean(dblforflankdblC)
dsforflankdblW = take_mean(dsforflankdblW)
#dsforflankdblC = take_mean(dsforflankdblC)

#wtforflankconC = take_mean(wtforflankconC)
wtforflankconW = take_mean(wtforflankconW)
senforflankconW = take_mean(senforflankconW)
#senforflankconC = take_mean(senforflankconC)
dblforflankconW = take_mean(dblforflankconW)
#dblforflankconC = take_mean(dblforflankconC)
dsforflankconW = take_mean(dsforflankconW)
#dsforflankconC = take_mean(dsforflankconC)

#wtforflankdsC = take_mean(wtforflankdsC)
wtforflankdsW = take_mean(wtforflankdsW)
senforflankdsW = take_mean(senforflankdsW)
#senforflankdsC = take_mean(senforflankdsC)
dblforflankdsW = take_mean(dblforflankdsW)
#dblforflankdsC = take_mean(dblforflankdsC)
dsforflankdsW = take_mean(dsforflankdsW)
#dsforflankdsC = take_mean(dsforflankdsC)

wtforflanksel10W = take_mean(wtforflanksel10W)
wtforflanksel10C = take_mean(wtforflanksel10C)
senforflanksel10W = take_mean(senforflanksel10W)
senforflanksel10C = take_mean(senforflanksel10C)
dblforflanksel10W = take_mean(dblforflanksel10W)
dblforflanksel10C = take_mean(dblforflanksel10C)
dsforflanksel10W = take_mean(dsforflanksel10W)
dsforflanksel10C = take_mean(dsforflanksel10C)

wtforflanksel200W = take_mean(wtforflanksel200W)
wtforflanksel200C = take_mean(wtforflanksel200C)
senforflanksel200W = take_mean(senforflanksel200W)
senforflanksel200C = take_mean(senforflanksel200C)
dblforflanksel200W = take_mean(dblforflanksel200W)
dblforflanksel200C = take_mean(dblforflanksel200C)
dsforflanksel200W = take_mean(dsforflanksel200W)
dsforflanksel200C = take_mean(dsforflanksel200C)







def plotall(ww, wc, sw, sc, bw, bc, dw,dc, title):
    fi,((ax1)) = plt.subplots(1,1, sharex=True, sharey=True) 

    ax1.plot(ww.index.values, ww['smoo_mean'], color = 'black', alpha= 0.5)
    ax1.plot(wc.index.values, wc['smoo_mean'], color = 'black', alpha= 0.5)
    ax1.set_title(title)
    ax1.set_ylabel('WT')


    ax1.plot(sw.index.values, sw['smoo_mean'], color = 'steelblue', alpha= 0.5)
    ax1.plot(sc.index.values, sc['smoo_mean'], color = 'steelblue', alpha= 0.5)
    ax1.set_ylabel('Sen1∆')
    

    ax1.plot(bw.index.values, bw['smoo_mean'], color = 'orange', alpha= 0.5)
    ax1.plot(bc.index.values, bc['smoo_mean'], color = 'orange', alpha= 0.5)
    ax1.set_ylabel('Dbl8∆')
        

    ax1.plot(dw.index.values, dw['smoo_mean'], color = 'tomato', alpha= 0.5)

    ax1.plot(dc.index.values, dc['smoo_mean'], color = 'tomato', alpha= 0.5)
    ax1.set_ylabel('Sen1∆Dbl8∆')

    ax1.set_xticks([0,100,1100])
    ax1.set_xticklabels(['+0.1Kb','TSS','-1kb'])
    

sen1forstalls = plotall(wtforflanksenW, wtforflanksenC, senforflanksenW, senforflanksenC, dblforflanksenW, dblforflanksenC, dsforflanksenW, dsforflanksenC, 'sen1d  stall genes tss flank')
dbl8forstalls = plotall(wtforflankdblW, wtforflankdblC, senforflankdblW, senforflankdblC, dblforflankdblW, dblforflankdblC, dsforflankdblW, dsforflankdblC, 'dbld  stall genes TSS flank')
dsstalls = plotall(wtforflankdsW,wtforflankdsC,senforflankdsW, senforflankdsC, dblforflankdsW, dblforflankdsC, dsforflankdsW, dsforflankdsC, 'double mutant stall genes TSS flank' )

control = plotall(wtforflankconW,wtforflankconC,senforflankconW, senforflankconC, dblforflankconW, dblforflankconC, dsforflankconW, dsforflankconC, 'control genes TSS flank' )

sel10 = plotall(wtforflanksel10W, wtforflanksel10C, senforflanksel10W, senforflanksel10C, dblforflanksel10W, dblforflanksel10C, dsforflanksel10W, dsforflanksel10C, 'random selection of 9 genes tss flank')
sel200 = plotall(wtforflanksel200W, wtforflanksel200C, senforflanksel200W, senforflanksel200C, dblforflanksel200W, dblforflanksel200C, dsforflanksel200W, dsforflanksel200C, 'random selection of 195 genes tss flank')







def plotall(ww, sw, bw, dw, title, intwt,intsen, intdbl, intds):
    fi,((ax1)) = plt.subplots(1,1, sharex=True, sharey=True) 

    ax1.plot(ww.index.values, ww['smoo_mean'], color = 'black', alpha= 0.5)
 #   ax1.plot(wc.index.values, wc['smoo_mean'], color = 'black', alpha= 0.5)
    ax1.set_title(title)
    ax1.set_ylabel('WT')
    x = np.arange(0, 6001, 1)
    y1 = np.array(intwt['lower']).flatten()
    y2= np.array(intwt['upper']).flatten()
    ax1.fill_between(x, y1, y2, facecolor='grey', alpha =0.2)


    ax1.plot(sw.index.values, sw['smoo_mean'], color = 'steelblue', alpha= 0.5)
#    ax1.plot(sc.index.values, sc['smoo_mean'], color = 'steelblue', alpha= 0.5)
    ax1.set_ylabel('Sen1∆')
    y1s =np.array(intsen['lower']).flatten()
    y2s=np.array(intsen['upper']).flatten()
    ax1.fill_between(x, y1s, y2s, facecolor='deepskyblue', alpha =0.2, )
    

    ax1.plot(bw.index.values, bw['smoo_mean'], color = 'orange', alpha= 0.5)
 #   ax1.plot(bc.index.values, bc['smoo_mean'], color = 'orange', alpha= 0.5)
    ax1.set_ylabel('Dbl8∆')
    y1b =np.array(intdbl['lower']).flatten()
    y2b=np.array(intdbl['upper']).flatten()
    ax1.fill_between(x, y1b, y2b, facecolor='gold', alpha =0.2)
        

    ax1.plot(dw.index.values, dw['smoo_mean'], color = 'tomato', alpha= 0.5)
    y1ds =np.array(intds['lower']).flatten()
    y2ds=np.array(intds['upper']).flatten()
    ax1.fill_between(x, y1ds, y2ds, facecolor='indianred', alpha =0.1)

 #   ax1.plot(dc.index.values, dc['smoo_mean'], color = 'tomato', alpha= 0.5)
    ax1.set_ylabel('Sen1∆Dbl8∆')

    ax1.set_xticks([0,3000,6000])
    ax1.set_xticklabels(['3Kb','TSS','-3kb'])
    

sen1forstalls = plotall(wtforflanksenW, senforflanksenW, dblforflanksenW, dsforflanksenW, 'sen1d  stall genes tss flank',senWw,senWs,senWb, senWd)
dbl8forstalls = plotall(wtforflankdblW, senforflankdblW, dblforflankdblW, dsforflankdblW, 'dbld  stall genes TSS flank',dblWw, dblWs, dblWb, dblWd)
dsstalls = plotall(wtforflankdsW,senforflankdsW, dblforflankdsW, dsforflankdsW, 'double mutant stall genes TSS flank', dsWw, dsWs,dsWb, dsWd)

control = plotall(wtforflankconW,senforflankconW, dblforflankconW, dsforflankconW, 'control genes TSS flank',conW,conWs, conWb, conWd )


sel10 = plotall(wtforflanksel10W, senforflanksel10W, dblforflanksel10W, dsforflanksel10W, 'random selection of 9 genes tss flank')




figgy, (ax1) =plt.subplots(1, sharey=True)
ax1.plot(awline.index, awline['SmoothedData'], color = 'black', alpha =0.8, linewidth=1)
ax1.plot(assline.index, assline['SmoothedData'], color = 'deepskyblue', alpha=0.8, linewidth=1)
ax1.plot(adline.index, adline['SmoothedData'], color = 'orange', alpha=0.8, linewidth=1)
ax1.plot(adsline.index, adsline['SmoothedData'], color = 'tomato', alpha=0.8, linewidth=1)
x = np.arange(0, 1000, 1)
y1 = np.array(saCON['lower']).flatten()
y2= np.array(saCON['upper']).flatten()
ax1.fill_between(x, y1, y2, facecolor='grey', alpha =0.2)
y1s =np.array(sasCON['lower']).flatten()
y2s=np.array(sasCON['upper']).flatten()
ax1.fill_between(x, y1s, y2s, facecolor='deepskyblue', alpha =0.2, )
y1b =np.array(sabCON['lower']).flatten()
y2b=np.array(sabCON['upper']).flatten()
#ax1.fill_between(x, y1b, y2b, facecolor='gold', alpha =0.2)
y1ds =np.array(sadsCON['lower']).flatten()
y2ds=np.array(sadsCON['upper']).flatten()
ax1.fill_between(x, y1ds, y2ds, facecolor='indianred', alpha =0.1)
ax1.set_xticks([0,1000])
ax1.set_xticklabels(['TSS','TES'])
#ax1.set_ylim([-4, 124])
ax1.set_ylabel('normalised reads')
ax1.legend(loc='best')

ax1.set_xticks([0,1000])
ax1.set_xticklabels(['TSS','TES'])
ax1.legend(loc='best')
ax1.set_ylabel('normalised reads')


wtforflanksel200W = take_mean(wtforflanksel200W)
wtforflanksel200C = take_mean(wtforflanksel200C)
senforflanksel200W = take_mean(senforflanksel200W)
senforflanksel200C = take_mean(senforflanksel200C)
dblforflanksel200W = take_mean(dblforflanksel200W)
dblforflanksel200C = take_mean(dblforflanksel200C)
dsforflanksel200W = take_mean(dsforflanksel200W)
dsforflanksel200C = take_mean(dsforflanksel200C)


rolling_avg = df.rolling(window=window_size, axis=1).mean()

# Print the resulting DataFrame
print(rolling_avg)


#%%
#writ  asubst function for centromeres 
aaaah = {'start':[3753687,1602264,1070904], 'end': [3789421,1644747,1137003], 'Chromosome': [1,2,3]}
comcentro
ah1 = {'start'}


def centrocheck(file,frame1):
    data_output = pd.DataFrame()

    for i in range(len(file)):
        tempstart = file.iloc[i]["start"]
        tempend = file.iloc[i]["end"]
      #  tempchro = file.iloc[i]["Chromosome"]
        
        length = tempend - tempstart


        tempsubsetfor = frame1.loc[(frame1['Pos']>= tempstart) & (frame1['Pos'] <= tempend)]

        tempsubsetfor = tempsubsetfor.reset_index(drop = True)
        
        line = pd.DataFrame({"Pos": tempend}, index=[length])
        line2 = pd.DataFrame({"Pos": tempstart}, index=[0])
        
        
        tempsubsetfor = tempsubsetfor.append(line, ignore_index=False)
        tempsubsetfor = tempsubsetfor.append(line2, ignore_index=False)
         
        tempsubsetfor = tempsubsetfor.sort_values('Pos', ascending=True)


        for x in range(len(tempsubsetfor)):
            pos = tempsubsetfor.iloc[x]['Pos']
               # print(pos)
            tpos = tempsubsetfor.iloc[0]['Pos']

            tempsubsetfor.loc[tempsubsetfor.index[x], 'aPos'] = (pos - tpos)

        tempsubsetfor = tempsubsetfor.set_index('aPos')
            
            

        data_blank = pd.DataFrame({
    
            "aPos" : np.arange(0.0, length, 1)}).merge(
        
                tempsubsetfor, on = "aPos", how = "left"

                )

        data_blank = data_blank[["aPos", "normscore"]]
            
        data_output = pd.concat([data_output, data_blank['normscore']], axis =1, ignore_index=True)
          
            
        data_output = data_output.fillna(0)



    return data_output
centrowt1 = centrocheck(centro1,wtchr1)
centrowt2 = centrocheck(centro2,wtchr2)
centrowt3 = centrocheck(centro3,wtchr3)

centot2t1 = centrocheck(centro1,t2tchr1)
centot2t2 = centrocheck(centro2,t2tchr2)
centot2t3 = centrocheck(centro3,t2tchr3)

centot2m1 = centrocheck(centro1,t2mchr1)
centot2m2 = centrocheck(centro2,t2mchr2)
centot2m3 = centrocheck(centro3,t2mchr3)

centot2r1 = centrocheck(centro1,t2rchr1)
centot2r2 = centrocheck(centro2,t2rchr2)
centot2r3 = centrocheck(centro3,t2rchr3)

centot2e1 = centrocheck(centro1,t2echr1)
centot2e2 = centrocheck(centro2,t2echr2)
centot2e3 = centrocheck(centro3,t2echr3)

#t2tcen = pd.concat([centot2t1, centot2t2, centot2t3], axis=1)

# Rename the columns if needed
merged_df.columns = ['column1', 'column2', 'column3']


def cen_plot (wt1, wt2, wt3, t1, t2, t3, m1, m2, m3):
    ff, (ax2, ax3, ax4) = plt.subplots(3,1, sharex=True)
    #ax1.set_title('WT')
    ax2.plot(cc1['Pos'], cc1['Wnorm'], color ='firebrick', alpha=0.8)
    ax2.plot(cc1['Pos'], cc1['Cnorm'], color ='steelblue', alpha=0.8)
   # ax2.set_ylim(-0.5,0.5)
    ax2.set_ylabel('WT (HpM)')



#now i need to check by a plot
def plotall(wt1, wt2, wt3, t1, t2, t3, m1, m2, m3):
    fi,((ax1, ax2, ax3)) = plt.subplots(3,1, sharex=True)  
    for column in wt1:
        ax1.plot(wt1.index.values, wt1[column], color = 'dodgerblue', alpha= 0.5)
    for column in t1:
        ax1.plot(t1.index.values, t1[column], color = 'goldenrod', alpha= 0.5)
    for column in m1:
        ax1.plot(m1.index.values, m1[column], color = 'maroon', alpha= 0.3)
      #  ax1.set_title(title)
      #  ax1.set_title('Chromosome 1 centromere')
        ax1.set_ylim(-0.04, 0.4)
    
    for column in wt2:
        ax2.plot(wt2.index.values, wt2[column], color = 'dodgerblue', alpha= 0.5)
    for column in t2:
        ax2.plot(t2.index.values, t2[column], color = 'goldenrod', alpha= 0.5)
    for column in m2:
        ax2.plot(m2.index.values, m2[column], color = 'maroon', alpha= 0.3)
      #  ax1.set_title(title)
      #  ax2.set_title('Chromosome 2 centromere')
        ax2.set_ylim(-0.04, 0.6)

    for column in wt3:
        ax3.plot(wt3.index.values, wt3[column], color = 'dodgerblue', alpha= 0.5)
    for column in t3:
        ax3.plot(t3.index.values, t3[column], color = 'goldenrod', alpha= 0.5)
    for column in m3:
        ax3.plot(m3.index.values, m3[column], color = 'maroon', alpha= 0.3)
      #  ax1.set_title(title)
      #  ax3.set_title('Chromosome 3 centromere')
        ax3.set_ylim(-0.04, 0.6)


        
    ax1.axes.xaxis.set_ticks([])
    ax2.axes.xaxis.set_ticks([])
    
    ax3.axes.xaxis.set_ticks([])
   # ax4.set_xticklabels(['-TSS','-1kb'])
    return fi

sen1forstalls = plotall(centrowt1, centrowt2, centrowt3, centot2t1, centot2t2, centot2t3, centot2m1, centot2m2, centot2m3)

#%%




# Assuming your DataFrames are named df1, df2, and df3
dataframes = [centrowt1, centrowt2, centrowt3]
dataframest = [centot2t1, centot2t2, centot2t3]

dataframesr = [centot2r1, centot2r2, centot2r3]
dataframeste = [centot2e1, centot2e2, centot2e3]

# Calculate the middle index for each DataFrame
middle_indices = [df.shape[0] // 2 for df in dataframes]
middle_indicest = [dft.shape[0] // 2 for dft in dataframest]
middle_indicesr = [dfr.shape[0] // 2 for dfr in dataframesr]
middle_indiceste = [dfte.shape[0] // 2 for dfte in dataframeste]

# Calculate the maximum length of the DataFrames
max_length = max([df.shape[0] for df in dataframes])
max_lengtht = max([dft.shape[0] for dft in dataframest])
max_lengthr = max([dfr.shape[0] for dfr in dataframesr])
max_lengthte = max([dfte.shape[0] for dfte in dataframeste])

# Define colors and alphas for the lines
colors = ['royalblue', 'forestgreen', 'teal']
colors2 = ['firebrick', 'orangered', 'orange']

alphas = [0.8, 0.8, 0.7]




####FOR WT

middle_indices = [df.shape[0] // 2 for df in [centrowt1, centrowt2, centrowt3]]

# Calculate the maximum length of the DataFrames
max_length = max([df.shape[0] for df in [centrowt1, centrowt2, centrowt3]])


# Create new DataFrames with centered indices
centered_dfs = []
for df, middle_index in zip([centrowt1, centrowt2, centrowt3], middle_indices):
    offset = max_length // 2 - middle_index
    centered_df = df.copy()
    centered_df.index = range(offset, offset + df.shape[0])
    centered_dfs.append(centered_df)

# Merge the centered DataFrames on index
wtcentroo = pd.concat(centered_dfs, axis=1)

# Rename the columns if needed
wtcentroo.columns = ['column1', 'column2', 'column3']

mean_values = wtcentroo.mean(axis=1)

# Create a new DataFrame with the mean values
mean_wtcentroo = pd.DataFrame(mean_values, columns=['Mean'])

# Display the mean DataFrame
print(mean_wtcentroo)


window_size = 500

# Calculate the rolling window average
rolling_average_wt = mean_wtcentroo['Mean'].rolling(window=window_size, center=True).mean()

# Create a new DataFrame with the rolling average values
rolling_average_df_WT = pd.DataFrame({'Rolling Average': rolling_average_wt})




####FOR WT

middle_indices = [df.shape[0] // 2 for df in [centot2t1, centot2t2, centot2t3]]

# Calculate the maximum length of the DataFrames
max_length = max([df.shape[0] for df in [centot2t1, centot2t2, centot2t3]])


# Create new DataFrames with centered indices
centered_dfs = []
for df, middle_index in zip([centot2t1, centot2t2, centot2t3], middle_indices):
    offset = max_length // 2 - middle_index
    centered_df = df.copy()
    centered_df.index = range(offset, offset + df.shape[0])
    centered_dfs.append(centered_df)

# Merge the centered DataFrames on index
t2tcentro = pd.concat(centered_dfs, axis=1)

# Rename the columns if needed
t2tcentro.columns = ['column1', 'column2', 'column3']

mean_values = t2tcentro.mean(axis=1)

# Create a new DataFrame with the mean values
mean_t2tcentro = pd.DataFrame(mean_values, columns=['Mean'])

# Display the mean DataFrame
print(mean_t2tcentro)

#window_size = 3

# Calculate the rolling window average
rolling_average_t2t = mean_t2tcentro['Mean'].rolling(window=window_size, center=True).mean()

# Create a new DataFrame with the rolling average values
rolling_average_df_t2t = pd.DataFrame({'Rolling Average': rolling_average_t2t})


####for t2 repeatt
# Create new DataFrames with centered indices

middle_indicesr = [df.shape[0] // 2 for df in [centot2r1, centot2r2, centot2r3]]

# Calculate the maximum length of the DataFrames
max_lengthr = max([df.shape[0] for df in [centot2r1, centot2r2, centot2r3]])

centered_dfsr = []
for dfr, middle_index in zip([centot2r1, centot2r2, centot2r3], middle_indicesr):
    offset = max_lengthr // 2 - middle_index
    centered_df = dfr.copy()
    centered_df.index = range(offset, offset + dfr.shape[0])
    centered_dfsr.append(centered_df)

# Merge the centered DataFrames on index
t2rcentro = pd.concat(centered_dfsr, axis=1)

# Rename the columns if needed
t2rcentro.columns = ['column1', 'column2', 'column3']

mean_values = t2rcentro.mean(axis=1)

# Create a new DataFrame with the mean values
mean_t2rcentro = pd.DataFrame(mean_values, columns=['Mean'])

# Display the mean DataFrame
print(mean_t2rcentro)

#window_size = 3

# Calculate the rolling window average
rolling_average_t2r = mean_t2rcentro['Mean'].rolling(window=window_size, center=True).mean()

# Create a new DataFrame with the rolling average values
rolling_average_df_t2r = pd.DataFrame({'Rolling Average': rolling_average_t2r})


#####for ethanopl


middle_indiceste = [df.shape[0] // 2 for df in [centot2e1, centot2e2, centot2e3]]
max_lengthte = max([df.shape[0] for df in [centot2e1, centot2e2, centot2e3]])

# Calculate the maximum length of the DataFrames
max_lengthte = max([df.shape[0] for df in [centot2e1, centot2e2, centot2e3]])
centered_dfse = []
for dfe, middle_index in zip([centot2e1, centot2e2, centot2e3], middle_indiceste):
    offset = max_lengthte // 2 - middle_index
    centered_df = dfe.copy()
    centered_df.index = range(offset, offset + dfe.shape[0])
    centered_dfse.append(centered_df)

# Merge the centered DataFrames on index
t2ecentro = pd.concat(centered_dfse, axis=1)

# Rename the columns if needed
t2ecentro.columns = ['column1', 'column2', 'column3']

mean_values = t2ecentro.mean(axis=1)

# Create a new DataFrame with the mean values
mean_t2ecentro = pd.DataFrame(mean_values, columns=['Mean'])

# Display the mean DataFrame
print(mean_t2ecentro)

#window_size = 3

# Calculate the rolling window average
rolling_average_t2e = mean_t2ecentro['Mean'].rolling(window=window_size, center=True).mean()

# Create a new DataFrame with the rolling average values
rolling_average_df_t2e = pd.DataFrame({'Rolling Average': rolling_average_t2e})




# Create a figure and axis
fig, ((ax1, ax2)) = plt.subplots(2,1, sharex=True)

# Plot each DataFrame centered on the shared axis
for i, df in enumerate(dataframes):
    offset = max_length // 2 - middle_indices[i]
    x = range(offset, offset + df.shape[0])
    y = df.iloc[:, 0]
    ax1.plot(x, y, color=colors[i], alpha=alphas[i])
    

for b, dft in enumerate(dataframest):
    offsett = max_lengtht // 2 - middle_indicest[b]
    xt = range(offsett, offsett + dft.shape[0])
    yt = dft.iloc[:, 0]
    ax2.plot(xt, yt, color=colors2[b], alpha=alphas[b])
   # ax2.plot(mean_t2tcentro.index.values, mean_t2tcentro['Mean'], color='black')
    

# Set labels and title
#ax1.set_xlabel('centromeric signal')

ax2.set_xlabel('position relative to centromere midpoint (bp)')
#ax2.set_xlabel('centromeric signal')
#ax2.set_ylabel('Aggregated CC-seq Signal (HpM)')
ax1.set_xticks([])
ax2.set_xticks([0,16525,33050,49575,66100])
ax2.set_xticklabels(['-33050','-16525','0','+16525','+33050'])
ax1.set_title('Wild Type')
ax2.set_title('Top2-191t')

# Show the plot
plt.show()



fig, ((ax1)) = plt.subplots(1,1, sharex=True)
ax1.plot(rolling_average_df_WT.index.values, rolling_average_df_WT['Rolling Average'], color='black')
ax1.plot(rolling_average_df_t2t.index.values, rolling_average_df_t2t['Rolling Average'], color='tomato')
ax1.plot(rolling_average_df_t2r.index.values, rolling_average_df_t2r['Rolling Average'], color='crimson')
ax1.plot(rolling_average_df_t2e.index.values, rolling_average_df_t2e['Rolling Average'], color='deepskyblue')

#ax1.set_xlabel('midpoint alligned centromeric signal (relative bp)')
#ax1.set_ylabel('Aggregated CC-seq Signal (100bp smoothed)(HpM)')
ax1.set_xlabel('midpoint alligned centromeric signal (relative bp)')
#ax2.set_ylabel('Aggregated CC-seq Signal (100bp smoothed)(HpM)')
ax1.set_xticks([0,16525,33050,49575,66100])
ax1.set_xticklabels(['-33050','-16525','0','+16525','+33050'])
#ax1.set_title('Wild Type')
#ax2.set_title('Top2-191t')

#%%

import scipy.stats as stats
from scipy.stats import t

def conint(array):
    confidence = 0.95
    trythis = []
    x = pd.DataFrame(array)
    #print(x)
    for column in x:

        m = (x[column]).mean()
        s = x[column].std()

        dof = len(x[column])-1 
        t_crit = np.abs(t.ppf((1-confidence)/2,dof))
        interval = (m-s*t_crit/np.sqrt(len(x)), m+s*t_crit/np.sqrt(len(x)))
        
        trythis.append(interval)
        saint = pd.DataFrame(trythis, columns=['lower', 'upper'])
        #oft = saint.T
    return saint

saCON = conint(wtforflanksenW)
sasCON = conint(senforflanksenW)
sabCON = conint(dblforflanksenW)
sadsCON = conint(dsforflanksenW)

saCONc = conint(wtforflanksenC)
sasCONc = conint(senforflanksenC)
sabCONc = conint(dblforflanksenC)
sadsCONc = conint(dsforflanksenC)


def plotall(ww, wc, sw, sc, bw, bc, dw,dc, title):
    fi,((ax1, ax2, ax3,ax4)) = plt.subplots(4,1, sharex=True, sharey=True)  
    for column in ww:
        ax1.plot(ww.index.values, ww[column], color = 'black', alpha= 0.3)
    for column in wc:
        ax1.plot(wc.index.values, wc[column], color = 'black', alpha= 0.3)
        ax1.set_title(title)
        ax1.set_ylabel('WT')

    for column in sw:
        ax2.plot(sw.index.values, sw[column], color = 'steelblue', alpha= 0.3)
    for column in sc:
        ax2.plot(sc.index.values, sc[column], color = 'steelblue', alpha= 0.3)
        ax2.set_ylabel('Sen1∆')
    
    for column in bw:
        ax3.plot(bw.index.values, bw[column], color = 'orange', alpha= 0.3)
    for column in bc:
        ax3.plot(bc.index.values, bc[column], color = 'orange', alpha= 0.3)
        ax3.set_ylabel('Dbl8∆')
        
    for column in dw:
        ax4.plot(dw.index.values, dw[column], color = 'tomato', alpha= 0.3)
    for column in dc:
        ax4.plot(dc.index.values, dc[column], color = 'tomato', alpha= 0.3)
        ax4.set_ylabel('Sen1∆Dbl8∆')

    ax1.axes.xaxis.set_ticks([])
    ax2.axes.xaxis.set_ticks([])
    ax3.axes.xaxis.set_ticks([])
    ax4.set_xticks([0,1000])
    ax4.set_xticklabels(['-TSS','-1kb'])
    return fi

sen1forstalls = plotall(wtforflanksenW, wtforflanksenC, senforflanksenW, senforflanksenC, dblforflanksenW, dblforflanksenC, dsforflanksenW, dsforflanksenC, 'sen1d forward strand stall genes tss flank')
dbl8forstalls = plotall(wtforflankdblW, wtforflankdblC, senforflankdblW, senforflankdblC, dblforflankdblW, dblforflankdblC, dsforflankdblW, dsforflankdblC, 'dbld forward strand stall genes TSS flank')
control = plotall(wtforflankconW,wtforflankconC,senforflankconW, senforflankconC, dblforflankconW, dblforflankconC, dsforflankconW, dsforflankconC, 'control genes TSS flank' )

import seaborn as sns
import matplotlib.pyplot as plt

def plot_heatmap(data, vmin=None, vmax=None):
    sns.heatmap(data, cmap='coolwarm', vmin=vmin, vmax=vmax)
    plt.xlabel('Columns')
    plt.ylabel('Rows')
    plt.title('Heatmap')
    plt.show()
    
v = plot_heatmap(wtforflankdblW, vmin= 0, vmax=0.1)

#sen1 stall sties flank minus 


#dbl8 stall sties flank minus 
fi,((ax1, ax2, ax3,ax4)) = plt.subplots(4,1, sharex=True, sharey=True)  
ax1.plot(wtforflankdblW.index.values, wtforflankdblW['Wmean'], color = 'blue', alpha= 0.8)
ax1.plot(wtforflankdblC.index.values, wtforflankdblC["Cmean"], color = 'tomato', alpha= 0.8)
ax1.set_title('Dbl8∆ stall sites')
ax1.set_ylabel('WT')
ax2.plot(senforflankdblW.index.values, senforflankdblW['Wmean'], color = 'blue', alpha= 0.8)
ax2.plot(senforflankdblC.index.values, senforflankdblC["Cmean"], color = 'tomato', alpha= 0.8)
ax2.set_ylabel('Sen1∆')
ax3.plot(dblforflankdblW.index.values, dblforflankdblW['Wmean'], color = 'blue', alpha= 0.8)
ax3.plot(dblforflankdblC.index.values, dblforflankdblC["Cmean"], color = 'tomato', alpha= 0.8)
ax3.set_ylabel('Dbl8∆')
ax4.plot(dsforflankdblW.index.values, dsforflankdblW['Wmean'], color = 'blue', alpha= 0.8)
ax4.plot(dsforflankdblC.index.values, dsforflankdblC["Cmean"], color = 'tomato', alpha= 0.8)
ax4.set_ylabel('Sen1∆Dbl8∆')

ax1.axes.xaxis.set_ticks([])
ax2.axes.xaxis.set_ticks([])
ax3.axes.xaxis.set_ticks([])
ax4.set_xticks([0,1000])
ax4.set_xticklabels(['-1kb','TSS'])



#%%


def Chromosome_plot (data1, data2, data3, centro, genee, telo):
    ff, (ax1, ax2, ax3, ax4) = plt.subplots(4,1, sharex=True)
    #ax1.set_title('WT')
    ax1.plot(data1['Pos'], data1['Wnorm'], color ='firebrick', alpha=0.8)
    ax1.plot(data1['Pos'], data1['Cnorm'], color ='steelblue', alpha=0.8)
    ax1.set_ylim(-10,20)
  #  ax1.set_xlim(3753687,3789421)
    ax1.set_ylabel('WT (HpM)')
    
    #ax2.set_title('top2-191t')
    ax2.plot(data2['Pos'], data2['Wnorm'], color ='coral', alpha=0.8)
    ax2.plot(data2['Pos'], data2['Cnorm'], color ='teal', alpha=0.8)
    ax2.set_ylim(-5,5)
    ax2.set_ylabel('top2-191t (HpM)')
    
    #ax3.set_title('top2-191m')
    ax3.plot(data3['Pos'], data3['Wnorm'], color ='coral', alpha=0.8)
    ax3.plot(data3['Pos'], data3['Cnorm'], color ='teal', alpha=0.8)
    ax3.set_ylim(-5,5)
    ax3.set_ylabel('top2-191m (HpM)')
    ax4.set_xlabel('Chromosome position')

                  
  #  for fe in featurex.itertuples(index=False, name=None):
   #     if fe[6] == 'reverse':
    #        ax4.axvspan(fe[3],fe[4],0.3,0.5,color="teal",alpha=0.3)
            
     #   elif fe[6] == 'forward':
      #          ax4.axvspan(fe[3],fe[4],0.5,0.8,color="coral",alpha=0.3)
       #         ax4.set_ylabel('Gene annotations')
                
    
 #   for n in no.itertuples(index=False, name=None):
  #      if n[6] == 'reverse':
   #         ax4.axvspan(n[3],n[4],0.3,0.5,color="lightsteelblue", alpha=0.3)
            
    #    elif n[6] == 'forward':
     #           ax4.axvspan(n[3],n[4],0.5,0.8,color="mistyrose", alpha=0.3)
                
    
    for c in centro.itertuples(index=False, name=None):
            ax4.axvspan(c[0],c[1],0.1,0.15,color="black",alpha=0.8)
            
    for t in telo.itertuples(index=False, name=None):
            ax4.axvspan(t[0],t[1],0.1,0.15,color="grey",alpha=0.8)
            
    for ge in genee.itertuples(index=False, name=None):
        if ge[1] == 'protein coding':
            if ge[5] == 'reverse':
                ax4.axvspan(ge[2],ge[3],0.2,0.3,color="lightsteelblue",alpha=0.5)
            elif ge[5] == 'forward':
                    ax4.axvspan(ge[2],ge[3],0.8,0.9,color="mistyrose",alpha=0.5)
                    '''
    posit =[]
    cow = [] # too slow to plot each pont, so first make 2 lists
    for feat in Net.itertuples(index=False, name=None):
        if feat[1] == 'chr1' and feat[3] == 'forward':
            #log_2 = np.log2(feat[3])
            posit.append(feat[2])
            cow.append(feat[5])
        elif feat[1] == 'chr1' and feat[3] == 'reverse':
            #log_2 = np.log2(feat[3])
            posit.append(feat[2])
            cow.append(feat[5]*-1)

    ax5.plot(posit,cow)
    #ax5.set_ylim(-5, +5)
'''
    return ff

chromosome1 = Chromosome_plot(wtchr1, t2tchr1, t2mchr1, centro1, gene1, telo1)
chromosome2 = Chromosome_plot(wtchr2, t2tchr2, t2mchr2, centro2, gene2, telo2)
chromosome3 = Chromosome_plot(wtchr3, t2tchr3, t2mchr3, centro3, gene3, telo3)
