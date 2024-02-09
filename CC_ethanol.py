#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 11:56:39 2023

@author: joannafernandez
"""

""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



def Collect_data(file1,fact):
    full = pd.read_csv(file1, sep="	")
    #reverse = pd.read_csv(file2)
    full['score'] = full['Watson'] + full['Crick']
  #  print(full)
    print (full['score'].sum())
    full['Wnorm']= (full['Watson'] / full['score'].sum())*fact
    full['Cnorm']= (full['Crick'] / full['score'].sum())*fact
    full['Cnorm']= full['Cnorm']*-1
    full['normscore']= (full['score']/ full['score'].sum())*fact

   # print(full)
    chr1= full[(full["Chr"] == 1)]
    chr2= full[(full["Chr"] == 2)]
    chr3= full[(full["Chr"] == 3)]
    return chr1,chr2,chr3,full

#original
wtchr1, wtchr2, wtchr3, wt = Collect_data('FullMap.pombase220208-WG141_WT6mins.txt', 10000)
t2tchr1, t2tchr2, t2tchr3, t2t = Collect_data('FullMap.pombase220208-WG142_top2-191t.txt',10000)
t2mchr1, t2mchr2, t2mchr3, t2m = Collect_data('FullMap.pombase220208-WG143_top2-191m.txt',10000)

#old
wtochr1, wtochr2, wtochr3, wto = Collect_data('FullMap.pombase220208-WT-old.txt', 10000)
t2rchr1, t2rchr2, t2rchr3, t2r = Collect_data('FullMap.pombase220208-T2-old.txt',10000)

#ethanol
wtETchr1, wtETchr2, wtETchr3, wtET = Collect_data('FullMap.pombase220208-WT-Et.txt', 10000)
wtETBchr1, wtETBchr2, wtETBchr3, wtETB = Collect_data('FullMap.pombase220208-WT-Et-b.txt', 10000)
wtETCchr1, wtETCchr2, wtETCchr3, wtETC = Collect_data('FullMap.pombase220208-WT-Et-c.txt', 10000)

t2Echr1, t2Echr2, t2Echr3, t2E = Collect_data('FullMap.pombase220208-T2E.txt',10000)
t2EBchr1, t2EBchr2, t2EBchr3, t2EB = Collect_data('FullMap.pombase220208-T2E-b.txt',10000)

formchr1, formchr2, formchr3, form = Collect_data('FullMap.pombase220208-WT-FORM.txt',10000)


#

#averaged
wtETchr1, wtETchr2, wtETchr3, wtET = Collect_data('FullMap.pombase220208-WT.txt', 10000)
t2Echr1, t2Echr2, t2Echr3, t2E = Collect_data('FullMap.pombase220208-T2.txt',10000)

#offset filtered
wtETchr1, wtETchr2, wtETchr3, wtET = Collect_data('FullMap.pombase220208-WT_FinalAverage_3bpOffset_with_strandminTH_0_p0.txt', 10000)
t2Echr1, t2Echr2, t2Echr3, t2E = Collect_data('FullMap.pombase220208-T2_FinalAverage_3bpOffset_with_strandminTH_0_p0.txt',10000)



#####
#averaged final
wtETchr1, wtETchr2, wtETchr3, wtET = Collect_data('FullMap.7.pombase220208-WT-Et_8.pombase220208-WT-Et-b_Average.txt', 10000)
t2Echr1, t2Echr2, t2Echr3, t2E = Collect_data('FullMap.3.pombase220208-T2E_4.pombase220208-T2E-b_Average.txt',10000)

wtET35chr1, wtET35chr2, wtET35chr3, wtET35 = Collect_data('FullMap.5.pombase220208-wt35_S11_6.pombase220208-wt35_B_S13_Average.txt', 10000)
t2E30chr1, t2E30chr2, t2E30chr3, t2E30 = Collect_data('FullMap.1.pombase220208-t2_p_S14_2.pombase220208-t2_p2_S16_Average.txt',10000)



#averaged BINNED
bwtETchr1, bwtETchr2, bwtETchr3, bwtET = Collect_data('FullMap.7.pombase220208-WT-Et_8.pombase220208-WT-Et-b_Average_0.1Kbp_binned_OS3_full_strand.txt', 10000)
bt2Echr1, bt2Echr2, bt2Echr3, bt2E = Collect_data('FullMap.3.pombase220208-T2E_4.pombase220208-T2E-b_Average_0.1Kbp_binned_OS3_full_strand.txt',10000)

bwtET35chr1, bwtET35chr2, bwtET35chr3, bwtET35 = Collect_data('FullMap.5.pombase220208-wt35_S11_6.pombase220208-wt35_B_S13_Average_0.1Kbp_binned_OS3_full_strand.txt', 10000)
bt2E30chr1, bt2E30chr2, bt2E30chr3, bt2E30 = Collect_data('FullMap.1.pombase220208-t2_p_S14_2.pombase220208-t2_p2_S16_Average_0.1Kbp_binned_OS3_full_strand.txt',10000)



#binned individual for controls:
bt2pchr1, bt2pchr2, bt2pchr3, bt2p = Collect_data('FullMap.pombase220208-t2_p_S14_0.1Kbp_binned_OS0_full_strand.txt', 10000)
bt2p2chr1, bt2p2chr2, bt2p2chr3, bt2p2 = Collect_data('FullMap.pombase220208-t2_p2_S16_0.1Kbp_binned_OS0_full_strand.txt', 10000)

bwt35chr1, bwt35chr2, bwt35chr3, bwt35 = Collect_data('FullMap.pombase220208-wt35_S11_0.1Kbp_binned_OS0_full_strand.txt', 10000)
bwt35bchr1,bwt35bchr2,bwt35bchr3, bwt35b = Collect_data('FullMap.pombase220208-wt35_B_S13_0.1Kbp_binned_OS0_full_strand.txt', 10000)

    

#binned
wtETchr1, wtETchr2, wtETchr3, wtET = Collect_data('FullMap.pombase220208-WT-Et_0.1Kbp_binned_OS3_full_strand.txt', 10000)
wtETBchr1, wtETBchr2, wtETBchr3, wtETB = Collect_data('FullMap.pombase220208-WT-Et-b_0.1Kbp_binned_OS3_full_strand.txt', 10000)


t2Echr1, t2Echr2, t2Echr3, t2E = Collect_data('FullMap.pombase220208-T2E_0.1Kbp_binned_OS3_full_strand.txt',10000)
t2EBchr1, t2EBchr2, t2EBchr3, t2EB = Collect_data('FullMap.pombase220208-T2E-b_0.1Kbp_binned_OS3_full_strand.txt',10000)
#%%
bwtET['source'] = 'WT'
bt2E['source'] = 'top2-191'
bwtET35['source'] = 'WT-35.5'
bt2E30['source'] = 'top2-191-30'

# Concatenate the dataframes
combined_averaged = pd.concat([wtET, t2E, wtET35, t2E30], ignore_index=True)
combined_averaged['log2_score'] = np.log2(combined_averaged['normscore'])
import seaborn as sns
plt.figure(figsize=(10, 8))
sns.violinplot(data=combined_averaged, x="source", y="log2_score", hue='source',inner="quart")



# Print the result dataframe
print(result_df)



from scipy.stats import wilcoxon

# Assuming df1, df2, df3, df4 are your dataframes

# Create a list of dataframes
dataframes = [bwtET, bt2E, bwtET35, bt2E30]

wilcoxon_results = []

# Iterate through pairs of dataframes and conduct Wilcoxon signed-rank test
for i in range(len(dataframes)-1):
    for j in range(i+1, len(dataframes)):
        # Extract scores from the dataframes
        scores1 = dataframes[i]['normscore']
        scores2 = dataframes[j]['normscore']
        
        # Conduct Wilcoxon signed-rank test
        stat, p_value = wilcoxon(scores1, scores2)
        
        # Store the results
        result = {
            'series 1': f'{dataframes[i]["source"].iloc[0]}',
            'series 2': f'{dataframes[j]["source"].iloc[0]}',
            'Statistic': stat,
            'P-value': p_value
        }
        wilcoxon_results.append(result)

# Create a dataframe from the results
wilcoxon_df = pd.DataFrame(wilcoxon_results)

# Print the results dataframe
print(wilcoxon_df)



pivoted_df = wilcoxon_df.pivot(index='series 1', columns='series 2', values='P-value')



# Create a heatmap using seaborn
plt.figure(figsize=(10, 4))
sns.heatmap(pivoted_df, annot=True, cmap='RdBu', fmt=".1g")
plt.title("Correlogram for p-values")
plt.show()









#%%

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

feat1 = ffeat[(ffeat['chro'] == 'chr1')]
feat2 = ffeat[(ffeat['chro'] == 'chr2')]
feat3 = ffeat[(ffeat['chro'] == 'chr3')]



aaaah = {'start':[3753687,1602264,1070904], 'end': [3789421,1644747,1137003], 'Chromosome': [1,2,3]}
comcentro = pd.DataFrame(data=aaaah)
comcentro['Sbinpos'] = comcentro['start']/50
comcentro['Sbinpos'] = comcentro['Sbinpos'].astype(int)
comcentro['Sbinpos'] = comcentro['Sbinpos']*50 +25
comcentro['Ebinpos'] = comcentro['end']/50
comcentro['Ebinpos'] = comcentro['Ebinpos'].astype(int)
comcentro['Ebinpos'] = comcentro['Ebinpos']*50 +25

#for trna id, remove 20kb either side


d = {'start': [3733687], 'end': [3809421], 'chro': ['chr1']}
centro1 = pd.DataFrame(data=d)
dd = {'start': [1582264], 'end': [1664747], 'chro': ['chr2']}
centro2 = pd.DataFrame(data=dd)

ddd = {'start': [1050904], 'end': [1157003], 'chro': ['chr3']}
centro3 = pd.DataFrame(data=ddd)


d = {'start': [3753687], 'end': [3789421], 'chro': ['chr1']}
centro1 = pd.DataFrame(data=d)

dd = {'start': [1602264], 'end': [1644747], 'chro': ['chr2']}
centro2 = pd.DataFrame(data=dd)

ddd = {'start': [1070904], 'end': [1137003], 'chro': ['chr3']}
centro3 = pd.DataFrame(data=ddd)

centro_merged = pd.concat([centro1, centro2, centro3], ignore_index=True)



te1 = {'start': [1, 5554844], 'end': [29663,5579133], 'chro': ['chr1', 'chr1']}
telo1 = pd.DataFrame(data=te1)
te2 ={'start': [1, 4500619], 'end': [39186,4539800], 'chro': ['chr2', 'chr2']}
telo2 = pd.DataFrame(data=te2)
te3 ={'start': [], 'end': []}
telo3 = pd.DataFrame(data=te3)

r3 = {'start': [1], 'end': [28297], 'chro': ['chr3'], 'chro': ['chr3']}
tr33 = pd.DataFrame(data=r3)

ends_merged = pd.concat([telo1, telo2, tr33], ignore_index=True)


def Chromosome_plot (cc1, cc2, cc3, cc4, centro, genee, telo):
    ff, (ax2, ax3, ax4, ax5, ax6) = plt.subplots(5,1, sharey=True)
    #ax1.set_title('WT')
    ax2.plot(cc1['Pos'], cc1['Wnorm'], color ='firebrick', alpha=0.8, linewidth=1.0)
    ax2.plot(cc1['Pos'], cc1['Cnorm'], color ='steelblue', alpha=0.8, linewidth=1.0)
    ax2.set_ylim(-2,2)
    ax2.set_ylabel('WT (HpM)')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    #ax2.set_ylim(-5,5)
   # ax2.set_ylabel('WT')

    ax3.plot(cc2['Pos'], cc2['Wnorm'], color ='firebrick', alpha=0.8, linewidth=1.0)
    ax3.plot(cc2['Pos'], cc2['Cnorm'], color ='steelblue', alpha=0.8, linewidth=1.0)
  #  ax3.set_ylim(-1,1)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.spines['bottom'].set_visible(False)
   # ax3.set_ylabel('SEN1d (HpM)')
    ax3.set_ylabel('WT-35 (HpM)')

    ax6.set_xlabel('Chromosome position')
    ax4.plot(cc3['Pos'], cc3['Wnorm'], color ='firebrick', alpha=0.8, linewidth=1.0)
    ax4.plot(cc3['Pos'], cc3['Cnorm'], color ='steelblue', alpha=0.8, linewidth=1.0)
 #   ax4.set_ylim(-0.5,0.5)
    ax4.set_ylabel('top2-191 (HpM)')
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)
    ax4.spines['bottom'].set_visible(False)

    ax5.plot(cc4['Pos'], cc4['Wnorm'], color ='firebrick', alpha=0.8, linewidth=1.0)
    ax5.plot(cc4['Pos'], cc4['Cnorm'], color ='steelblue', alpha=0.8, linewidth=1.0)
  #  ax5.set_ylim(-2,2)
    ax5.set_ylabel('top2-191-30 (HpM)')
    ax5.spines['top'].set_visible(False)
    ax5.spines['right'].set_visible(False)
    ax5.spines['bottom'].set_visible(False)
    ax6.spines['top'].set_visible(False)
    ax6.spines['right'].set_visible(False)
    
    ax2.get_xaxis().set_ticks([])
    ax3.get_xaxis().set_ticks([])
    ax4.get_xaxis().set_ticks([])
    ax5.get_xaxis().set_ticks([])
  #  ax6.get_xaxis().set_ticks([])
 #   ax4.set_xlabel('Chromosome position')

 
    for c in centro.itertuples(index=False, name=None):
            ax6.axvspan(c[0],c[1],0.1,0.15,color="black",alpha=0.8)
            
    for t in telo.itertuples(index=False, name=None):
            ax6.axvspan(t[0],t[1],0.1,0.15,color="grey",alpha=0.8)
            
    for ge in genee.itertuples(index=False, name=None):
     #   ax6.annotate(ge[0], xy = [ge[4],0.45])  
        if ge[3] == 'protein coding':
            if ge[7] == 'reverse':
                ax6.axvspan(ge[4],ge[5],0.2,0.4,color="firebrick",alpha=0.3)
            elif ge[7] == 'forward':
                ax6.axvspan(ge[4],ge[5],0.7,0.9,color="steelblue",alpha=0.3)

    return ff

chromosome3 = Chromosome_plot(wtETchr3, wtET35chr3, t2Echr3, t2E30chr3, centro3, feat3, telo3)
chromosome1 = Chromosome_plot(t2tchr1, t2rchr1, t2Echr1, t2EBchr1, centro1, feat1, telo1)



def Chromosome_plot (cc1):
    ff, (ax2) = plt.subplots(1,1, sharex=True)
    #ax1.set_title('WT')
    ax2.bar(cc1['Pos'], cc1['Wnorm'], color ='firebrick', alpha=0.8)
    ax2.bar(cc1['Pos'], cc1['Cnorm'], color ='steelblue', alpha=0.8)
    ax2.set_xlim(120000,120010)
    ax2.set_ylabel('WT (HpM)')
    #ax2.set_ylim(-5,5)
    ax2.set_ylabel('WT')


    return ff

chromosome1 = Chromosome_plot(wtETchr1)





def Chromosome_plot (cc1, cc2, cc3, cc4, centro):
    ff, (ax2, ax3, ax4, ax5) = plt.subplots(4,1)
    #ax1.set_title('WT')
    ax2.plot(cc1['Pos'], cc1['Wnorm'], color ='firebrick', alpha=0.8)
    ax2.plot(cc1['Pos'], cc1['Cnorm'], color ='steelblue', alpha=0.8)
    ax2.set_ylim(-2,2)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
   # ax2.set_ylabel('WT')
    #ax2.set_ylim(-5,5)
    ax2.set_ylabel('WT_A')

    ax3.plot(cc2['Pos'], cc2['Wnorm'], color ='firebrick', alpha=0.8)
    ax3.plot(cc2['Pos'], cc2['Cnorm'], color ='steelblue', alpha=0.8)
    ax3.set_ylim(-2,2)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.spines['bottom'].set_visible(False)
   # ax3.set_ylabel('SEN1d (HpM)')
    ax3.set_ylabel('WT_B')


    ax4.plot(cc3['Pos'], cc3['Wnorm'], color ='firebrick', alpha=0.8)
    ax4.plot(cc3['Pos'], cc3['Cnorm'], color ='steelblue', alpha=0.8)
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)
    ax4.spines['bottom'].set_visible(False)
    ax4.set_ylim(-2,2)
    ax4.set_ylabel('Top2-191_A')

    ax5.plot(cc4['Pos'], cc4['Wnorm'], color ='firebrick', alpha=0.8)
    ax5.plot(cc4['Pos'], cc4['Cnorm'], color ='steelblue', alpha=0.8)
    ax5.set_ylim(-2,2)
    ax5.set_ylabel('Top2-191_B')
   # ax6.set_ylabel('Gene Annotations')
    ax5.set_xlabel('Chromosome position')
    ax5.spines['top'].set_visible(False)
    ax5.spines['right'].set_visible(False)

    ax2.get_xaxis().set_ticks([])
    ax3.get_xaxis().set_ticks([])
    ax4.get_xaxis().set_ticks([])

 
    for c in centro.itertuples(index=False, name=None):
            ax2.axvspan(c[0],c[1],0.1,0.15,color="black",alpha=0.8)
            ax3.axvspan(c[0],c[1],0.1,0.15,color="black",alpha=0.8)
            ax4.axvspan(c[0],c[1],0.1,0.15,color="black",alpha=0.8)
            ax5.axvspan(c[0],c[1],0.1,0.15,color="black",alpha=0.8)
            
   
    return ff

chromosome1 = Chromosome_plot(wtETchr1, wtETBchr1, t2Echr1, t2EBchr1, centro1)
chromosome2 = Chromosome_plot(wtETchr2, wtETBchr2, t2Echr2, t2EBchr2, centro2)
chromosome1 = Chromosome_plot(wtETchr3, wtETBchr3, t2Echr3, t2EBchr3, centro3)



#For the averages
def Chromosome_plot (cc1, cc2, centro, genee, telo):
    ff, (ax2, ax3, ax6) = plt.subplots(3,1, sharex=True)
    #ax1.set_title('WT')
    ax2.plot(cc1['Pos'], cc1['Wnorm'], color ='firebrick', alpha=0.8, linewidth=1.0)
    ax2.plot(cc1['Pos'], cc1['Cnorm'], color ='steelblue', alpha=0.8 , linewidth=1.0)
    #ax2.set_ylim(-2,2)
    ax2.set_ylabel('WT')
    ax2.set_ylim(-0.25,0.25)
    ax2.set_ylabel('WT')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
  #  ax2.get_yaxis().set_ticks([])
  #  ax2.spines['left'].set_visible(False)
  #E55807'
  #E55807'

    ax3.plot(cc2['Pos'], cc2['Wnorm'], color ='firebrick', alpha=0.8, linewidth=1.0)
    ax3.plot(cc2['Pos'], cc2['Cnorm'], color ='steelblue', alpha=0.8, linewidth=1.0)
  #  ax3.set_ylim(-2,2)
    ax3.set_ylim(-0.25,0.25)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.spines['bottom'].set_visible(False)
 #   ax3.get_yaxis().set_ticks([])
  #  ax3.spines['left'].set_visible(False)
   # ax3.set_ylabel('SEN1d (HpM)')
    ax3.set_ylabel('Top2-191')

    ax6.set_xlabel('Chromosome position')


    ax6.set_ylabel('Gene Annotations')
    ax6.set_xlabel('Chromosome position')
 #   ax2.get_xaxis().set_ticks([])
 #   ax3.get_xaxis().set_ticks([])
    ax6.spines['top'].set_visible(False)
    ax6.spines['right'].set_visible(False)
 #   ax6.axvspan(660000,760000,0,1,color="pink",alpha=0.6)

 #   ax2.set_xlim(1550000,1700000)
  #  ax3.set_xlim(1550000,1700000)
   # ax6.set_xlim(1550000,1700000)
    

 
    for c in centro.itertuples(index=False, name=None):
            ax6.axvspan(c[0],c[1],0.1,0.15,color="black",alpha=0.8)
            
    for t in telo.itertuples(index=False, name=None):
            ax6.axvspan(t[0],t[1],0.1,0.15,color="grey",alpha=0.8)
            
    for ge in genee.itertuples(index=False, name=None):
        
        if ge[3] == 'protein coding':
            if ge[7] == 'reverse':
                ax6.axvspan(ge[4],ge[5],0.2,0.4,color="steelblue",alpha=0.3)
                #ax6.annotate(ge[0], xy = [ge[4],0.45]) 
            elif ge[7] == 'forward':
                ax6.axvspan(ge[4],ge[5],0.7,0.9,color="firebrick",alpha=0.3)
                #ax6.annotate(ge[0], xy = [ge[4],0.65]) 

    return ff

chromosome1 = Chromosome_plot(wtETchr1, t2Echr1, centro1, feat1, telo1)
chromosome2 = Chromosome_plot(wtETchr2, t2Echr2, centro2, feat2, telo2)
chromosome1 = Chromosome_plot(wtETchr3, t2Echr3, centro3, feat3, telo3)

"#1D5B79", "#468B97", "#EF6262", "#F3AA60"


#try to make a plot for centromeres of all three chromosomes
#%%

#For the averages
def Chromosome_plot (cc1, cc1t, cc2, cc2t, cc3, cc3t, centroo, centroo2, centroo3):
    ff, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3,2)
    
    #chromosome 1    
    ax1.plot(cc1['Pos'], cc1['Wnorm'], color ='firebrick', alpha=0.8, linewidth=1.0)
    ax1.plot(cc1['Pos'], cc1['Cnorm'], color ='steelblue', alpha=0.8, linewidth=1.0)
    #ax2.set_ylim(-2,2)
  #  ax1.set_ylabel('WT')
    ax1.set_ylim(-0.5,0.5)
    ax1.set_title('WT')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    
    ax2.plot(cc1t['Pos'], cc1t['Wnorm'], color ='firebrick', alpha=0.8, linewidth=1.0)
    ax2.plot(cc1t['Pos'], cc1t['Cnorm'], color ='steelblue', alpha=0.8, linewidth=1.0)
    #ax2.set_ylim(-2,2)
  #  ax2.set_ylabel('WT')
    ax2.set_ylim(-0.5,0.5)
    ax2.set_title('top2-191')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
 #   ax2.spines['left'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    #ax2.spines['bottom'].set_visible(False)
 #   ax2.get_yaxis().set_ticks([])
  #  ax2.spines['left'].set_visible(False)
  
  
  
    

    
    ax3.plot(cc2['Pos'], cc2['Wnorm'], color ='firebrick', alpha=0.8, linewidth=1.0)
    ax3.plot(cc2['Pos'], cc2['Cnorm'], color ='steelblue', alpha=0.8, linewidth=1.0)
  #  ax3.set_ylim(-2,2)
    ax3.set_ylim(-0.5,0.5)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.spines['bottom'].set_visible(False)
  #  ax3.spines['bottom'].set_visible(False)
 #   ax3.get_yaxis().set_ticks([])
  #  ax3.spines['left'].set_visible(False)
   # ax3.set_ylabel('SEN1d (HpM)')
  #  ax3.set_ylabel('Top2-191')
  
    ax4.plot(cc2t['Pos'], cc2t['Wnorm'], color ='firebrick', alpha=0.8, linewidth=1.0)
    ax4.plot(cc2t['Pos'], cc2t['Cnorm'], color ='steelblue', alpha=0.8, linewidth=1.0)
  #  ax4.get_yaxis().set_ticks([])
  #  ax3.set_ylim(-2,2)
    ax4.set_ylim(-0.5,0.5)
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)
  #  ax4.spines['left'].set_visible(False)
    ax4.spines['bottom'].set_visible(False)
  #  ax4.spines['bottom'].set_visible(False)
    
    
    ax5.plot(cc3['Pos'], cc3['Wnorm'], color ='firebrick', alpha=0.8, linewidth=1.0)
    ax5.plot(cc3['Pos'], cc3['Cnorm'], color ='steelblue', alpha=0.8, linewidth=1.0)
  #  ax3.set_ylim(-2,2)
    ax5.set_ylim(-0.5,0.5)
    ax5.spines['top'].set_visible(False)
    ax5.spines['right'].set_visible(False)
    ax5.spines['bottom'].set_visible(False)
    
  #  ax5.spines['bottom'].set_visible(False)
 #   ax3.get_yaxis().set_ticks([])
  #  ax3.spines['left'].set_visible(False)
   # ax3.set_ylabel('SEN1d (HpM)')
  #  ax3.set_ylabel('Top2-191')
  
    ax6.plot(cc3t['Pos'], cc3t['Wnorm'], color ='firebrick', alpha=0.8, linewidth=1.0)
    ax6.plot(cc3t['Pos'], cc3t['Cnorm'], color ='steelblue', alpha=0.8, linewidth=1.0)
  #  ax3.set_ylim(-2,2)
    ax6.set_ylim(-0.5,0.5)
    ax6.spines['top'].set_visible(False)
    ax6.spines['right'].set_visible(False)
  #  ax6.spines['left'].set_visible(False)
    ax6.spines['bottom'].set_visible(False)
  #  ax6.get_yaxis().set_ticks([])
   # ax6.spines['bottom'].set_visible(False)
    for c in centroo.itertuples(index=False, name=None):
        ax1.axvspan(c[0],c[1],0.1,0.11,color="black",alpha=0.5)
        ax2.axvspan(c[0],c[1],0.1,0.11,color="black",alpha=0.5)

    for c2 in centroo2.itertuples(index=False, name=None):
        ax3.axvspan(c2[0],c2[1],0.1,0.11,color="black",alpha=0.5)
        ax4.axvspan(c2[0],c2[1],0.1,0.11,color="black",alpha=0.5)
    for c3 in centroo3.itertuples(index=False, name=None):
        ax5.axvspan(c3[0],c3[1],0.1,0.11,color="black",alpha=0.5)
        ax6.axvspan(c3[0],c3[1],0.1,0.11,color="black",alpha=0.5)





   # ax6.set_xlabel('Chromosome position')
 #   ax2.get_xaxis().set_ticks([])
 #   ax3.get_xaxis().set_ticks([])
    ax1.get_xaxis().set_ticks([])
    ax2.get_xaxis().set_ticks([])
    ax3.get_xaxis().set_ticks([])
    ax4.get_xaxis().set_ticks([])
    ax5.get_xaxis().set_ticks([])
    
    ax6.get_xaxis().set_ticks([])

    ax2.set_xlim(3718504,3824604)
    ax1.set_xlim(3718504,3824604)
    ax3.set_xlim(1570455,1676556)
    ax4.set_xlim(1570455,1676556)
    ax5.set_xlim(1050904,1157003)
    ax6.set_xlim(1050904,1157003)
    


    return ff

chromosome1 = Chromosome_plot(wtETchr1, t2Echr1, wtETchr2, t2Echr2, wtETchr3, t2Echr3, centro1, centro2, centro3)
chromosome2 = Chromosome_plot(wtETchr2, t2Echr2, centro2, feat2, telo2)
chromosome1 = Chromosome_plot(wtETchr3, t2Echr3, centro3, feat3, telo3)


def score(file, frame1, frame2, frame3, frame4):
    for i in range(len(file)):
        tempstart = file.iloc[i]["start"]
        tempend = file.iloc[i]['end']
        tempchro = file.iloc[i]["Chromosome"]
    #    rna = file.iloc[i]['rna']

        tempsubset = frame1.loc[(frame1['Pos'] >= tempstart) & (frame1['Pos'] <= tempend) & (frame1['Chr'] == tempchro)]
        tempsubsett = frame2.loc[(frame2['Pos'] >= tempstart) & (frame2['Pos'] <= tempend) & (frame2['Chr'] == tempchro)]
        tempsubset35 = frame3.loc[(frame3['Pos'] >= tempstart) & (frame3['Pos'] <= tempend) & (frame3['Chr'] == tempchro)]
        tempsubsett30 = frame4.loc[(frame4['Pos'] >= tempstart) & (frame4['Pos'] <= tempend) & (frame4['Chr'] == tempchro)]


        file.loc[file.index[i], 'thing'] = tempsubset['normscore'].sum()
        file.loc[file.index[i], 'epb'] = tempsubset['normscore'].sum() / len(tempsubset)
        
        file.loc[file.index[i], 'thingt'] = tempsubsett['normscore'].sum()
        file.loc[file.index[i], 'epbt'] = tempsubsett['normscore'].sum() / len(tempsubsett)
        
        file.loc[file.index[i], 'thing35'] = tempsubset35['normscore'].sum()
        file.loc[file.index[i], 'epb35'] = tempsubset35['normscore'].sum() / len(tempsubset35)
        
        file.loc[file.index[i], 'thingt30'] = tempsubsett30['normscore'].sum()
        file.loc[file.index[i], 'epbt30'] = tempsubsett30['normscore'].sum() / len(tempsubsett30)

        # Add log2 of epb and epbt as new columns
        file.loc[file.index[i], 'log2_epb'] = np.log2(file.loc[file.index[i], 'epb'])
        file.loc[file.index[i], 'log2_epbt'] = np.log2(file.loc[file.index[i], 'epbt'])

        # Calculate the length of the region and add log2_length column
        length = tempend - tempstart
        file.loc[file.index[i], 'log2_length'] = np.log2(length)
      #  file.loc[file.index[i], 'log2_rna'] = np.log2(rna)

    return file

comcentro = score(comcentro, wtET, t2E, wtET35, t2E30)

#need to subset first 
comcentrosubset = comcentro[["Chromosome", "thing", "thingt", 'thing35', 'thingt30']]
comcentrosubset = pd.melt(comcentrosubset, id_vars=['Chromosome'], value_vars=['thing', 'thingt', 'thing35', 'thingt30'])
import seaborn as sns
custom_palette = ["#1D5B79", "#EF6262", "#A6CF98", "#557C55" ]
sns.set_palette(custom_palette)
ax = sns.barplot(x="Chromosome", y="value", hue="variable", data=comcentrosubset)

#%%



normscore_df1 = bwtET['score']
normscore_df2 = bwtET35['score']

normscore_df3 = bt2E['score']
normscore_df4 = bt2E30['score']

# Create a new DataFrame with the selected 'normscore' columns
new_df = pd.DataFrame({'WT_permissive': normscore_df1, 'WT_restrictive': normscore_df2, 'top2-191_restrictive': normscore_df3, 'top2-191_permissive': normscore_df4})

iris_corr_matrix = new_df.corr()
print(iris_corr_matrix)
import seaborn as sns
# Create the heatmap using the `heatmap` function of Seaborn
sns.heatmap(iris_corr_matrix, cmap='coolwarm', annot=True, vmin= 0.60, vmax = 1)

# Display the heatmap using the `show` method of the `pyplot` module from matplotlib.
plt.show()


new_df['WT_A'] = np.log2(new_df['normscore_df1'])
new_df['WT_B'] = np.log2(new_df['normscore_df2'])
new_df['T2_A'] = np.log2(new_df['normscore_df3'])
new_df['T2_B'] = np.log2(new_df['normscore_df4'])

df_dropped = new_df.drop("normscore_df1", axis=1)
df_dropped = df_dropped.drop("normscore_df2", axis=1)
df_dropped = df_dropped.drop("normscore_df3", axis=1)
df_dropped = df_dropped.drop("normscore_df4", axis=1)

import seaborn as sns
sns.pairplot(new_df, diag_kind="hist")

correlation = wtET[column_df1].corr(t2E[column_df2])

print(new_df)

g = sns.PairGrid(df_dropped)
#g.map_diag(sns.hist, hue=None, color=".3")
g.map_offdiag(sns.scatterplot)
g.add_legend()

f = sns.pairplot(df_dropped, plot_kws=dict(marker="+", linewidth=0.5), diag_kind="hist",corner=True)




from sklearn.linear_model import LinearRegression

x = np.array(wtET['normscore']).reshape((-1,1))
y = np.array(wtETB['normscore'])
#model = LinearRegression()
#model.fit(x, y)
model = LinearRegression().fit(x, y)
r_sq = model.score(x, y)
print('coefficient of determination:', r_sq)

xx = np.array(t2E['normscore']).reshape((-1,1))
yy = np.array(t2EB['normscore'])
#model = LinearRegression()
#model.fit(x, y)
modell = LinearRegression().fit(xx, yy)
r_sqC = model.score(xx, yy)
print('coefficient of determination:', r_sqC)



xt = np.array(wtET['normscore']).reshape((-1,1))
yt = np.array(t2E['normscore'])
#model = LinearRegression()
#model.fit(x, y)
modelt = LinearRegression().fit(xt, yt)
r_sqt = model.score(xt, yt)
print('coefficient of determination:', r_sqt)


xxt = np.array(t2EB['normscore']).reshape((-1,1))
yyt = np.array(wtETB['normscore'])
#model = LinearRegression()
#model.fit(x, y)
modellt = LinearRegression().fit(xxt, yyt)
r_sqCt = model.score(xxt, yyt)
print('coefficient of determination:', r_sqCt)


#%%


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
centrowtET1 = centrocheck(centro1,wtETchr1)
centrowtET2 = centrocheck(centro2,wtETchr2)
centrowtET3 = centrocheck(centro3,wtETchr3)

centrowtETB1 = centrocheck(centro1,wtETBchr1)
centrowtETB2 = centrocheck(centro2,wtETBchr2)
centrowtETB3 = centrocheck(centro3,wtETBchr3)

centrowtETC1 = centrocheck(centro1,wtETCchr1)
centrowtETC2 = centrocheck(centro2,wtETCchr2)
centrowtETC3 = centrocheck(centro3,wtETCchr3)


centot2E1 = centrocheck(centro1,t2Echr1)
centot2E2 = centrocheck(centro2,t2Echr2)
centot2E3 = centrocheck(centro3,t2Echr3)

centot2EB1 = centrocheck(centro1,t2EBchr1)
centot2EB2 = centrocheck(centro2,t2EBchr2)
centot2EB3 = centrocheck(centro3,t2EBchr3)

cenform1 = centrocheck(centro1,formchr1)
cenform2 = centrocheck(centro2,formchr2)
cenform3 = centrocheck(centro3,formchr3)



#%%

def calculate_mean_and_rolling_average(centot2t1, centot2t2, centot2t3, window_size=500):
    # Calculate middle indices for each DataFrame
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

    # Calculate the rolling window average
    rolling_average_t2t = mean_t2tcentro['Mean'].rolling(window=window_size, center=True).mean()

    # Create a new DataFrame with the rolling average values
    rolling_average_df_t2t = pd.DataFrame({'Rolling Average': rolling_average_t2t})

    return mean_t2tcentro, rolling_average_df_t2t

# Example usage:
# Assuming you have defined the required DataFrames (centot2t1, centot2t2, centot2t3) and optionally window_size
# mean_t2tcentro, rolling_average_df_t2t = calculate_mean_and_rolling_average(centot2t1, centot2t2, centot2t3, window_size=3)
# print(mean_t2tcentro)
# print(rolling_average_df_t2t)

mean_t2Ecentro, rolling_t2Ecentro = calculate_mean_and_rolling_average(centot2E1, centot2E2, centot2E3, window_size=500)
mean_t2EBcentro, rolling_t2EBcentro = calculate_mean_and_rolling_average(centot2EB1, centot2EB2, centot2EB3, window_size=500)

mean_formcentrro, rolling_formcentrro = calculate_mean_and_rolling_average(cenform1, cenform2, cenform3, window_size=500)

mean_WTETcentro, rolling_WTETcentro = calculate_mean_and_rolling_average(centrowtET1, centrowtET2, centrowtET3, window_size=500)
mean_WTETBcentro, rolling_WTETBcentro = calculate_mean_and_rolling_average(centrowtETB1, centrowtETB2, centrowtETB3, window_size=500)


"#1D5B79", "#468B97", "#EF6262", "#F3AA60"

fig, ((ax1)) = plt.subplots(1,1, sharex=True)
#ax1.plot(rolling_formcentrro.index.values, rolling_formcentrro['Rolling Average'], color='black')
ax1.plot(rolling_t2Ecentro.index.values, rolling_t2Ecentro['Rolling Average'], color='#E14D2A')
#ax1.plot(rolling_t2EBcentro.index.values, rolling_t2EBcentro['Rolling Average'], color='crimson')
ax1.plot(rolling_WTETcentro.index.values, rolling_WTETcentro['Rolling Average'], color='#001253')
#ax1.plot(rolling_WTETBcentro.index.values, rolling_WTETBcentro['Rolling Average'], color='steelblue')

#ax1.set_xlabel('midpoint alligned centromeric signal (relative bp)')
ax1.set_ylabel('Aggregated CC-seq Signal (500bp smoothed)(HpM)')
ax1.set_xlabel('midpoint alligned centromeric signal (relative bp)')
#ax2.set_ylabel('Aggregated CC-seq Signal (100bp smoothed)(HpM)')
ax1.set_xticks([0,16525,33050,49575,66100])
ax1.set_xticklabels(['-33050','-16525','0','+16525','+33050'])
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
#ax1.set_title('Wild Type')
#ax2.set_title('Top2-191t')



#%%

#Let's try to adapt our feat list

def Findfeat(file):
    genes = pd.read_csv(file, delimiter="\t")
    print(genes)
    genes.loc[genes['Chromosome'] == "I", 'chro'] = 'chr1'
    genes.loc[genes['Chromosome'] == "II", 'chro'] = 'chr2'
    genes.loc[genes['Chromosome'] == "III", 'chro'] = 'chr3'
    
    featfor = genes[(genes['Strand'] == 'forward')]
    featrev = genes[(genes['Strand'] == 'reverse')]
    
    return featfor, featrev, genes

listfor, istrev, listy = Findfeat('protein_coding_gene_list.tsv')

listy.drop('Product description', axis=1, inplace=True)
listy.drop('Chromosome', axis=1, inplace=True)
listy.drop('Feature type', axis=1, inplace=True)
  

listone = listy[(listy['chro'] == 'chr1')]
list1s = list1.sort_values(by='star')
list2 = listy[(listy['chro'] == 'chr2')]
list3 = listy[(listy['chro'] == 'chr3')]





def check_plot (centro, genee, telo, ggr):
    ff, (ax6) = plt.subplots(1,1, sharex=True)

    ax6.set_ylabel('Gene Annotations')
    ax6.set_xlabel('Chromosome position')

 
    for c in centro.itertuples(index=False, name=None):
            ax6.axvspan(c[0],c[1],0.1,0.15,color="black",alpha=0.8)
            
    for t in telo.itertuples(index=False, name=None):
            ax6.axvspan(t[0],t[1],0.1,0.15,color="grey",alpha=0.8)
            
    for ge in genee.itertuples(index=False, name=None):
        
        if ge[3] == 'protein coding':
            if ge[7] == 'reverse':
                ax6.axvspan(ge[4],ge[5],0.2,0.4,color="firebrick",alpha=0.3)
                #ax6.annotate(ge[0], xy = [ge[4],0.45]) 
            elif ge[7] == 'forward':
                ax6.axvspan(ge[4],ge[5],0.7,0.9,color="steelblue",alpha=0.3)
                #ax6.annotate(ge[0], xy = [ge[4],0.65]) 
                
    for gg in ggr.itertuples(index=False, name=None):

        ax6.axvspan(gg[0],gg[1],0.4,0.6,color="green",alpha=0.3)

    return ff

chromosome1 = check_plot(centro1, feat1, telo1,igr)





#%%%



def check_plot (centro, genee, telo, ggr):
    ff, (ax6) = plt.subplots(1,1, sharex=True)

    ax6.set_ylabel('Gene Annotations')
    ax6.set_xlabel('Chromosome position')

 
    for c in centro.itertuples(index=False, name=None):
            ax6.axvspan(c[0],c[1],0.1,0.15,color="black",alpha=0.8)
            
    for t in telo.itertuples(index=False, name=None):
            ax6.axvspan(t[0],t[1],0.1,0.15,color="grey",alpha=0.8)
            
    for ge in genee.itertuples(index=False, name=None):
        
        if ge[3] == 'protein coding':
            if ge[7] == 'reverse':
                ax6.axvspan(ge[4],ge[5],0.2,0.4,color="firebrick",alpha=0.3)
                #ax6.annotate(ge[0], xy = [ge[4],0.45]) 
            elif ge[7] == 'forward':
                ax6.axvspan(ge[4],ge[5],0.7,0.9,color="steelblue",alpha=0.3)
                #ax6.annotate(ge[0], xy = [ge[4],0.65]) 
                
    for gg in ggr.itertuples(index=False, name=None):

        ax6.axvspan(gg[0],gg[1],0.4,0.6,color="green",alpha=0.3)

    return ff

chromosome1 = check_plot(centro1, feat1, telo1,igr_all)



#%%





def check_plot (centro, genee, telo, ggr):
    ff, (ax6) = plt.subplots(1,1, sharex=True)

    ax6.set_ylabel('Gene Annotations')
    ax6.set_xlabel('Chromosome position')

 
    for c in centro.itertuples(index=False, name=None):
            ax6.axvspan(c[0],c[1],0.1,0.15,color="black",alpha=0.8)
            
    for t in telo.itertuples(index=False, name=None):
            ax6.axvspan(t[0],t[1],0.1,0.15,color="grey",alpha=0.8)
            
    for ge in genee.itertuples(index=False, name=None):
        
        if ge[3] == 'protein coding':
            if ge[7] == 'reverse':
                ax6.axvspan(ge[4],ge[5],0.2,0.4,color="firebrick",alpha=0.3)
                #ax6.annotate(ge[0], xy = [ge[4],0.45]) 
            elif ge[7] == 'forward':
                ax6.axvspan(ge[4],ge[5],0.7,0.9,color="steelblue",alpha=0.3)
                #ax6.annotate(ge[0], xy = [ge[4],0.65]) 
                
    for gg in ggr.itertuples(index=False, name=None):
       # if gg[3] == 'Tandem':
        ax6.axvspan(gg[0],gg[1],0.4,0.6,color="pink",alpha=0.3)
   #     if gg[3] == 'Divergent':
    #            ax6.axvspan(gg[0],gg[1],0.4,0.6,color="yellow",alpha=0.3)
     #   if gg[3] == 'Convergent':
      #          ax6.axvspan(gg[0],gg[1],0.4,0.6,color="green",alpha=0.3)

      #  ax6.axvspan(gg[0],gg[1],0.4,0.6,color="green",alpha=0.3)

    return ff

chromosome1 = check_plot(centro1, feat1, telo1,igr_all)




def find_intergenic_regions_extra(list1, centro, ends):
    # Initialize lists to store IGR data
    igr_starts = []
    igr_ends = []
    igr_chros = []
    igr_types = []  # List to store IGR types (tandem, convergent, divergent)

    # Group the DataFrame by 'chro' (chromosome) to handle multiple chromosomes
    grouped_list1 = list1.groupby('chro')

    # Iterate through each chromosome group to find IGRs
    for chro, group in grouped_list1:
        # Sort the DataFrame by 'Start position' in ascending order for the current chromosome
        group_sorted = group.sort_values(by='Start position')

        # Iterate through the sorted DataFrame to find IGRs for the current chromosome
        prev_end = -1
        for i, row in group_sorted.iterrows():
            start_pos = row['Start position']
            end_pos = row['End position']
            strand = row['Strand']

            # Check if there is a gap between the previous gene and the current gene
            if start_pos  prev_end + 1:
                igr_starts.append(prev_end + 1)
                igr_ends.append(start_pos - 1)
                igr_chros.append(chro)
                igr_types.append('Tandem')  # IGR between two genes on the same strand is tandem

            prev_end = max(prev_end, end_pos)

            # Check if the current gene and the next gene are on different strands
            if i < len(group_sorted) - 1:
                next_row = group_sorted.iloc[i + 1]
                next_strand = next_row['Strand']
                if strand == 'forward' and next_strand == 'reverse':
                    igr_types[-1] = 'Convergent'  # IGR between forward and reverse genes is convergent
                elif strand == 'reverse' and next_strand == 'forward':
                    igr_types[-1] = 'Divergent'  # IGR between reverse and forward genes is divergent

            # Modify the IGR end position based on IGR type
            if igr_types[-1] == 'Tandem':
                # IGR between two genes on the same strand (tandem) - from stop codon to the next start codon
                igr_ends[-1] = group_sorted.iloc[i + 1]['Start position'] - 1 if i + 1 < len(group_sorted) else end_pos
            elif igr_types[-1] == 'Convergent':
                # IGR between forward and reverse genes (convergent) - from start codon to the next start codon
                igr_ends[-1] = group_sorted.iloc[i + 1]['Start position'] - 1 if i + 1 < len(group_sorted) else end_pos
            elif igr_types[-1] == 'Divergent':
                # IGR between reverse and forward genes (divergent) - from stop codon to the start codon
                igr_ends[-1] = group_sorted.iloc[i + 1]['End position'] if i + 1 < len(group_sorted) else end_pos

        # Check if the last gene ends before the chromosome ends
        if prev_end < group_sorted['End position'].max():
            igr_starts.append(prev_end + 1)
            igr_ends.append(group_sorted['End position'].max())
            igr_chros.append(chro)
            igr_types.append('Tandem')  # IGR between the last gene and chromosome end is tandem

    # Create a new DataFrame with IGR data
    igr_df = pd.DataFrame({
        'Start position': igr_starts,
        'End position': igr_ends,
        'chro': igr_chros,
        'IGR Type': igr_types  # Add the IGR Type column
    })

    # Exclude IGRs that overlap with centromeres
    if not centro.empty:
        valid_igr = ~igr_df.apply(lambda row: centro.apply(lambda cent_row:
                                    cent_row['start'] <= row['End position'] and
                                    cent_row['end'] >= row['Start position'] and
                                    cent_row['chro'] == row['chro'], axis=1).any(), axis=1)
        igr_df = igr_df[valid_igr]

    # Exclude IGRs that overlap with "ends" DataFrame
    if not ends.empty:
        valid_igr = ~igr_df.apply(lambda row: ends.apply(lambda end_row:
                                    end_row['start'] <= row['End position'] and
                                    end_row['end'] >= row['Start position'] and
                                    end_row['chro'] == row['chro'], axis=1).any(), axis=1)
        igr_df = igr_df[valid_igr]

    return igr_df

igr_all_new = find_intergenic_regions_extra(listy, centro_merged, ends_merged)
#%%



def calculate_igr_length_and_quartile(igr_df):
    # Calculate the length of each IGR
    igr_df['Length'] = igr_df['End position'] - igr_df['Start position'] + 1

    # Find the quartile each IGR lies in based on its length
    quartile_25 = igr_df['Length'].quantile(0.25)
    quartile_50 = igr_df['Length'].quantile(0.50)
    quartile_75 = igr_df['Length'].quantile(0.75)

    def find_quartile_category(length):
        if length <= quartile_25:
            return 'Q1'
        elif length <= quartile_50:
            return 'Q2'
        elif length <= quartile_75:
            return 'Q3'
        else:
            return 'Q4'

    igr_df['Quartile'] = igr_df['Length'].apply(find_quartile_category)
    igr_df.loc[igr_df['chro'] == "chr1", 'Chr'] = 1
    igr_df.loc[igr_df['chro'] == "chr2", 'Chr'] = 2
    igr_df.loc[igr_df['chro'] == "chr3", 'Chr'] = 3
    igr_df['Chr'] = igr_df['Chr'].astype(int)

    return igr_df

igrst = calculate_igr_length_and_quartile(igr_all)




qqqQ1IGR = igrst[(igrst['Quartile'] == 'Q1')]
qqqQ2IGR = igrst[(igrst['Quartile'] == 'Q2')]
qqqQ3IGR = igrst[(igrst['Quartile'] == 'Q3')]
qqqQ4IGR = igrst[(igrst['Quartile'] == 'Q4')]

#%%

def centrocheck(file,frame1):
    data_output = pd.DataFrame()

    for i in range(len(file)):
        tempstart = file.iloc[i]["Start position"]
        
        tempend = file.iloc[i]["End position"]
        tempchro = file.iloc[i]["Chr"]
       # print(tempchro)
        
        lenny = tempend - tempstart
       # print(frame1.dtypes['Pos'])  # Print the data type of the 'Pos' column
        #print(frame1.dtypes['Chr'])  # Print the data type of the 'Chr' column

# Assuming tempstart, tempend, and tempchro are variables
        #print(type(tempstart))  # Print the data type of tempstart
       # print(type(tempend))    # Print the data type of tempend
       # print(type(tempchro)) 

        tempsubsetfor = frame1.loc[(frame1['Pos']>= tempstart) & (frame1['Pos'] <= tempend) & (frame1['Chr'] == tempchro)]
        print(tempsubsetfor)
        tempsubsetfor = tempsubsetfor.reset_index(drop = True)
        
        
        line = pd.DataFrame({"Pos": tempend}, index=[lenny])
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
    
            "aPos" : np.arange(0.0, lenny, 1)}).merge(
        
                tempsubsetfor, on = "aPos", how = "left"

                )

        data_blank = data_blank[["aPos", "normscore"]]
            
        data_output = pd.concat([data_output, data_blank['normscore']], axis =1, ignore_index=True)
          
            
        data_output = data_output.fillna(0)



    return data_output
centrowtET1Q1 = centrocheck(qqqQ1IGR,wtET)
centrowtET1Q2 = centrocheck(qqqQ2IGR,wtET)
centrowtET1Q3 = centrocheck(qqqQ3IGR,wtET)
centrowtET1Q4 = centrocheck(qqqQ4IGR,wtET)


def calculate_mean_and_rolling_average2(dataframe, window_size=10):
    # Calculate middle index for the DataFrame
    middle_index = dataframe.shape[0] // 2
    print(middle_index)

    # Calculate the maximum length of the DataFrame
    max_length = dataframe.shape[0]

    # Create a new DataFrame with centered index
    offset = max_length // 2 - middle_index
    centered_df = dataframe.copy()
    centered_df.index = range(offset, offset + dataframe.shape[0])

    # Calculate the mean values
    mean_values = centered_df.mean(axis=1)

    # Create a new DataFrame with the mean values
    mean_t2tcentro = pd.DataFrame(mean_values, columns=['Mean'])

    # Display the mean DataFrame
    print(mean_t2tcentro)

    # Calculate the rolling window average
    rolling_average_t2t = mean_t2tcentro['Mean'].rolling(window=window_size, center=True).mean()

    # Create a new DataFrame with the rolling average values
    rolling_average_df_t2t = pd.DataFrame({'Rolling Average': rolling_average_t2t})

    return mean_t2tcentro, rolling_average_df_t2t

mean_centrowtET1Q1, rolling_centrowtET1Q1 = calculate_mean_and_rolling_average2 (centrowtET1Q1, window_size =10)
mean_centrowtET1Q2, rolling_centrowtET1Q2 = calculate_mean_and_rolling_average2 (centrowtET1Q2, window_size =10)
mean_centrowtET1Q3, rolling_centrowtET1Q3 = calculate_mean_and_rolling_average2 (centrowtET1Q3, window_size =10)
mean_centrowtET1Q4, rolling_centrowtET1Q4 = calculate_mean_and_rolling_average2 (centrowtET1Q4, window_size =10)





fig, ((ax1)) = plt.subplots(1,1, sharex=True)
ax1.plot(rolling_centrowtET1Q1.index.values, rolling_centrowtET1Q1['Rolling Average'], color='darkslategray')
ax1.plot(rolling_centrowtET1Q2.index.values, rolling_centrowtET1Q2['Rolling Average'], color='teal')
ax1.plot(rolling_centrowtET1Q3.index.values, rolling_centrowtET1Q3['Rolling Average'], color='firebrick')
ax1.plot(rolling_centrowtET1Q4.index.values, rolling_centrowtET1Q4['Rolling Average'], color='lightcoral')


#ax1.set_xlabel('midpoint alligned centromeric signal (relative bp)')
ax1.set_ylabel('Aggregated CC-seq Signal (100bp smoothed)(HpM)')
ax1.set_xlabel('midpoint alligned IGR signal (relative bp)')
#ax2.set_ylabel('Aggregated CC-seq Signal (100bp smoothed)(HpM)')
ax1.set_xticks([0,9106])
ax1.set_xticklabels(['IGR start','+IGR end'])
#ax1.set_title('Wild Type')
#ax2.set_title('Top2-191t')


fig, ax1 = plt.subplots(1, 1, sharex=True)

# Use the offset to center the data on the plot
ax1.plot(rolling_centrowtET1Q1.index.values - len(rolling_centrowtET1Q1) // 2, rolling_centrowtET1Q1['Rolling Average'], color='darkslategray')
ax1.plot(rolling_centrowtET1Q2.index.values - len(rolling_centrowtET1Q2) // 2, rolling_centrowtET1Q2['Rolling Average'], color='teal')
ax1.plot(rolling_centrowtET1Q3.index.values - len(rolling_centrowtET1Q3) // 2, rolling_centrowtET1Q3['Rolling Average'], color='firebrick')
ax1.plot(rolling_centrowtET1Q4.index.values - len(rolling_centrowtET1Q4) // 2, rolling_centrowtET1Q4['Rolling Average'], color='lightcoral')

ax1.set_ylabel('Aggregated CC-seq Signal (100bp smoothed)(HpM)')
ax1.set_xlabel('midpoint aligned IGR signal (relative bp)')
ax1.set_xticks([-len(rolling_centrowtET1Q1) // 2, -len(rolling_centrowtET1Q1) // 4, 0, len(rolling_centrowtET1Q1) // 4, len(rolling_centrowtET1Q1) // 2])
#ax1.set_xticklabels(['-33050', '-16525', '0', '+16525', '+33050'])

# Add legends and other customizations as needed
plt.show()

#%%

def score(file, frame1, frame2):
    for i in range(len(file)):
        tempstart = file.iloc[i]["Start position"]
        tempend = file.iloc[i]['End position']
        tempchro = file.iloc[i]["Chr"]
        rna = file.iloc[i]['rna']

        tempsubset = frame1.loc[(frame1['Pos'] >= tempstart) & (frame1['Pos'] <= tempend) & (frame1['Chr'] == tempchro)]
        tempsubsett = frame2.loc[(frame2['Pos'] >= tempstart) & (frame2['Pos'] <= tempend) & (frame2['Chr'] == tempchro)]

        file.loc[file.index[i], 'thing'] = tempsubset['normscore'].sum()
        file.loc[file.index[i], 'epb'] = tempsubset['normscore'].sum() / len(tempsubset)
        
        file.loc[file.index[i], 'thingt'] = tempsubsett['normscore'].sum()
        file.loc[file.index[i], 'epbt'] = tempsubsett['normscore'].sum() / len(tempsubsett)

        # Add log2 of epb and epbt as new columns
        file.loc[file.index[i], 'log2_epb'] = np.log2(file.loc[file.index[i], 'epb'])
        file.loc[file.index[i], 'log2_epbt'] = np.log2(file.loc[file.index[i], 'epbt'])

        # Calculate the length of the region and add log2_length column
        length = tempend - tempstart
        file.loc[file.index[i], 'log2_length'] = np.log2(length)
        file.loc[file.index[i], 'log2_rna'] = np.log2(rna)

    return file

dataIGR = score(megeed_dye, wtET, t2E)
#


dataIGR['epb_quartile'] = pd.qcut(dataIGR['epb'], q=4, labels=False)
dataIGR['length_quartile'] = pd.qcut(dataIGR['log2_length'], q=4, labels=False)

dataQ1IGR = dataIGR[(dataIGR['Quartile'] == 'Q1')]
dataQ2IGR = dataIGR[(dataIGR['Quartile'] == 'Q2')]
dataQ3IGR = dataIGR[(dataIGR['Quartile'] == 'Q3')]
dataQ4IGR = dataIGR[(dataIGR['Quartile'] == 'Q4')]

import seaborn as sns

melted_df = pd.melt(dataIGR, id_vars=['IGR Type'], value_vars=['epb', 'epbt'], var_name='source', value_name='score')
melted_df['log2_score'] = np.log2(melted_df['score'])
# Print the new DataFrame
print(melted_df)

sns.catplot(data=melted_df, x="IGR Type", y="log2_score", kind="box", hue= 'source')

sns.jointplot(
    data=dataIGR,
    x="log2_length", y="log2_epb", hue="IGR Type",
    kind = 'kde'
)

ff, axes = plt.subplots(4,1)




sns.jointplot(
    data=dataIGR,
    x="log2_length", y="log2_epb", hue="epb_quartile",
    kind = 'kde'
)

ff, axes = plt.subplots(4,1)

sns.jointplot(data=dataQ1IGR, x="log2_length", y="epb", ax=axes[0])
sns.jointplot(data=dataQ2IGR, x="log2_length", y="log2_epb", ax=axes[1])
sns.jointplot(data=dataQ3IGR, x="log2_length", y="log2_epb", ax=axes[2])
sns.jointplot(data=dataQ4IGR, x="log2_length", y="log2_epb", ax=axes[3])


#girly pop needs to do the same thing for genes
#seperate by size or by transcription?
#positions of g4s or some shit? 
#do I consider the replication and transcription direction?
#write into the pipeline convergent vs divergent vs concordant igrs

#%%%

    igr_df.loc[igr_df['chro'] == "chr1", 'Chr'] = 1
    igr_df.loc[igr_df['chro'] == "chr2", 'Chr'] = 2
    igr_df.loc[igr_df['chro'] == "chr3", 'Chr'] = 3

def Findfeat(file):
    genes = pd.read_csv(file, delimiter="\t")
    print(genes)
    genes.loc[genes['Chromosome'] == "I", 'Chr'] = 1
    genes.loc[genes['Chromosome'] == "II", 'Chr'] = 2
    genes.loc[genes['Chromosome'] == "III", 'Chr'] = 3
    genes.dropna(inplace=True)

    genes['Chr'] = genes['Chr'].astype(int)
    
    featfor = genes[(genes['Strand'] == 'forward')]
    featrev = genes[(genes['Strand'] == 'reverse')]
    
    return featfor, featrev, genes

featfor, featrev, ffeat = Findfeat('protein_coding_gene_list.tsv')




Genesy = pd.read_csv('Schizosaccharomyces_pombe_all_chromosomes.gff3',sep='\t', names=('1', '2','3', '4','5', '6','7', '8', '9'),skiprows=1)
unique_items_array = Genesy['3'].unique()
print(unique_items_array)



mask = Genesy['3'] == 'gene'
Genes = Genesy[mask]

masko = Genesy['3'] == 'origin_of_replication'
ori = Genesy[masko]

maskprro = Genesy['3'] == 'promoter'
promoter = Genesy[maskprro]


mask5u = Genesy['3'] == 'five_prime_UTR'
utr5 = Genesy[mask5u]


maskpr = Genesy['3'] == 'rRNA'
rRNA = Genesy[maskpr]

maskpnc = Genesy['3'] == 'ncRNA'
ncRNA = Genesy[maskpnc]


maskpsn = Genesy['3'] == 'snRNA'
snRNA = Genesy[maskpsn]


masksno = Genesy['3'] == 'snoRNA'
snoRNA = Genesy[masksno]

maskt = Genesy['3'] == 'tRNA'
tRNA = Genesy[maskt]




for line in Genes.itertuples():
     index = line[0]
     infor = line[9]
     position = infor.rfind(';')
    # print(position)
     if position >0:
         Genes.at[index, '9'] = infor[3:position]
     else:
         Genes.at[index, '9'] = infor[3:]

Genes.drop(labels = ['2', '3', '6', '8'], inplace =True, axis = 1)
Genes.rename(columns={'9':'Systematic ID', '4':'Start position', '5':'End position', '7':'strand', '1':'Chroo'}, inplace = True)
Genes.loc[Genes['Chroo'] == "I", 'Chr'] = 1
Genes.loc[Genes['Chroo'] == "II", 'Chr'] = 2
Genes.loc[Genes['Chroo'] == "III", 'Chr'] = 3
Genes.loc[Genes['strand'] == "-", 'Strand'] = 'reverse'
Genes.loc[Genes['strand'] == "+", 'Strand'] = 'forward'




for line in tRNA.itertuples():
     index = line[0]
     infor = line[9]
     position = infor.rfind(';')
    # print(position)
     if position >0:
         tRNA.at[index, '9'] = infor[3:position]
     else:
         tRNA.at[index, '9'] = infor[3:]

tRNA.drop(labels = ['2', '3', '6', '8'], inplace =True, axis = 1)
tRNA.rename(columns={'9':'ID', '4':'Start position', '5':'End position', '7':'strand', '1':'Chroo'}, inplace = True)
tRNA.loc[tRNA['Chroo'] == "I", 'Chr'] = 1
tRNA.loc[tRNA['Chroo'] == "II", 'Chr'] = 2
tRNA.loc[tRNA['Chroo'] == "III", 'Chr'] = 3
tRNA.loc[tRNA['strand'] == "-", 'Strand'] = 'reverse'
tRNA.loc[tRNA['strand'] == "+", 'Strand'] = 'forward'


for line in promoter.itertuples():
     index = line[0]
     infor = line[9]
     position = infor.rfind(';')
    # print(position)
     if position >0:
         promoter.at[index, '9'] = infor[3:position]
     else:
         promoter.at[index, '9'] = infor[3:]

promoter.drop(labels = ['2', '3', '6', '8'], inplace =True, axis = 1)
promoter.rename(columns={'9':'ID', '4':'Start position', '5':'End position', '7':'strand', '1':'Chroo'}, inplace = True)
promoter.loc[promoter['Chroo'] == "I", 'Chr'] = 1
promoter.loc[promoter['Chroo'] == "II", 'Chr'] = 2
promoter.loc[promoter['Chroo'] == "III", 'Chr'] = 3
promoter.loc[promoter['strand'] == "-", 'Strand'] = 'reverse'
promoter.loc[promoter['strand'] == "+", 'Strand'] = 'forward'



for line in utr5.itertuples():
     index = line[0]
     infor = line[9]
     position = infor.rfind(';')
    # print(position)
     if position >0:
         utr5.at[index, '9'] = infor[3:position]
     else:
         utr5.at[index, '9'] = infor[3:]

utr5.drop(labels = ['2', '3', '6', '8'], inplace =True, axis = 1)
utr5.rename(columns={'9':'ID', '4':'Start position', '5':'End position', '7':'strand', '1':'Chroo'}, inplace = True)
utr5.loc[utr5['Chroo'] == "I", 'Chr'] = 1
utr5.loc[utr5['Chroo'] == "II", 'Chr'] = 2
utr5.loc[utr5['Chroo'] == "III", 'Chr'] = 3
utr5.loc[utr5['strand'] == "-", 'Strand'] = 'reverse'
utr5.loc[utr5['strand'] == "+", 'Strand'] = 'forward'

###############

for line in ncRNA.itertuples():
     index = line[0]
     infor = line[9]
     position = infor.rfind(';')
    # print(position)
     if position >0:
         ncRNA.at[index, '9'] = infor[3:position]
     else:
         ncRNA.at[index, '9'] = infor[3:]

ncRNA.drop(labels = ['2', '3', '6', '8'], inplace =True, axis = 1)
ncRNA.rename(columns={'9':'ID', '4':'Start position', '5':'End position', '7':'strand', '1':'Chroo'}, inplace = True)
ncRNA.loc[ncRNA['Chroo'] == "I", 'Chr'] = 1
ncRNA.loc[ncRNA['Chroo'] == "II", 'Chr'] = 2
ncRNA.loc[ncRNA['Chroo'] == "III", 'Chr'] = 3
ncRNA.loc[ncRNA['strand'] == "-", 'Strand'] = 'reverse'
ncRNA.loc[ncRNA['strand'] == "+", 'Strand'] = 'forward'



for line in snRNA.itertuples():
     index = line[0]
     infor = line[9]
     position = infor.rfind(';')
    # print(position)
     if position >0:
         snRNA.at[index, '9'] = infor[3:position]
     else:
         snRNA.at[index, '9'] = infor[3:]

snRNA.drop(labels = ['2', '3', '6', '8'], inplace =True, axis = 1)
snRNA.rename(columns={'9':'ID', '4':'Start position', '5':'End position', '7':'strand', '1':'Chroo'}, inplace = True)
snRNA.loc[snRNA['Chroo'] == "I", 'Chr'] = 1
snRNA.loc[snRNA['Chroo'] == "II", 'Chr'] = 2
snRNA.loc[snRNA['Chroo'] == "III", 'Chr'] = 3
snRNA.loc[snRNA['strand'] == "-", 'Strand'] = 'reverse'
snRNA.loc[snRNA['strand'] == "+", 'Strand'] = 'forward'


for line in snoRNA.itertuples():
     index = line[0]
     infor = line[9]
     position = infor.rfind(';')
    # print(position)
     if position >0:
         snoRNA.at[index, '9'] = infor[3:position]
     else:
         snoRNA.at[index, '9'] = infor[3:]

snoRNA.drop(labels = ['2', '3', '6', '8'], inplace =True, axis = 1)
snoRNA.rename(columns={'9':'ID', '4':'Start position', '5':'End position', '7':'strand', '1':'Chroo'}, inplace = True)
snoRNA.loc[snoRNA['Chroo'] == "I", 'Chr'] = 1
snoRNA.loc[snoRNA['Chroo'] == "II", 'Chr'] = 2
snoRNA.loc[snoRNA['Chroo'] == "III", 'Chr'] = 3
snoRNA.loc[snoRNA['strand'] == "-", 'Strand'] = 'reverse'
snoRNA.loc[snoRNA['strand'] == "+", 'Strand'] = 'forward'

for line in rRNA.itertuples():
     index = line[0]
     infor = line[9]
     position = infor.rfind(';')
    # print(position)
     if position >0:
         rRNA.at[index, '9'] = infor[3:position]
     else:
         rRNA.at[index, '9'] = infor[3:]

rRNA.drop(labels = ['2', '3', '6', '8'], inplace =True, axis = 1)
rRNA.rename(columns={'9':'ID', '4':'Start position', '5':'End position', '7':'strand', '1':'Chroo'}, inplace = True)
rRNA.loc[rRNA['Chroo'] == "I", 'Chr'] = 1
rRNA.loc[rRNA['Chroo'] == "II", 'Chr'] = 2
rRNA.loc[rRNA['Chroo'] == "III", 'Chr'] = 3
rRNA.loc[rRNA['strand'] == "-", 'Strand'] = 'reverse'
rRNA.loc[rRNA['strand'] == "+", 'Strand'] = 'forward'

###############


for line in ori.itertuples():
     index = line[0]
     infor = line[9]
     position = infor.rfind(';')
    # print(position)
     if position >0:
         ori.at[index, '9'] = infor[3:position]
     else:
         ori.at[index, '9'] = infor[3:]

ori.drop(labels = ['2', '3', '6', '8'], inplace =True, axis = 1)
ori.rename(columns={'9':'ID', '4':'Start position', '5':'End position', '7':'strand', '1':'Chroo'}, inplace = True)
ori.loc[ori['Chroo'] == "I", 'Chr'] = 1
ori.loc[ori['Chroo'] == "II", 'Chr'] = 2
ori.loc[ori['Chroo'] == "III", 'Chr'] = 3
ori.loc[ori['strand'] == "-", 'Strand'] = 'reverse'
ori.loc[ori['strand'] == "+", 'Strand'] = 'forward'


centro_merged.loc[centro_merged['chro'] == "chr1", 'Chr'] = 1
centro_merged.loc[centro_merged['chro'] == "chr2", 'Chr'] = 2
centro_merged.loc[centro_merged['chro'] == "chr3", 'Chr'] = 3

df_overlap_tRNAs = pd.DataFrame()
df_within_centromeres = pd.DataFrame()

tRNA_grouped = tRNA.groupby('Chr')
centromere_grouped = centro_merged.groupby('Chr')

# Step 2: Iterate through each chromosome group to identify tRNAs within or overlapping centromeres
overlap_tRNAs = []
within_centromeres = []

for chro, tRNA_chro_group in tRNA_grouped:
    centromere_chro_group = centromere_grouped.get_group(chro) if chro in centromere_grouped.groups else None

    if centromere_chro_group is not None:
        for _, tRNA_row in tRNA_chro_group.iterrows():
            start_tRNA, end_tRNA = tRNA_row['Start position'], tRNA_row['End position']
            overlapping = False
            within_centromere = False

            for _, centromere_row in centromere_chro_group.iterrows():
                start_centromere, end_centromere = centromere_row['start'], centromere_row['end']

                # Check if tRNA overlaps with centromere
                if start_tRNA <= end_centromere and end_tRNA >= start_centromere:
                    overlapping = True

                    # Check if tRNA is entirely within centromere
                    if start_tRNA >= start_centromere and end_tRNA <= end_centromere:
                        within_centromere = True
                        break

            if overlapping:
                if within_centromere:
                    within_centromeres.append(tRNA_row)
                else:
                    overlap_tRNAs.append(tRNA_row)

# Step 3: Create two separate DataFrames for tRNAs that overlap with centromeres and those within centromeres
df_overlap_tRNAs = pd.DataFrame(overlap_tRNAs)
df_within_centromeres = pd.DataFrame(within_centromeres)

# Print the two DataFrames
print("tRNAs that overlap with centromeres:")
print(df_overlap_tRNAs)

print("\ntRNAs that are within centromeres:")
print(df_within_centromeres)

df_tRNA


xls = pd.ExcelFile('S_pombe_expression_data_PMID_23101633.xlsx')
df1 = pd.read_excel(xls, 'Table_S4')
df1.columns=df1.iloc[5]

df1 = df1.iloc[6: , :]

rob = df1[["Systematic.name", "MM.mRNA.cpc"]]


rob.rename({'Systematic.name': 'Systematic ID'}, axis=1, inplace=True)
rob.rename({"MM.mRNA.cpc": 'rna'}, axis=1, inplace=True)

joined = rob.merge(ffeat, how='left', on='Systematic ID')
joinedd = joined.dropna(how='any', axis = 0)





#featfor = joinedd[(joinedd['Strand'] == 'forward')]
#featrev = joinedd[(joinedd['Strand'] == 'reverse')]

#joinedd1 = joinedd[(joinedd['Chr'] == 1)]
#joinedd2 = joinedd[(joinedd['Chr'] == 2)]
#joinedd3 = joinedd[(joinedd['Chr'] == 3)]



GenesandRob = rob.merge(Genes, how='left', on='Systematic ID')
GenesandRoby = GenesandRob.dropna(how='any', axis = 0)
joinedd1 = GenesandRoby[(GenesandRoby['Chr'] == 1)]
joinedd2 = GenesandRoby[(GenesandRoby['Chr'] == 2)]
joinedd3 = GenesandRoby[(GenesandRoby['Chr'] == 3)]




featfor = GenesandRob[(GenesandRob['Strand'] == 'forward')]
featrev = GenesandRob[(GenesandRob['Strand'] == 'reverse')]


def score(file,frame1, frame2):
    for i in range(len(file)):
        tempstart = file.iloc[i]["Start position"]
        tempend = file.iloc[i]['End position']
        tempchro = file.iloc[i]["Chr"]
        rna = file.iloc[i]['rna']

        tempsubset = frame1.loc[(frame1['Pos']>= tempstart) & (frame1['Pos']<=tempend) & (frame1['Chr'] == tempchro)]
        tempsubsett = frame2.loc[(frame2['Pos']>= tempstart) & (frame2['Pos']<=tempend) & (frame2['Chr'] == tempchro)]
       # tempsubsetm = frame3.loc[(frame1['Pos']>= tempstart) & (frame3['Pos']<=tempend) & (frame3['Chr'] == tempchro)]
        file.loc[file.index[i], 'thing'] = tempsubset['normscore'].sum()
        file.loc[file.index[i], 'epb'] = tempsubset['normscore'].sum()/len(tempsubset)
        
        file.loc[file.index[i], 'thingt'] = tempsubsett['normscore'].sum()
        file.loc[file.index[i], 'epbt'] = tempsubsett['normscore'].sum()/len(tempsubsett)
        
        
        file.loc[file.index[i], 'log2_epb'] = np.log2(file.loc[file.index[i], 'epb'])
        file.loc[file.index[i], 'log2_epbt'] = np.log2(file.loc[file.index[i], 'epbt'])

        # Calculate the length of the region and add log2_length column
        length = tempend - tempstart
        file.loc[file.index[i], 'log2_length'] = np.log2(length)
        file.loc[file.index[i], 'log2_rna'] = np.log2(rna)
        

    return file

featfor = score(featfor, wtET, t2E)
  

def scorerev(file,frame, frame2):
    for i in range(len(file)):
        tempstart = file.iloc[i]["Start position"]
        tempend = file.iloc[i]['End position']
        tempchro = file.iloc[i]["Chr"]
        rna = file.iloc[i]['rna']

        tempsubset = frame.loc[(frame['Pos']>= tempend) & (frame['Pos']<=tempstart) & (frame['Chr'] == tempchro)]
        tempsubsett = frame2.loc[(frame2['Pos']>= tempend) & (frame2['Pos']<=tempstart) & (frame2['Chr'] == tempchro)]

        file.loc[file.index[i], 'thing'] = tempsubset['normscore'].sum()
        file.loc[file.index[i], 'epb'] = tempsubset['normscore'].sum()/len(tempsubset)
        
        file.loc[file.index[i], 'thingt'] = tempsubsett['normscore'].sum()
        file.loc[file.index[i], 'epbt'] = tempsubsett['normscore'].sum()/len(tempsubsett)
        
        file.loc[file.index[i], 'log2_epb'] = np.log2(file.loc[file.index[i], 'epb'])
        file.loc[file.index[i], 'log2_epbt'] = np.log2(file.loc[file.index[i], 'epbt'])

        # Calculate the length of the region and add log2_length column
        length = tempend - tempstart
        file.loc[file.index[i], 'log2_length'] = np.log2(length)
        file.loc[file.index[i], 'log2_rna'] = np.log2(rna)

    return file
    
featrev = scorerev(featrev, wtET, t2E) 


merged_feat = pd.concat([featfor, featrev], axis=0)

merged_feat['length_quartile'] = pd.qcut(merged_feat['log2_length'], q=4, labels=False)
merged_feat['epb_quartile'] = pd.qcut(merged_feat['log2_epb'], q=4, labels=False)
merged_feat['expression_quartile'] = pd.qcut(merged_feat['rna'], q=4, labels=False)

    merged_feat.loc[merged_feat['expression_quartile'] == 0, 'expression_quartile'] = 'Q1'
    merged_feat.loc[merged_feat['expression_quartile'] == 1, 'expression_quartile'] = 'Q2'
    merged_feat.loc[merged_feat['expression_quartile'] == 2, 'expression_quartile'] = 'Q3'
    merged_feat.loc[merged_feat['expression_quartile'] == 3, 'expression_quartile'] = 'Q4'


merged_featy = merged_feat.copy()
merged_featy["expression_quartile"] = pd.Categorical(merged_featy["expression_quartile"],
                                                 ordered = True,
                                                 categories = ["Q1", "Q2", "Q3", "Q4"])

print(merged_featy['expression_quartile'])
# Custom color palette
custom_palette = ["#1D5B79", "#468B97", "#EF6262", "#F3AA60"]  # Replace these colors with your desired ones

# Set the custom color palette
sns.set_palette(custom_palette)

custom_palette = ["#1D5B79", "#EF6262" ]  # Replace these colors with your desired ones

# Set the custom color palette
sns.set_palette(custom_palette)

import seaborn as sns
sns.jointplot(
    data=merged_featy,
    x="log2_rna", y="log2_epb", hue="expression_quartile"
    
)

sns.jointplot(
    data=merged_featy,
    x="log2_rna", y="log2_epb", hue="expression_quartile",
    ylim= [-8,-5]
)

sns.lmplot(data=merged_featy, x="log2_rna", y="log2_epb", hue="expression_quartile")




sns.lmplot(
    data=merged_featy, x="log2_rna", y="log2_epb",
    col="expression_quartile", height=3,
    facet_kws=dict(sharex=False, sharey=False),
)





dataIGRy = dataIGR.copy()
dataIGRy["Quartile"] = pd.Categorical(dataIGRy["Quartile"],
                                                 ordered = True,
                                                 categories = ["Q1", "Q2", "Q3", "Q4"])
print(dataIGRy['Quartile'])

sns.jointplot(
    data=dataIGRy,
    x="log2_length", y="log2_epb", hue="Quartile",
    kind = 'kde', ylim= [-1,-10]
)

#%%
#I'm going to try and write my igr funnction mysellf using joinedd

def last_resort (file):
    #data_output = pd.DataFrame()
    igr_start = []
    igr_end = []
    rnalist = []

    igr_types = [] 
    name1 = []
    name2 = []
    name3 = []
    
    for i in range(len(file)-2):
        tempname = file.iloc[i]['Systematic ID']
        tempstart = file.iloc[i]["Start position"]
        tempend = file.iloc[i]["End position"]
        tempstrand = file.iloc[i]['Strand']
        
        nextname = file.iloc[i+1]['Systematic ID']
        nextstart = file.iloc[i+1]["Start position"]
        nextend = file.iloc[i+1]["End position"]
        nextstrand = file.iloc[i+1]['Strand']
        rnale = file.iloc[i+1]['rna']


        lastname = file.iloc[i+2]['Systematic ID']
        laststart = file.iloc[i+2]["Start position"]
        lastend = file.iloc[i+2]["End position"]
        laststrand = file.iloc[i+2]['Strand']
        
        if nextstart - tempend < 0:
            #this is an overlap

            name2.append(nextname)
            rnalist.append(rnale)
            igr_types.append('overlap')
            igr_start.append(tempend)
            igr_end.append(nextstart)
            
                    
        elif nextstart > tempend and tempstrand == nextstrand:
            #this is a tandem igr 
            tandem_len = nextstart - tempend
            rnalist.append(rnale)
            name2.append(nextname)
            igr_types.append('tandem')
            igr_start.append(tempend)
            igr_end.append(nextstart)
        
           # print(tandem_len)
            #asign this to a new df. 

            
        elif nextstart > tempend and tempstrand == 'reverse' and nextstrand == 'forward':
                #this is a divergent igr
            name2.append(nextname)
            rnalist.append(rnale)
            divergentlen = nextstart - tempend
            igr_types.append('divergent')
            igr_start.append(tempend)
            igr_end.append(nextstart)

            
            
            
               # print(divergentlen)
                
        elif nextstart > tempend and tempstrand == 'forward' and nextstrand == 'reverse':
            #this is a convergent igr
            name2.append(nextname)
            rnalist.append(rnale)
            convergentlen = nextstart - tempend
            igr_types.append('convergent')
            igr_start.append(tempend)
            igr_end.append(nextstart)

            
            
            
    print(len(igr_start))
    print(len(igr_end))
    print(len(igr_types))
    print(len(nextname))
    print(len(rnalist))
    igr_df = pd.DataFrame({
        'Start position': igr_start,
        'End position': igr_end,
        'IGR Type': igr_types,
        'IGR name': name2,
        'rna': rnalist
        
        
        # Add the IGR Type column
    })
        
    return igr_df

dye = last_resort(joinedd1)
dye2 = last_resort(joinedd2)     
dye3 = last_resort(joinedd3)     
        
dye['Chr'] = 1
dye2['Chr'] = 2
dye3['Chr'] = 3

megeed_dye = pd.concat([dye, dye2, dye3], ignore_index=True)

def check_plot (centro, genee, telo, frig):
    ff, (ax6) = plt.subplots(1,1, sharex=True)

    ax6.set_ylabel('Gene Annotations')
    ax6.set_xlabel('Chromosome position')

 
    for c in centro.itertuples(index=False, name=None):
            ax6.axvspan(c[0],c[1],0.1,0.15,color="black",alpha=0.8)
            
    for t in telo.itertuples(index=False, name=None):
            ax6.axvspan(t[0],t[1],0.1,0.15,color="grey",alpha=0.8)
            
    for ge in genee.itertuples(index=False, name=None):
        
        if ge[4] == 'protein coding':
            if ge[8] == 'reverse':
                ax6.axvspan(ge[5],ge[6],0.2,0.4,color="firebrick",alpha=0.3)
                #ax6.annotate(ge[0], xy = [ge[4],0.45]) 
            elif ge[8] == 'forward':
                ax6.axvspan(ge[5],ge[6],0.7,0.9,color="steelblue",alpha=0.3)
                #ax6.annotate(ge[0], xy = [ge[4],0.65]) 
                

    
    for fr in frig.itertuples(index=False, name=None):
        if fr[2] == 'tandem':
            ax6.axvspan(fr[0],fr[1],0.4,0.6,color="yellow",alpha=0.3)
        if fr[2] == 'convergent':
            ax6.axvspan(fr[0],fr[1],0.4,0.6,color="red",alpha=0.3)
        if fr[2] == 'divergent':
            ax6.axvspan(fr[0],fr[1],0.4,0.6,color="blue",alpha=0.3)
        if fr[2] == 'overlap':
            ax6.axvspan(fr[0],fr[1],0.4,0.6,color="black",alpha=0.3)

    return ff

chromosome1 = check_plot(centro1, joinedd1, telo1, dye)



#%%
def Create_df(df,dr,ef,er,w):
    df = pd.read_csv(df) # Reads a csv file
    df['count'].replace(0,1, inplace = True) # This replaces zeros with ones. inplace = true saves the value to the dataframe (not useing a copy)
    df.rename(columns = {"count" : "df_count"}, inplace = True) # renames the colunm to appropriate counts
    
    dr = pd.read_csv(dr, usecols=[2]) # Read only the counts column from the next file
    dr['count'].replace(0,1, inplace = True)
    dr.rename(columns = {'count' : 'dr_count'}, inplace = True) # renames the colunm to appropriate counts 
    
    ef = pd.read_csv(ef, usecols=[2])
    ef['count'].replace(0,1, inplace = True)   
    ef.rename(columns = {'count' : 'ef_count'}, inplace = True)
    
    er = pd.read_csv(er, usecols=[2])
    er['count'].replace(0,1, inplace = True)
    er.rename(columns = {'count' : 'er_count'}, inplace = True)
    
    all_data = pd.concat([df, dr, ef, er], axis=1, join='outer') # Create a single dataframe by merging the 4 dataframes

  
# Here we add new colums to the dataframe
# First we normalise the data. ie each line in the column is simply the value divided by the sum of the column
# Note: pandas automatically itterate through the rows.    
    all_data['norm_df'] = all_data['df_count']/all_data['df_count'].sum()
    all_data['norm_dr'] = all_data['dr_count']/all_data['dr_count'].sum()
    all_data['norm_ef'] = all_data['ef_count']/all_data['ef_count'].sum()
    all_data['norm_er'] = all_data['er_count']/all_data['er_count'].sum()

# Next we calculate the ratios for each strand and assign to a new colunm
    all_data['ratio_delta_f'] = all_data['norm_df']/(all_data['norm_df'] + all_data['norm_ef'])
    all_data['ratio_delta_r'] = all_data['norm_dr']/(all_data['norm_dr'] + all_data['norm_er'])
    all_data['ratio_epsilon_f'] = all_data['norm_ef']/(all_data['norm_ef'] + all_data['norm_df'])
    all_data['ratio_epsilon_r'] = all_data['norm_er']/(all_data['norm_er'] + all_data['norm_dr'])


# Now we have column for pol delta useage for the duplex
    all_data['d_usage'] = (all_data['ratio_delta_f'] + all_data['ratio_delta_r']) / 2

# now we a column for the percentage of right-moving forks
    all_data['right_forks']  = all_data['ratio_epsilon_f']*2 - 1
    
    
# now we will produce a new colum for each sliding window average for each of the calculated columns
# Note: centre = True, means that the data is summed from both left and right. False means its the sum of the last of the number of values.
    
    all_data['smoo_ratio_d_f'] = all_data['ratio_delta_f'].rolling(window = w, center=True).mean()
    all_data['smoo_ratio_d_r'] = all_data['ratio_delta_r'].rolling(window = w, center=True).mean()
    all_data['smoo_ratio_e_f'] = all_data['ratio_epsilon_f'].rolling(window = w, center=True).mean()
    all_data['smoo_ratio_e_r'] = all_data['ratio_epsilon_r'].rolling(window = w, center=True).mean()
    all_data['smoo_d_usage'] = all_data['d_usage'].rolling(window = w, center=True).mean()
    all_data['smoo_right_forks'] = all_data['right_forks'].rolling(window = w, center=True).mean()
    

# now create a differential collunm and a discreate orgin (centre = 1 bin) column.
# Note the script removes differentials not present in both strands and which are singular 
    
    # take unsmoothed pol usage data into two pandas arrays
    ser_etop = all_data['ratio_epsilon_f']
    ser_ebot = all_data['ratio_epsilon_r']

    # Calculate the differentials 
    ser_topd = ser_etop.diff()
    ser_botd = ser_ebot.diff()
    # Reverse the sign on the bottom strand differentials to match top strand
    ser_botd = ser_botd*-1

    # curently dont use a roling average, but here they are if needed
    #ser_topd = ser_topd.rolling(window = 3, center=True).mean()
    #ser_botd = ser_botd.rolling(window = 3, center=True).mean()

    # Removing all the negative valuse to zero
    ser_topd.clip(0, 1, inplace=True)
    ser_botd.clip(0, 1, inplace=True)
    
    # If the value in one is zero, then seting it in zero in both datasets (this removes a lot of noise)
    for i in range(len(ser_topd)):
        if ser_topd.iloc[i] == 0:
            ser_botd.iloc[i] = 0
    for i in range(len(ser_botd)):
        if ser_botd.iloc[i] == 0:
            ser_topd.iloc[i] = 0
    
    # Now we add the two things together and divide by two - i.e we have an average origin activity
    ser_cumd = (ser_topd + ser_botd)/2     

    # Now we want to calculate the quantile of all the positive data.
    templist = np.array([])
    for i in range(1,len(ser_cumd)):
        if ser_cumd.iloc[i] != 0: 
            templist = np.append(templist, ser_cumd.iloc[i])
    # set a cutoff threshold at of the top 10% of all positive values
    cutoff = np.quantile(templist, 0.9)

    # now if a single +ve value (i.e at least one zero either side) is below this cutoff, then set to zero. This again removes noise.
    for i in range(1, len(ser_cumd)-1):
        if ser_cumd.iloc[i] != 0:
            if ser_cumd.iloc[i-1] == 0 and ser_cumd.iloc[i+1] == 0 and ser_cumd.iloc[i]< cutoff:
                ser_cumd.iloc[i] = 0

    # Now we save the data back to the daframe. This is a list of differentials (origin activity) per bin.
    all_data['differentials'] = ser_cumd

    # Next we want to position "zonal" origins to a single point (bin) and give them a cumulative value (efficiency of the zone)
    start=0
    stop=0
    origins = ser_cumd.clip(0, 0) # creates a new pd.series of zero values same length of ser_cumd
    for i in range(len(ser_cumd)):
        if i <= stop: # simply prevents itterating over the non-zero block
            continue # continue goes back to the for loop
        else:
            if ser_cumd.iloc[i] != 0:
                start = i        
                for z in range(start+1, len(ser_cumd)): # find the ned of the non-zero block
                    if ser_cumd.iloc[z] == 0:
                        stop = z
                        tot = 0
                        tem = 0
                        tem_loc = 0
                        for v in range(i,z): # adds the non-zero block numbers together and identifies the location of the highest value
                            tot = tot + ser_cumd.iloc[v]
                            if ser_cumd.iloc[v] > tem:
                                tem = ser_cumd.iloc[v]
                                tem_loc = v
                        origins.iloc[tem_loc] = tot # adds the total of the non-zero block to the bin with the highest individula differential
                        break

    all_data['origin'] = origins

# create three separate dataframes for the three chromosomes and return these too
 
    chrI = all_data.loc[all_data['chro']=='chr1']
    chrII = all_data.loc[all_data['chro']=='chr2']
    chrIII = all_data.loc[all_data['chro']=='chr3']
    return chrI, chrII, chrIII

 
owtchrI, owtchrII, owtchrIII = Create_df('RZ259-RhArepr-d.e1.f-w300.count.csv', 'RZ259-RhArepr-d.e1.r-w300.count.csv', 'RZ267-RhArepr-e.e1.f-w300.count.csv', 'RZ267-RhArepr-e.e1.r-w300.count.csv', 5)


owtpu = pd.concat([owtchrI, owtchrII, owtchrIII], ignore_index=True)
owtpu['log2_origin'] = np.log2(owtpu['origin'])


owtpufilt = owtpu[owtpu['origin'].notna() & (owtpu['origin'] > 0)]
owtpufilt['origin efficiency'] = pd.qcut(owtpufilt['log2_origin'], q=4, labels=False)

    
seventy_fifth_percentile = np.percentile(owtpufilt['origin'], 75)

# Filter the DataFrame to include values above the 75th percentile
owtpufilt_above_75th = owtpufilt[owtpufilt['origin'] > seventy_fifth_percentile]
owtpufilt_above_75th['high efficiency origins'] = pd.qcut(owtpufilt_above_75th['log2_origin'], q=4, labels=False)


plt.figure(figsize=(10, 4))
sns.kdeplot(
   data=owtpufilt, x="log2_origin", hue='origin efficiency',
   fill=True, common_norm=False,
   alpha=.5, linewidth=1,
)
sns.kdeplot(
   data=owtpufilt, x="log2_origin",
   fill=True, common_norm=False,
   alpha=.5, linewidth=1,
)
sns.kdeplot(
   data=owtpufilt_above_75th, x="log2_origin",
   fill=True, common_norm=False,hue = 'high efficiency origins',
   alpha=.5, linewidth=1,
)

owtpufilt_above_75th.loc[owtpufilt_above_75th['chro'] == "chr1", 'Chr'] = 1
owtpufilt_above_75th.loc[owtpufilt_above_75th['chro'] == "chr2", 'Chr'] = 2
owtpufilt_above_75th.loc[owtpufilt_above_75th['chro'] == "chr3", 'Chr'] = 3



merged_dforiigr = pd.merge(owtpufilt_above_75th, dataIGR
                           , on='Chr')

# Subset the DataFrame to select rows where 'pos' falls between 'Start position' and 'Stop position'
oinigr = merged_dforiigr[(merged_dforiigr['pos'] >= merged_dforiigr['Start position']) & (merged_dforiigr['pos'] <= merged_dforiigr['End position'])]

# Print the 'oinigr' DataFrame
print(oinigr)


other_rows = merged_dforiigr[~((merged_dforiigr['pos'] >= merged_dforiigr['Start position']) & (merged_dforiigr['pos'] <= merged_dforiigr['End position']))]

# Print the 'oinigr' DataFrame
print(oinigr)

# Print the 'other_rows' DataFrame
print(other_rows)


melted_dfoinigr = pd.melt(oinigr, id_vars=['IGR Type'], value_vars=['epb', 'epbt'], var_name='source', value_name='score')
melted_df['log2_score'] = np.log2(melted_df['score'])
# Print the new DataFrame
print(melted_df)


sns.catplot(data=oinigr, x="IGR Type", y="log2_epb", kind="box", hue= 'origin')
oinigr['IGR Type'] = 'IGR containing origin'
combined_df = pd.concat([oinigr, megeed_dye], ignore_index=True)

# Drop duplicates based on all columns (including 'IGR Type')
combined_df = combined_df.drop_duplicates()


megeed_dyeb = megeed_dye.copy()
megeed_dyeb['IGR Type'] = 'total IGR'

# Print the combined DataFrame
print(combined_df)
dfyyy = combined_df[~combined_df.duplicated(subset='Start position', keep=False)]

sns.catplot(data=combined_df, x="IGR Type", y="log2_epb", kind="box")

melted_dfoinigrcombined_df = pd.melt(combined_df, id_vars=['IGR Type'], value_vars=['epb', 'epbt'], var_name='source', value_name='score')
melted_dfoinigrcombined_df['log2_score'] = np.log2(melted_dfoinigrcombined_df['score'])


sns.catplot(data=melted_dfoinigrcombined_df, x="IGR Type", y="log2_score", kind="box", hue= 'source', showfliers=False)

combined_df2 = pd.concat([combined_df, megeed_dyeb], ignore_index=True)




melted_dfoinigrcombined_df2 = pd.melt(combined_df2, id_vars=['IGR Type'], value_vars=['epb', 'epbt'], var_name='source', value_name='score')
melted_dfoinigrcombined_df2['log2_score'] = np.log2(melted_dfoinigrcombined_df2['score'])



sns.catplot(data=melted_dfoinigrcombined_df2, x="IGR Type", y="log2_score", kind="box", hue= 'source', showfliers=False, legend=False)



#now to combine meged_feat with IGR stuff

df_cleanedg['IGR Type'] = 'intragenic'
combiny = pd.concat([df_cleanedg, combined_df2], ignore_index=True)
intrainterr = pd.melt(combiny, id_vars=['IGR Type'], value_vars=['epb', 'epbt'], var_name='source', value_name='score')

intrainterr['log2_score'] = np.log2(intrainterr['score'])

intrainterr["IGR Type"] = pd.Categorical(intrainterr["IGR Type"],
                                                 ordered = True,
                                                 categories = ["intragenic", "total IGR", "IGR containing origin", "divergent", "tandem", 'convergent', "overlap"])

unique_items_array = intrainterr['IGR Type'].unique()
epbonlyintrainterr = intrainterr[intrainterr['source'] != 'epbt']


plt.figure(figsize=(10, 4))
sns.set_palette("Set2")
sns.kdeplot(
   data=epbonlyintrainterr, x="log2_score", hue="IGR Type",
   fill=True, common_norm=False,
   alpha=.5, linewidth=1,
)

wintragenic = intrainterr[(intrainterr['IGR Type'] == 'intragenic')]
wtotaligr = intrainterr[(intrainterr['IGR Type'] == 'total IGR')]
wigrco = intrainterr[(intrainterr['IGR Type'] == 'IGR containing origin')]
wdivegent = intrainterr[(intrainterr['IGR Type'] == 'divergent')]
wtandem = intrainterr[(intrainterr['IGR Type'] == 'tandem')]
wconveegent = intrainterr[(intrainterr['IGR Type'] == 'convergent')]

#now need them to be seperated into epb vs epbt

wintragenicw = wintragenic[(wintragenic['source'] == 'epb')]
wintragenict = wintragenic[(wintragenic['source'] == 'epbt')]

wtotaligrw = wtotaligr[(wtotaligr['source'] == 'epb')]
wtotaligrt = wtotaligr[(wtotaligr['source'] == 'epbt')]

wigrcow = wigrco[(wigrco['source'] == 'epb')]
wigrcot = wigrco[(wigrco['source'] == 'epbt')]

wdivegentw = wdivegent[(wdivegent['source'] == 'epb')]
wdivegentt = wdivegent[(wdivegent['source'] == 'epbt')]

wtandemw = wtandem[(wtandem['source'] == 'epb')]
wtandemt = wtandem[(wtandem['source'] == 'epbt')]

wconveegentw = wtandem[(wtandem['source'] == 'epb')]
wconveegentt = wtandem[(wtandem['source'] == 'epbt')]

from scipy.stats import wilcoxon

# List of pairs of DataFrames
pairs = [
    (wintragenicw, wintragenict),
    (wtotaligrw, wtotaligrt),
    (wigrcow, wigrcot),
    (wdivegentw, wdivegentt),
    (wtandemw, wtandemt),
    (wconveegentw, wconveegentt)
]

# Perform Wilcoxon signed-rank test for each pair
results = []
for df1, df2 in pairs:
    
    
    statistic, p_value = wilcoxon(df1['log2_score'], df2['log2_score'])
    results.append({'DataFrame1': df1.name, 'DataFrame2': df2.name, 'Statistic': statistic, 'p-value': p_value})

# Create a DataFrame to store the results
results_dfwicox = pd.DataFrame(results)

# Display the results
print(results_df)


# Define a list of pairs and their corresponding names
pairs_with_names = [
    (wintragenicw, wintragenict, 'wintragenicw', 'wintragenict'),
    (wtotaligrw, wtotaligrt, 'wtotaligrw', 'wtotaligrt'),
    (wigrcow, wigrcot, 'wigrcow', 'wigrcot'),
    (wdivegentw, wdivegentt, 'wdivegentw', 'wdivegentt'),
    (wtandemw, wtandemt, 'wtandemw', 'wtandemt'),
    (wconveegentw, wconveegentt, 'wconveegentw', 'wconveegentt')
    # Add more pairs in a similar format
]

results = []

for df1, df2, df1_name, df2_name in pairs_with_names:
    # Perform the Wilcoxon test
    statistic, p_value = wilcoxon(df1['log2_score'], df2['log2_score'])
    
    # Append results to the list
    results.append({
        'DataFrame1': df1_name,
        'DataFrame2': df2_name,
        'Statistic': statistic,
        'p-value': p_value
    })

# Create a DataFrame to store the results
results_dfwicox = pd.DataFrame(results)





series_list = [data_q1, data_q2, data_q3, data_q4, data_q5, data_q6]

resultsks = []

for i in range(len(series_list)):
    for j in range(i + 1, len(series_list)):
        ks_statistic, p_value = ks_2samp(series_list[i], series_list[j])
        resultsks.append({
            'Series1': f'Series {i+1}',
            'Series2': f'Series {j+1}',
            'KS Statistic': ks_statistic,
            'p-value': p_value
        })

# Create a DataFrame to store the results
df = pd.DataFrame(resultsks)

# Display the results
print(df)

pivoted_df = df.pivot(index='Series1', columns='Series2', values='p-value')



# Create a heatmap using seaborn
plt.figure(figsize=(10, 8))
sns.heatmap(pivoted_df, annot=True, cmap='RdBu', fmt=".1g", vmin= 1e-107, vmax=1e-11)
plt.title("Correlogram for p-values")
plt.show()


sns.scatterplot(data=df, x="Series1", y="Series2", size="p-value",, legend="full", vmin= 1e-107, vmax=1e-3)

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors

# Create your correlation matrix
correlation_matrix = df.corr()

# Set your vmin and vmax
vmin = 1e-107
vmax = 1e-11

# Define the number of colors for a smoother transition
n_colors = 256

# Create a custom colormap
color_map = plt.get_cmap('coolwarm', n_colors)

# Normalize the values based on vmin and vmax
norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

# Create a scalar mappable with the colormap and normalization
scalar_mappable = plt.cm.ScalarMappable(cmap=color_map, norm=norm)

# Get the colors for the range of values
colors = [scalar_mappable.to_rgba(value) for value in np.linspace(vmin, vmax, n_colors)]

# Create the heatmap
plt.figure(figsize=(10, 8))
sns.heatmap(pivoted_df, annot=True, cmap=colors, fmt=".1g", vmin=vmin, vmax=vmax)
plt.title("Correlogram for p-values")

# Add a color bar with the custom colormap
cbar = plt.colorbar(scalar_mappable)
cbar.set_label("p-value")

plt.show()












data_q1 = wintragenicw['log2_score']
data_q2 = wtotaligrw['log2_score']
data_q3 = wigrcow['log2_score']
data_q4 = wdivegentw['log2_score']
data_q5 = wtandemw['log2_score']
data_q6 = wconveegentw['log2_score']


##KS stats
def ecdf(data):
    x = np.sort(data)
    y = np.arange(1, len(data) + 1) / len(data)
    return x, y

x1, y1 = ecdf(data_q1)
x2, y2 = ecdf(data_q2)
x3, y3 = ecdf(data_q3)
x4, y4 = ecdf(data_q4)
x5, y5 = ecdf(data_q5)
x6, y6 = ecdf(data_q6)

# Perform KS-2samp test (IGR)
ks_statistic12, p_value12 = ks_2samp(data_q1, data_q2)
ks_statistic13, p_value13 = ks_2samp(data_q1, data_q3)
ks_statistic45, p_value45 = ks_2samp(data_q4, data_q5)
ks_statistic5, p_value13 = ks_2samp(data_q1, data_q6)

#for genes
ks_statistic14, p_value14 = ks_2samp(data_q1, data_q4)
ks_statistic24, p_value24 = ks_2samp(data_q2, data_q4)
ks_statistic34, p_value34 = ks_2samp(data_q3, data_q4)





# Create ECDF plot GENE
plt.figure(figsize=(10, 4))

plt.subplot(1, 3, 1)
plt.step(x1, y1, label='Q1', color='#EF6262')
plt.step(x3, y3, label='Q4', color='#F3AA60')
plt.xlabel('log2_epb_gene')
plt.ylabel('ECDF')
plt.title('Q1 vs. Q3')
plt.annotate("KS = {:.3f}".format(ks_statistic13), (-6,0.07))
#plt.annotate("p-value = {:.5f}".format(p_value34), (-6,0.03))
plt.legend()

plt.subplot(1, 3, 2)
plt.step(x2, y2, label='Q2', color='#468B97')
plt.step(x4, y4, label='Q4', color='#F3AA60')
plt.annotate("KS = {:.3f}".format(ks_statistic24), (-6,0.07))
#plt.annotate("p-value = {:.9f}".format(p_value24), (-7,0.03))
plt.xlabel('log2_epb_gene')
plt.ylabel('ECDF')
plt.title('Q2 vs. Q4')
plt.legend()

plt.subplot(1, 3, 3)
plt.step(x1, y1, label='Q1', color='#1D5B79')
plt.step(x4, y4, label='Q4', color='#F3AA60')
plt.annotate("KS = {:.3f}".format(ks_statistic14), (-6,0.07))
#plt.annotate("p-value = {:.14f}".format(p_value14), (-9,0.03))
plt.xlabel('log2_epb_gene')
plt.ylabel('ECDF')
plt.title('Q1 vs. Q4')
plt.legend()

##ECDF PLOT IGR

ks_statistic23, p_value23 = ks_2samp(data_q2, data_q3)
ks_statistic13, p_value13 = ks_2samp(data_q1, data_q3)




sns.catplot(data=intrainterr, x="IGR Type", y="log2_score", kind="box", hue= 'source', showfliers=False)

dataIGR.rename(columns={'IGR name': 'Systematic ID'}, inplace=True)

joined_df_of_friends = GenesandRoby.join(dataIGR.set_index('Systematic ID'), on='Systematic ID')

joined_df_of_friends = pd.merge(merged_featy, dataIGR, on='Systematic ID', suffixes=('_gene', '_IGR'))
joined_df_of_friends['rna'] = joined_df_of_friends['rna'].astype(int)
joined_df_of_friends['log2_rna'] = np.log2(joined_df_of_friends['rna'])


sns.jointplot(
    data=joined_df_of_friends,
    x="log2_rna", y="log2_epb_gene", hue="expression_quartile", kind= 'kde'
    
)



#######
#####
#####
df_cleaned = dataIGR.dropna(subset=['log2_rna', 'log2_epb'])

df_cleaned




df_cleaned['expression_quartile'] = pd.qcut(df_cleaned['rna'], q=4, labels=False)

df_cleaned.loc[df_cleaned['expression_quartile'] == 0, 'expression_quartile'] = 'Q1'
df_cleaned.loc[df_cleaned['expression_quartile'] == 1, 'expression_quartile'] = 'Q2'
df_cleaned.loc[df_cleaned['expression_quartile'] == 2, 'expression_quartile'] = 'Q3'
df_cleaned.loc[df_cleaned['expression_quartile'] == 3, 'expression_quartile'] = 'Q4'


df_cleaned = df_cleaned.copy()
df_cleaned["expression_quartile"] = pd.Categorical(df_cleaned["expression_quartile"],
                                                 ordered = True,
                                                 categories = ["Q1", "Q2", "Q3", "Q4"])

unique_items_array = df_cleaned['IGR Type'].unique()

df_cleaned_tandem = df_cleaned[(df_cleaned['IGR Type'] == 'tandem')]
df_cleaned_convergent = df_cleaned[(df_cleaned['IGR Type'] == 'divergent')]
df_cleaned_divergent = df_cleaned[(df_cleaned['IGR Type'] == 'convergent')]

g = sns.JointGrid()
x, y = df_cleaned_divergent["log2_rna"], df_cleaned_divergent["log2_epb"]
sns.regplot(x=x, y=y, scatter_kws={'s': 1}, line_kws={'color': 'black'}, ax=g.ax_joint)
sns.histplot(data=df_cleaned_divergent, x="log2_rna", hue="expression_quartile",
    element="step",  common_norm = False, ax=g.ax_marg_x, legend = False)
sns.histplot(data=df_cleaned_divergent, y="log2_epb", hue="expression_quartile",
    element="step",  common_norm = False,  ax=g.ax_marg_y, legend = False)



g = sns.JointGrid()
x, y = df_cleaned["log2_rna"], df_cleaned["log2_epb"]
sns.regplot(x=x, y=y, scatter_kws={'s': 1}, line_kws={'color': 'black'}, ax=g.ax_joint)
sns.histplot(data=df_cleaned, x="log2_rna", hue="expression_quartile",
    element="step",  common_norm = False, ax=g.ax_marg_x, legend = False)
sns.histplot(data=df_cleaned, y="log2_epb", hue="expression_quartile",
    element="step",  common_norm = False,  ax=g.ax_marg_y, legend = False)

plt.ylim(-7.5,-4.5)






from sklearn.linear_model import LinearRegression

x = np.array(df_cleaned['log2_rna']).reshape((-1,1))
y = np.array(df_cleaned['log2_epb'])
#model = LinearRegression()
#model.fit(x, y)
model = LinearRegression().fit(x, y)
r_sq = model.score(x, y)
print('coefficient of determination:', r_sq)



#%%
df_cleanedg = merged_feat.dropna(subset=['log2_rna', 'log2_epb'])
df_cleanedg['expression_quartiley'] = pd.qcut(df_cleanedg['rna'], q=4, labels=False)

df_cleanedg.loc[df_cleanedg['expression_quartiley'] == 0, 'expression_quartiley'] = 'Q1'
df_cleanedg.loc[df_cleanedg['expression_quartiley'] == 1, 'expression_quartiley'] = 'Q2'
df_cleanedg.loc[df_cleanedg['expression_quartiley'] == 2, 'expression_quartiley'] = 'Q3'
df_cleanedg.loc[df_cleanedg['expression_quartiley'] == 3, 'expression_quartiley'] = 'Q4'


df_cleanedgy = df_cleanedg.copy()
df_cleanedgy["expression_quartiley"] = pd.Categorical(df_cleanedgy["expression_quartiley"],
                                                 ordered = True,
                                                 categories = ["Q1", "Q2", "Q3", "Q4"])


g = sns.JointGrid()
x, y = df_cleanedgy["log2_rna"], df_cleanedgy["log2_epb"]
sns.regplot(x=x, y=y, scatter_kws={'s': 1}, line_kws={'color': 'black'}, ax=g.ax_joint)
sns.histplot(data=df_cleanedgy, x="log2_rna", hue="expression_quartiley",
    element="step",  common_norm = False, ax=g.ax_marg_x, legend = False)
sns.histplot(data=df_cleanedgy, y="log2_epb", hue="expression_quartiley",
    element="step",  common_norm = False,  ax=g.ax_marg_y, legend = False)

plt.ylim(-7.5,-4.5)



from sklearn.linear_model import LinearRegression



xx = np.array(df_cleanedgy['log2_rna']).reshape((-1, 1))
yy = np.array(df_cleanedgy['log2_epb'])
modell = LinearRegression().fit(xx, yy)
r_sqC = modell.score(xx, yy)
print('coefficient of determination:', r_sqC)



# Assuming you have defined 'xx' and 'yy' as in your previous code
modell = LinearRegression().fit(x, y)

# Get the slope (coefficient) and intercept from the model
slope = modell.coef_[0]
intercept = modell.intercept_

# Print the equation of the regression line
print(f'Equation of the regression line: y = {slope:.4f}x + {intercept:.4f}')


####
import statsmodels.api as sm

# Assuming you have defined 'xx' and 'yy' as in your previous code
xx = np.array(df_cleaned_divergent['log2_rna']).reshape((-1, 1))
yy = np.array(df_cleaned_divergent['log2_epb'])

# Add a constant term to the independent variable
xx = sm.add_constant(xx)

# Create a linear regression model
model = sm.OLS(yy, xx).fit()

# Get the p-values for the coefficients
p_values = model.pvalues

# Print the p-values
print('P-values:', p_values)

# Print the equation of the regression line
slope = model.params[1]
intercept = model.params[0]
print(f'Equation of the regression line: y = {slope:.4f}x + {intercept:.4f}')

#%%


######
jointplot = sns.jointplot(
    data=joined_df_of_friends,
    x="log2_rna", y="log2_epb_gene",
    kind="reg",scatter_kws={'s': 1}, line_kws={'color': 'black'}, ylim=(-7.5,-4.5)
)

# Calculate and display the regression equation
def regression_line(x, y):
    slope, intercept = np.polyfit(x, y, 1)
    return f'y = {slope:.2f}x + {intercept:.2f}'

sns.regplot(data=joined_df_of_friends, x="log2_rna", y="log2_epb_IGR", scatter=False, ax=jointplot.ax_joint)
jointplot.ax_joint.annotate(regression_line(joined_df_of_friends['log2_rna'], joined_df_of_friends['log2_epb_IGR']), 
                            xy=(0.5, 0.1), xycoords='axes fraction', fontsize=12, color='red')

plt.show()



from statsmodels.stats.outliers_influence import variance_inflation_factor

# Load your dataset into a DataFrame (df)
# For each predictor variable, calculate VIF
vif_data = pd.DataFrame()
vif_data["Variable"] = joined_df_of_friends.columns
vif_data["VIF"] = [variance_inflation_factor(joined_df_of_friends.values, i) for i in range(len(joined_df_of_friends.columns))]

print(vif_data)

# Sample data
x = joined_df_of_friends['log2_rna']
y = joined_df_of_friends['log2_epb_IGR']

# Fit a linear regression model
slope, intercept = np.polyfit(x, y, 1)

# Equation of the line: y = mx + b
equation = f'y = {slope:.2f}x + {intercept:.2f}'
print(equation)

#######
from scipy.stats import ks_2samp

joined_df_of_friendsQ1 = joined_df_of_friends[(joined_df_of_friends['expression_quartile'] == 0)]
joined_df_of_friendsQ2 = joined_df_of_friends[(joined_df_of_friends['expression_quartile'] == 1)]
joined_df_of_friendsQ3 = joined_df_of_friends[(joined_df_of_friends['expression_quartile'] == 2)]
joined_df_of_friendsQ4 = joined_df_of_friends[(joined_df_of_friends['expression_quartile'] == 3)]




# Perform KS test on two datasets (data1 and data2)
#Q1 v Q2
ks_statistic, p_value = ks_2samp(joined_df_of_friendsQ1['log2_epb_gene'], joined_df_of_friendsQ2['log2_epb_gene'])

print("KS statistic:", ks_statistic)
print("p-value:", p_value)

#Q2 v Q3
ks_statistic, p_value = ks_2samp(joined_df_of_friendsQ2['log2_epb_gene'], joined_df_of_friendsQ3['log2_epb_gene'])

print("KS statistic:", ks_statistic)
print("p-value:", p_value)


#Q3 v Q4
ks_statistic, p_value = ks_2samp(joined_df_of_friendsQ3['log2_epb_gene'], joined_df_of_friendsQ4['log2_epb_gene'])

print("KS statistic:", ks_statistic)
print("p-value:", p_value)




sns.lmplot(data=joined_df_of_friends, x="log2_rna_IGR", y="log2_epb_IGR")



# Sample data (replace this with your actual data)
data_q1 = joined_df_of_friendsQ1['log2_epb_gene']
data_q2 = joined_df_of_friendsQ2['log2_epb_gene']
data_q3 = joined_df_of_friendsQ3['log2_epb_gene']
data_q4 = joined_df_of_friendsQ4['log2_epb_gene']

# Compute KS statistic and p-value
ks_statistic, p_value = ks_2samp(data_q3, data_q4)

# Plot ECDFs
plt.figure(figsize=(10, 6))

sns.ecdfplot(data=data_q1, label='Q1', complementary=True)
sns.ecdfplot(data=data_q2, label='Q2', complementary=True)
sns.ecdfplot(data=data_q3, label='Q3', complementary=True)
sns.ecdfplot(data=data_q4, label='Q4', complementary=True)
plt.xlim(-7.5,-4.5)


# Highlight KS statistic with a vertical line
#plt.axvline(ks_statistic, color='red', linestyle='dashed', label=f'KS Statistic = {ks_statistic:.4f}')

plt.title('ECDF Plot with KS Statistic')
plt.xlabel('Data')
plt.ylabel('Cumulative Proportion')
plt.legend()

# Print KS statistic and p-value
plt.annotate(f'KS Statistic = {ks_statistic:.4f}\np-value = {p_value:.4f}',
             xy=(ks_statistic, 0.5), xytext=(ks_statistic + 0.5, 0.6),
             arrowprops=dict(facecolor='black', arrowstyle='->'))

plt.show()

quartiles = [joined_df_of_friendsQ1, joined_df_of_friendsQ2, joined_df_of_friendsQ3, joined_df_of_friendsQ4]

ks_statistics = []
p_values = []

for i in range(len(quartiles) - 1):
    ks_statistic, p_value = ks_2samp(quartiles[i]['log2_epb_gene'], quartiles[i + 1]['log2_epb_gene'])
    ks_statistics.append(ks_statistic)
    p_values.append(p_value)

# Create a figure and axis
custom_colors_bars = ['#87A2FB', '#9D84B7', '#FFAD60']
fig, ax = plt.subplots(1, 2, figsize=(10, 8))

# Plot KS statistics
ax[0].bar(range(1, len(quartiles)), ks_statistics, color=custom_colors_bars)
ax[0].set_xticks(range(1, len(quartiles)))
ax[0].set_xticklabels([f'Q{i + 1} vs Q{i + 2}' for i in range(len(quartiles) - 1)])
ax[0].set_xlabel('Quartile Comparison')
ax[0].set_ylabel('KS Statistic')
ax[0].set_title('KS Statistic for Each Quartile Comparison')
ax[0].spines['top'].set_visible(False)
ax[0].spines['right'].set_visible(False)


# Plot p-values
ax[1].bar(range(1, len(quartiles)), -np.log10(p_values), color=custom_colors_bars)
ax[1].set_xticks(range(1, len(quartiles)))
ax[1].set_xticklabels([f'Q{i + 1} vs Q{i + 2}' for i in range(len(quartiles) - 1)])
ax[1].set_xlabel('Quartile Comparison')
ax[1].set_ylabel('-log10(p-value)')
ax[1].set_title('-log10(p-value) for Quartile Comparisons')
ax[1].spines['top'].set_visible(False)
ax[1].spines['right'].set_visible(False)


#plt.tight_layout()
plt.show()


data = {
    'Quartile Comparison': [f'Q{i + 1} vs Q{i + 2}' for i in range(len(quartiles) - 1)],
    'KS Statistic': ks_statistics,
    '-log10(p-value)': p_values,
    'Category': ['KS'] * len(ks_statistics) + ['p val'] * len(p_values)
}

data = {
    'Quartile Comparison': [f'Q{i + 1} vs Q{i + 2}' for i in range(len(quartiles) - 1)],
    'KS Statistic': ks_statistics,
    '-log10(p-value)': p_values,
    'Category': ['KS'] * (len(quartiles) - 1) + ['p val'] * (len(quartiles) - 1)
}

3+3


df = pd.DataFrame(data)

# Create a barplot using Seaborn with hue
sns.set(style="whitegrid")
plt.figure(figsize=(12, 6))
sns.barplot(data=df, x='Quartile Comparison', y='KS Statistic', hue='Category', palette=['teal', 'orange'])
plt.xlabel('Quartile Comparison')
plt.ylabel('KS Statistic')
plt.title('KS Statistic and p-value for Each Quartile Comparison')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.show()





data_q1 = df_cleanedgy['log2_epb']
data_q2 = df_cleanedgy['log2_epb']
data_q3 = df_cleanedgy['log2_epb']
data_q4 = df_cleanedgy['log2_epb']

# Compute KS statistic and p-value
ks_statistic, p_value = ks_2samp(data_q3, data_q4)

# Plot ECDFs
plt.figure(figsize=(10, 6))

sns.ecdfplot(data=data_q1, label='Q1', complementary=True)
sns.ecdfplot(data=data_q2, label='Q2', complementary=True)
sns.ecdfplot(data=data_q3, label='Q3', complementary=True)
sns.ecdfplot(data=data_q4, label='Q4', complementary=True)



sns.histplot(
    data=joined_df_of_friends, y="log2_epb_gene", hue="expression_quartile",
    element="step",  common_norm = False)
plt.ylim(-7.5,-4.5)


plt.figure(figsize=(10, 4))
sns.kdeplot(
   data=df_cleanedgy, x="log2_epb", hue="expression_quartiley",
   fill=True, common_norm=False,
   alpha=.5, linewidth=1,
)
plt.xlim(-4.5,-7.5)



####
print('bitch')

sns.lmplot(
    data=joined_df_of_friends, x="log2_rna", y="log2_epb_IGR",
    col="expression_quartile", height=3,
    facet_kws=dict(sharex=False, sharey=False),
)

#%%%



df_cleaned1 = df_cleaned[(df_cleaned['expression_quartile'] == 'Q1')]
df_cleaned2 = df_cleaned[(df_cleaned['expression_quartile'] == 'Q2')]
df_cleaned3 = df_cleaned[(df_cleaned['expression_quartile'] == 'Q3')]
df_cleaned4 = df_cleaned[(df_cleaned['expression_quartile'] == 'Q4')]

data_q1 = df_cleaned1['log2_epb']
data_q2 = df_cleaned2['log2_epb']
data_q3 = df_cleaned3['log2_epb']
data_q4 = df_cleaned4['log2_epb']



df_cleanedgy1 = df_cleanedgy[(df_cleanedgy['expression_quartiley'] == 'Q1')]
df_cleanedgy2 = df_cleanedgy[(df_cleanedgy['expression_quartiley'] == 'Q2')]
df_cleanedgy3 = df_cleanedgy[(df_cleanedgy['expression_quartiley'] == 'Q3')]
df_cleanedgy4 = df_cleanedgy[(df_cleanedgy['expression_quartiley'] == 'Q4')]

data_q1 = df_cleanedgy1['log2_epb']
data_q2 = df_cleanedgy2['log2_epb']
data_q3 = df_cleanedgy3['log2_epb']
data_q4 = df_cleanedgy4['log2_epb']




# Generate sample data for two distributions
np.random.seed(123)
distribution1 = np.random.normal(0, 1, 1000)
distribution2 = np.random.normal(0.5, 1.2, 1000)

# Calculate ECDF for each distribution
def ecdf(data):
    x = np.sort(data)
    y = np.arange(1, len(data) + 1) / len(data)
    return x, y

x1, y1 = ecdf(data_q1)
x2, y2 = ecdf(data_q2)
x3, y3 = ecdf(data_q3)
x4, y4 = ecdf(data_q4)

# Perform KS-2samp test (IGR)
ks_statistic14, p_value14 = ks_2samp(data_q1, data_q4)
ks_statistic24, p_value24 = ks_2samp(data_q2, data_q4)
ks_statistic34, p_value34 = ks_2samp(data_q3, data_q4)
ks_statistic13, p_value13 = ks_2samp(data_q1, data_q3)

#for genes
ks_statistic14, p_value14 = ks_2samp(data_q1, data_q4)
ks_statistic24, p_value24 = ks_2samp(data_q2, data_q4)
ks_statistic34, p_value34 = ks_2samp(data_q3, data_q4)



# Create ECDF plot GENE
plt.figure(figsize=(10, 4))

plt.subplot(1, 3, 1)
plt.step(x1, y1, label='Q1', color='#EF6262')
plt.step(x3, y3, label='Q4', color='#F3AA60')
plt.xlabel('log2_epb_gene')
plt.ylabel('ECDF')
plt.title('Q1 vs. Q3')
plt.annotate("KS = {:.3f}".format(ks_statistic13), (-6,0.07))
#plt.annotate("p-value = {:.5f}".format(p_value34), (-6,0.03))
plt.legend()

plt.subplot(1, 3, 2)
plt.step(x2, y2, label='Q2', color='#468B97')
plt.step(x4, y4, label='Q4', color='#F3AA60')
plt.annotate("KS = {:.3f}".format(ks_statistic24), (-6,0.07))
#plt.annotate("p-value = {:.9f}".format(p_value24), (-7,0.03))
plt.xlabel('log2_epb_gene')
plt.ylabel('ECDF')
plt.title('Q2 vs. Q4')
plt.legend()

plt.subplot(1, 3, 3)
plt.step(x1, y1, label='Q1', color='#1D5B79')
plt.step(x4, y4, label='Q4', color='#F3AA60')
plt.annotate("KS = {:.3f}".format(ks_statistic14), (-6,0.07))
#plt.annotate("p-value = {:.14f}".format(p_value14), (-9,0.03))
plt.xlabel('log2_epb_gene')
plt.ylabel('ECDF')
plt.title('Q1 vs. Q4')
plt.legend()

##ECDF PLOT IGR

ks_statistic23, p_value23 = ks_2samp(data_q2, data_q3)
ks_statistic13, p_value13 = ks_2samp(data_q1, data_q3)




plt.figure(figsize=(10, 4))

plt.subplot(1, 3, 1)
plt.step(x1, y1, label='Q1', color='#1D5B79')

plt.step(x2, y2, label='Q2', color='#468B97')

plt.step(x4, y4, label='Q4', color='#F3AA60')
plt.xlabel('log2_epb_gene')
plt.ylabel('ECDF')
plt.title('Q1/Q2 vs. Q4')
plt.annotate("KS (Q1 vs. Q4) = {:.3f}".format(ks_statistic14), (-10,0.22))
#plt.annotate("p-value (Q1 vs. Q4) = {:.9f}".format(p_value14), (-9,0.18))
plt.annotate("KS (Q2 vs. Q4) = {:.3f}".format(ks_statistic24), (-10,0.14))
#plt.annotate("p-value (Q2 vs. Q4) = {:.14f}".format(p_value24), (-9,0.1))
#plt.xlim(-7.5,-4.5)
plt.legend()

plt.subplot(1, 3, 2)
plt.step(x1, y1, label='Q1', color='#1D5B79')
plt.step(x3, y3, label='Q3', color='#EF6262')
plt.annotate("KS = {:.3f}".format(ks_statistic13), (-10,0.22))
#plt.annotate("p-value = {:.5f}".format(p_value13), (-8.5,0.58))
plt.xlabel('log2_epb_gene')
plt.ylabel('ECDF')
plt.title('Q1 vs. Q3')
#plt.xlim(-7.5,-4.5)
plt.legend()

plt.subplot(1, 3, 3)
plt.step(x2, y2, label='Q2', color='#468B97')
plt.step(x3, y3, label='Q3', color='#EF6262')
plt.annotate("KS = {:.3f}".format(ks_statistic23), (-10,0.22))
#plt.annotate("p-value = {:.3f}".format(p_value23), (-6.5,0.58))
plt.xlabel('log2_epb_gene')
plt.ylabel('ECDF')
plt.title('Q2 vs. Q3')
#plt.xlim(-7.5,-4.5)
plt.legend()


custom_palette = ["#1D5B79", "#468B97", "#EF6262", "#F3AA60"]




# Create KS-2samp result plot
plt.subplot(1, 2, 2)
plt.bar(0, ks_statistic, color='blue', label='KS Statistic')
plt.bar(1, p_value, color='orange', label='p-value')
plt.xticks([0, 1], ['KS Statistic', 'p-value'])
plt.ylabel('Value')
plt.title('KS-2samp Test Result')
plt.legend()

plt.tight_layout()
plt.show()

print("KS Statistic:", ks_statistic)
print("p-value:", p_value)




#%%



#write a subset for let's say 1kb either side of a gene
def flank(file,frame1):
    data_output = pd.DataFrame()
    data_output_cnorm = pd.DataFrame()
    for i in range(len(file)):
        tempstart = file.iloc[i]["Start position"]
        tempchro = file.iloc[i]["Chr"]
        tempstrand = file.iloc[i]['Strand']
        
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

flankwtetfor, flankwtetforC = flank(ffeat, wtET)
flankt2efor, flankt2eforC = flank(ffeat, t2E)


#%%

merged_featQ1 = df_cleanedgy[(df_cleanedgy['expression_quartiley'] == 'Q1')]
merged_featQ2 = df_cleanedgy[(df_cleanedgy['expression_quartiley'] == 'Q2')]
merged_featQ3 = df_cleanedgy[(df_cleanedgy['expression_quartiley'] == 'Q3')]
merged_featQ4 = df_cleanedgy[(df_cleanedgy['expression_quartiley'] == 'Q4')]





#VERSION 2 INCLUDE REVERSE STRAND GENES
def flank(file,frame1):
    data_output = pd.DataFrame()
    data_output_cnorm = pd.DataFrame()
    
    for i in range(len(file)):
        tempstrand = file.iloc[i]['Strand']
        tempchro = file.iloc[i]["Chr"]
        if tempstrand == 'forward':
            tempstart = file.iloc[i]["Start position"]
            tempchro = file.iloc[i]["Chr"]

            startminus = (tempstart) - 1000
            tempstarty = (tempstart) + 500

        
            tempsubsetfor = frame1.loc[(frame1['Pos']>= startminus) & (frame1['Pos'] <= tempstarty) & (frame1['Chr'] == tempchro)]

            tempsubsetfor = tempsubsetfor.reset_index(drop = True)
            line = pd.DataFrame({"Pos": tempstarty}, index=[1500])
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
    
                "aPos" : np.arange(0.0, 1500, 1)}).merge(
        
                    tempsubsetfor, on = "aPos", how = "left"

                    )

            data_blank = data_blank[["aPos", "Wnorm", "Cnorm"]]
            
            data_output = pd.concat([data_output, data_blank['Wnorm']], axis =1, ignore_index=True)
            data_output_cnorm = pd.concat([data_output_cnorm, data_blank['Cnorm']], axis =1, ignore_index=True)
            
            data_output = data_output.fillna(0)

            data_output_cnorm = data_output_cnorm.fillna(0)
        
        
        if tempstrand == 'reverse':
            tempstartr = file.iloc[i]["End position"]            

            startminusr = (tempstartr) + 1000
            tempstartyr = (tempstartr) - 500

        
            tempsubsetrev = frame1.loc[(frame1['Pos']>= tempstartyr) & (frame1['Pos'] <= startminusr) & (frame1['Chr'] == tempchro)]

            tempsubsetrev = tempsubsetrev.reset_index(drop = True)
            line = pd.DataFrame({"Pos": tempstartyr}, index=[1500])
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
    
                "aPos" : np.arange(0.0, 1500, 1)}).merge(
        
                    tempsubsetrev, on = "aPos", how = "left"

                    )

            data_blankr = data_blankr[["aPos", "Wnorm", "Cnorm"]]
            
            data_output = pd.concat([data_output, data_blankr['Wnorm']], axis =1, ignore_index=True)
            data_output_cnorm = pd.concat([data_output_cnorm, data_blankr['Cnorm']], axis =1, ignore_index=True)
            
            data_output = data_output.fillna(0)

            data_output_cnorm = data_output_cnorm.fillna(0)

    return data_output,data_output_cnorm

flankwtetfor, flankwtetforC = flank(ffeat, wtET)
flankt2efor, flankt2eforC = flank(ffeat, t2E)


q1flank, q1flankC = flank(merged_featQ1, wtET)
q2flank, q2flankC = flank(merged_featQ2, wtET)
q3flank, q3flankC = flank(merged_featQ3, wtET)
q4flank, q4flankC = flank(merged_featQ4, wtET)





def Findfeat(file):
    genes = pd.read_csv(file, delimiter="\t")
    print(genes)


    return genes

ori = Findfeat('pombase_origins.tsv')

ori.rename(columns={'chr':'Chr', 'start':'Start position', 'end':'End position'}, inplace = True)


def flanko(file,frame1):
    data_output = pd.DataFrame()
   # data_output_cnorm = pd.DataFrame()
    
    for i in range(len(file)):

        tempstart = file.iloc[i]["pos"]
        tempchro = file.iloc[i]["Chr"]

        startminus = (tempstart) - 2000
        tempstarty = (tempstart) + 2000

        
        tempsubsetfor = frame1.loc[(frame1['Pos']>= startminus) & (frame1['Pos'] <= tempstarty) & (frame1['Chr'] == tempchro)]

        tempsubsetfor = tempsubsetfor.reset_index(drop = True)
        line = pd.DataFrame({"Pos": tempstarty}, index=[4000])
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
    
            "aPos" : np.arange(0.0, 4000, 1)}).merge(
        
                tempsubsetfor, on = "aPos", how = "left"

                )

        data_blank = data_blank[["aPos", "normscore"]]
            
        data_output = pd.concat([data_output, data_blank['normscore']], axis =1, ignore_index=True)
         #   data_output_cnorm = pd.concat([data_output_cnorm, data_blank['Cnorm']], axis =1, ignore_index=True)
            
        data_output = data_output.fillna(0)

         #   data_output_cnorm = data_output_cnorm.fillna(0)
        
     

    return data_output

oriwteflank = flanko(merged_ori, wtET)
orit2eflank = flanko(merged_ori, t2E)

pu_ori1 = owtchrI[owtchrI['origin'] > 0.2]
pu_ori2 = owtchrII[owtchrII['origin'] > 0.2]
pu_ori3 = owtchrIII[owtchrIII['origin'] > 0.2]




#VERSION 3 INCLUDE REVERSE STRAND GENES and ttake normscore
def flank(file,frame1):
    data_output = pd.DataFrame()
   # data_output_cnorm = pd.DataFrame()
    
    for i in range(len(file)):
        tempstrand = file.iloc[i]['Strand']
        tempchro = file.iloc[i]["Chr"]
        if tempstrand == 'forward':
            tempstart = file.iloc[i]["Start position"]
            tempchro = file.iloc[i]["Chr"]

            startminus = (tempstart) - 1000
            tempstarty = (tempstart) + 1000

        
            tempsubsetfor = frame1.loc[(frame1['Pos']>= startminus) & (frame1['Pos'] <= tempstarty) & (frame1['Chr'] == tempchro)]

            tempsubsetfor = tempsubsetfor.reset_index(drop = True)
            line = pd.DataFrame({"Pos": tempstarty}, index=[2000])
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
    
                "aPos" : np.arange(0.0, 2000, 1)}).merge(
        
                    tempsubsetfor, on = "aPos", how = "left"

                    )

            data_blank = data_blank[["aPos", "normscore"]]
            
            data_output = pd.concat([data_output, data_blank['normscore']], axis =1, ignore_index=True)
         #   data_output_cnorm = pd.concat([data_output_cnorm, data_blank['Cnorm']], axis =1, ignore_index=True)
            
            data_output = data_output.fillna(0)

         #   data_output_cnorm = data_output_cnorm.fillna(0)
        
        
        if tempstrand == 'reverse':
            tempstartr = file.iloc[i]["End position"]            

            startminusr = (tempstartr) + 1000
            tempstartyr = (tempstartr) - 1000

        
            tempsubsetrev = frame1.loc[(frame1['Pos']>= tempstartyr) & (frame1['Pos'] <= startminusr) & (frame1['Chr'] == tempchro)]

            tempsubsetrev = tempsubsetrev.reset_index(drop = True)
            line = pd.DataFrame({"Pos": tempstartyr}, index=[2000])
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
    
                "aPos" : np.arange(0.0, 2000, 1)}).merge(
        
                    tempsubsetrev, on = "aPos", how = "left"

                    )

            data_blankr = data_blankr[["aPos", "normscore"]]
            
            data_output = pd.concat([data_output, data_blankr['normscore']], axis =1, ignore_index=True)
           # data_output_cnorm = pd.concat([data_output_cnorm, data_blankr['Cnorm']], axis =1, ignore_index=True)
            
            data_output = data_output.fillna(0)

           # data_output_cnorm = data_output_cnorm.fillna(0)

    return data_output

#%%

#######rewrite without append

#VERSION 3 INCLUDE REVERSE STRAND GENES and ttake normscore
def flank(file,frame1):
    data_output = pd.DataFrame()
   # data_output_cnorm = pd.DataFrame()
    
    for i in range(len(file)):
        tempstrand = file.iloc[i]['Strand']
        tempchro = file.iloc[i]["Chr"]
        if tempstrand == 'forward':
            tempstart = file.iloc[i]["Start position"]
            tempchro = file.iloc[i]["Chr"]

            startminus = (tempstart) - 1000
            tempstarty = (tempstart) + 1000

        
            tempsubsetfor = frame1.loc[(frame1['Pos']>= startminus) & (frame1['Pos'] <= tempstarty) & (frame1['Chr'] == tempchro)]

            tempsubsetfor = tempsubsetfor.reset_index(drop = True)
            line = pd.DataFrame({"Pos": tempstarty}, index=[2000])
            line2 = pd.DataFrame({"Pos": startminus}, index=[0])
            tempsubsetfor = pd.concat([tempsubsetfor, line, line2], axis=0, ignore_index=False)
            
            tempsubsetfor = tempsubsetfor.sort_values('Pos', ascending=True)

            for x in range(len(tempsubsetfor)):
                pos = tempsubsetfor.iloc[x]['Pos']

                tpos = tempsubsetfor.iloc[0]['Pos']

                tempsubsetfor.loc[tempsubsetfor.index[x], 'aPos'] = (pos - tpos)

            tempsubsetfor = tempsubsetfor.set_index('aPos')
                        

            data_blank = pd.DataFrame({
    
                "aPos" : np.arange(0.0, 2000, 1)}).merge(
        
                    tempsubsetfor, on = "aPos", how = "left"

                    )

            data_blank = data_blank[["aPos", "normscore"]]
            
            data_output = pd.concat([data_output, data_blank['normscore']], axis =1, ignore_index=True)
         #   data_output_cnorm = pd.concat([data_output_cnorm, data_blank['Cnorm']], axis =1, ignore_index=True)
            
            data_output = data_output.fillna(0)

         #   data_output_cnorm = data_output_cnorm.fillna(0)
        
        
        if tempstrand == 'reverse':
            tempstartr = file.iloc[i]["End position"]            

            startminusr = (tempstartr) + 1000
            tempstartyr = (tempstartr) - 1000

        
            tempsubsetrev = frame1.loc[(frame1['Pos']>= tempstartyr) & (frame1['Pos'] <= startminusr) & (frame1['Chr'] == tempchro)]

            tempsubsetrev = tempsubsetrev.reset_index(drop = True)
            line = pd.DataFrame({"Pos": tempstartyr}, index=[2000])
            line2 = pd.DataFrame({"Pos": startminusr}, index=[0])
            tempsubsetrev = pd.concat([tempsubsetrev, line, line2], axis=0, ignore_index=False)
            
            tempsubsetrev = tempsubsetrev.sort_values('Pos', ascending=True)

            for x in range(len(tempsubsetrev)):
                pos = tempsubsetrev.iloc[0]['Pos']

                tpos = tempsubsetrev.iloc[x]['Pos']

                tempsubsetrev.loc[tempsubsetrev.index[x], 'aPos'] = (pos - tpos)

            tempsubsetrev = tempsubsetrev.set_index('aPos')
                        

            data_blankr = pd.DataFrame({
    
                "aPos" : np.arange(0.0, 2000, 1)}).merge(
        
                    tempsubsetrev, on = "aPos", how = "left"

                    )

            data_blankr = data_blankr[["aPos", "normscore"]]
            
            data_output = pd.concat([data_output, data_blankr['normscore']], axis =1, ignore_index=True)
           # data_output_cnorm = pd.concat([data_output_cnorm, data_blankr['Cnorm']], axis =1, ignore_index=True)
            
            data_output = data_output.fillna(0)

           # data_output_cnorm = data_output_cnorm.fillna(0)

    return data_output

flankwtetfor = flank(df_cleanedgy, wtET)
flankt2efor = flank(df_cleanedgy, t2E)
flanktrnawte = flank(tRNA, wtET)
flanktrnat2e = flank(tRNA, t2E)

flankpromoter5wte = flank(promoter, wtET)
flankpromoter5t2e = flank(promoter, t2E)

flankUTR5wte = flank(utr5, wtET)
flankUTR5t2e = flank(utr5, t2E)

ncRNAflankwte = flank(ncRNA, wtET)
ncRNAflankt2e = flank(ncRNA, t2E)

snRNAflankwte = flank(snRNA, wtET)
snRNAflankt2e = flank(snRNA, t2E)

snoRNAflankwte = flank(snoRNA, wtET)
snoRNAflankt2e = flank(snoRNA, t2E)

rRNAflankwte = flank(rRNA, wtET)
rRNAflankt2e = flank(rRNA, t2E)

promoterflankwte = flank(promoter, wtET)
promoterflankt2e = flank(promoter, t2E)



centroflanktrnawte = flank(df_within_centromeres, wtET)
centroflanktrnat2e = flank(df_within_centromeres, t2E)


q1flank = flank(merged_featQ1, wtET)
q2flank = flank(merged_featQ2, wtET)
q3flank = flank(merged_featQ3, wtET)
q4flank = flank(merged_featQ4, wtET)




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

wteW = conint(flankwtetfor)
wteC = conint(flankwtetforC)
t2eW = conint(flankt2efor)
t2eC = conint(flankt2eforC)
q1f = conint(q1flank)
q2f = conint(q2flank)
q3f = conint(q3flank)
q4f = conint(q4flank)

trnaWT = conint(flanktrnawte)
trnaT2E = conint(flanktrnat2e)





fixed, (ax1) =plt.subplots(1, sharey=True)
x = np.arange(0, 2001, 1)

y9ysc = np.array(q4f['upper']).flatten()
y10ysc= np.array(q4f['lower']).flatten()
ax1.fill_between(x, y9ysc, y10ysc, facecolor='#F3AA60', alpha =0.5)

y1sc = np.array(q3f['upper']).flatten()
y2sc= np.array(q3f['lower']).flatten()
ax1.fill_between(x, y1sc, y2sc, facecolor='#EF6262', alpha =0.5) 

y9sc = np.array(q2f['upper']).flatten()
y10sc= np.array(q2f['lower']).flatten()
ax1.fill_between(x, y9sc, y10sc, facecolor='#468B97', alpha =0.5)

y9s = np.array(q1f['upper']).flatten()
y10s= np.array(q1f['lower']).flatten()
ax1.fill_between(x, y9s, y10s, facecolor='#1D5B79', alpha =0.5)




fixed, (ax1) =plt.subplots(1, sharey=True)
x = np.arange(0, 2001, 1)
ax1.set_xticks([0,1000,2000])
ax1.set_xticklabels(['-1Kb','TSS','+1kb'])
    

y9ysc = np.array(wteW['upper']).flatten()
y10ysc= np.array(wteW['lower']).flatten()
ax1.fill_between(x, y9ysc, y10ysc, facecolor='#1D5B79', alpha = 0.5)
ax1.plot(flankwtetfor.index.values, flankwtetfor['smoo_mean'], color = 'midnightblue', alpha= 0.8)
y1sc = np.array(t2eW['upper']).flatten()
y2sc= np.array(t2eW['lower']).flatten()
ax1.fill_between(x, y1sc, y2sc, facecolor='#EF6262', alpha =0.5) 
ax1.plot(flankt2efor.index.values, flankt2efor['smoo_mean'], color = 'darkred', alpha= 0.8)






fixed, (ax1) =plt.subplots(1, sharey=True)
x = np.arange(0, 2001, 1)
ax1.set_xticks([0,1000,2000])
ax1.set_xticklabels(['-1Kb','TSS','+1kb'])
    

y9ysc = np.array(trnaWT['upper']).flatten()
y10ysc= np.array(trnaWT['lower']).flatten()
ax1.fill_between(x, y9ysc, y10ysc, facecolor='#1D5B79', alpha = 0.4)
ax1.plot(flanktrnawte.index.values, flanktrnawte['smoo_mean'], color = '#1D5B79', alpha= 0.95)
y1sc = np.array(trnaT2E['upper']).flatten()
y2sc= np.array(trnaT2E['lower']).flatten()
ax1.fill_between(x, y1sc, y2sc, facecolor='#EF6262', alpha =0.4) 
ax1.plot(flanktrnat2e.index.values, flanktrnat2e['smoo_mean'], color = 'darkred', alpha= 0.8)





def take_mean(all_data):
    all_data['mean'] = all_data.mean(axis=1)
    all_data['smoo_mean'] = all_data['mean'].rolling(window = 100, center=True).mean()
    
    return all_data

#wtforflanksenW = take_mean(wtforflanksenW, 100)
q1flank = take_mean(q1flank)
#q1flankC = take_mean(q1flankC)
q2flank = take_mean(q2flank)
#q2flankC = take_mean(q2flankC)
q3flank = take_mean(q3flank)
#q3flankC = take_mean(q3flankC)
q4flank = take_mean(q4flank)
#q4flankC = take_mean(q4flankC)
flankwtetfor = take_mean(flankwtetfor)
flankt2efor = take_mean(flankt2efor)

flanktrnawte = take_mean(flanktrnawte)
flanktrnat2e = take_mean(flanktrnat2e)

centroflanktrnawte = take_mean(centroflanktrnawte)
centroflanktrnat2e = take_mean(centroflanktrnat2e)
oriwteflank = take_mean(oriwteflank)
orit2eflank = take_mean(orit2eflank)

flankUTR5wte = take_mean(flankUTR5wte)
flankUTR5t2e = take_mean(flankUTR5t2e)


#flankpromoter5wte = take_mean(flankpromoter5wte)
#flankpromoter5t2e = take_mean(flankpromoter5t2e)


snRNAflankwte = take_mean(snRNAflankwte)
snRNAflankt2e = take_mean(snRNAflankt2e)

snoRNAflankwte = take_mean(snoRNAflankwte)
snoRNAflankt2e = take_mean(snoRNAflankt2e)

rRNAflankwte = take_mean(rRNAflankwte)
rRNAflankt2e = take_mean(rRNAflankt2e)



promoterflankwte = take_mean(promoterflankwte)
promoterflankt2e = take_mean(promoterflankt2e)





def plotall(ww, sw, title):
    fi,((ax1)) = plt.subplots(1,1, sharex=True, sharey=True) 

    ax1.plot(ww.index.values, ww['smoo_mean'], color = 'steelblue', alpha= 0.5)
  #  ax1.plot(wc.index.values, wc['smoo_mean'], color = 'steelblue', alpha= 0.5)
    ax1.set_title(title)
    ax1.set_ylabel('WT')


    ax1.plot(sw.index.values, sw['smoo_mean'], color = 'crimson', alpha= 0.5)
 #   ax1.plot(sc.index.values, sc['smoo_mean'], color = 'crimson', alpha= 0.5)
    ax1.set_ylabel('Top2-191')
    


    ax1.set_xticks([0,1000,2000])
    ax1.set_xticklabels(['-1Kb','TSS','+1kb'])
    

sen1forstalls = plotall(flankwtetfor, flankt2efor, 'Protein coding genes TSS flank')




def plotall(ww, wc, sw, sc, title):
    fi,((ax1)) = plt.subplots(1,1, sharex=True, sharey=True) 

    ax1.plot(ww.index.values, ww['smoo_mean'], color = 'steelblue', alpha= 0.5)
    ax1.plot(wc.index.values, wc['smoo_mean'], color = 'steelblue', alpha= 0.5)
    ax1.set_title(title)
    ax1.set_ylabel('WT')


    ax1.plot(sw.index.values, sw['smoo_mean'], color = 'crimson', alpha= 0.5)
    ax1.plot(sc.index.values, sc['smoo_mean'], color = 'crimson', alpha= 0.5)
    ax1.set_ylabel('Top2-191')
    


    ax1.set_xticks([0,500,1500])
    ax1.set_xticklabels(['+0.5Kb','TSS','-1kb'])
    

sen1forstalls = plotall(flankwtetfor, flankwtetforC, flankt2efor, flankt2eforC, 'Protein coding genes TSS flank')
dbl8forstalls = plotall(wtforflankdblW, wtforflankdblC, senforflankdblW, senforflankdblC, 'Protein coding genes TSS flank (Top2-191)')


# vww, vsw,
def plotall(ww, sw, title):
    fi,((ax1)) = plt.subplots(1,1, sharex=True, sharey=True) 

    ax1.plot(ww.index.values, ww['smoo_mean'], color = '#1D5B79')
   # ax1.plot(wc.index.values, wc['smoo_mean'], color = '#1D5B79', alpha= 0.5)
    ax1.set_title(title)
  #  ax1.set_ylabel('WT')


    ax1.plot(sw.index.values, sw['smoo_mean'], color = '#EF6262')
    #ax1.plot(sc.index.values, sc['smoo_mean'], color = '#468B97', alpha= 0.5)
   # ax1.set_ylabel('Top2-191')
   
  #  ax1.plot(vww.index.values, vww['smoo_mean'], color = '#EF6262')
    #ax1.plot(vwc.index.values, vwc['smoo_mean'], color = '#EF6262', alpha= 0.5)
    ax1.set_title(title)
  #  ax1.set_ylabel('WT')


  #  ax1.plot(vsw.index.values, vsw['smoo_mean'], color = '#F3AA60')
    #ax1.plot(vsc.index.values, vsc['smoo_mean'], color = '#F3AA60', alpha= 0.5)
   # ax1.set_ylabel('Top2-191')    

    ax1.set_ylabel('Aggregated CC-seq Signal (200bp smoothed)(HpM)')
    ax1.set_xlabel('relative TSS position (Kb)')
   # ax1.set_xticks([0,1000,2000, 3000, 4000])
   # ax1.set_xticklabels(['+2kb','+1Kb','origin','-1kb','-2kb'])
    ax1.set_xticklabels(['-1Kb','start','+1kb'])
    ax1.set_xticks([0,1000,2000])
  #  ax1.set_ylim(0.0003,0.00485)
   # 

sen1forstalls = plotall(q1flank, q1flankC, q2flank, q2flankC, q3flank, q3flankC, q4flank, q4flankC, 'Protein coding genes TSS flank')

sen1forstalls = plotall(q1flank, q2flank, q3flank, q4flank, 'Protein coding genes start site flank')

dbl8forstalls = plotall(wtforflankdblW, wtforflankdblC, senforflankdblW, senforflankdblC, 'Protein coding genes TSS flank (Top2-191)')

trnaflanks = plotall(flanktrnawte, flanktrnat2e, ' tRNA TSS flank')

promoterflanks = plotall(promoterflankwte, promoterflankt2e, 'promoter start flank')

originflanks = plotall(oriwteflank, oriwteflank,orit2eflank, orit2eflank, 'origin positions')

snRNAflanks = plotall(snRNAflankwte, snRNAflankt2e, 'snRNA start flank')

snoRNAflanks = plotall(snoRNAflankwte, snoRNAflankt2e, 'snoRNA start flank')

rRNAflanks = plotall(rRNAflankwte, rRNAflankt2e, 'rRNA start flank')


snoRNAflankwte = take_mean(snoRNAflankwte)
snoRNAflankt2e = take_mean(snoRNAflankt2e)

rRNAflankwte = take_mean(rRNAflankwte)
rRNAflankt2e = take_mean(rRNAflankt2e)






oriwteflank = take_mean(oriwteflank)
orit2eflank = take_mean(orit2eflank)

flankpromoter5wte = take_mean(flankpromoter5wte)
flankpromoter5t2e = take_mean(flankpromoter5t2e)


def plotall(ww, vww, title):
    fi,((ax1)) = plt.subplots(1,1, sharex=True, sharey=True) 

    ax1.plot(ww.index.values, ww['smoo_mean'], color = '#1D5B79')
   # ax1.plot(wc.index.values, wc['smoo_mean'], color = '#1D5B79', alpha= 0.5)
    ax1.set_title(title)
  #  ax1.set_ylabel('WT')
    ax1.plot(vww.index.values, vww['smoo_mean'], color = '#EF6262')
    #ax1.plot(vwc.index.values, vwc['smoo_mean'], color = '#EF6262', alpha= 0.5)
    ax1.set_title(title)
  #  ax1.set_ylabel('WT')


    ax1.set_ylabel('Aggregated CC-seq Signal (100bp smoothed)(HpM)')
    ax1.set_xlabel('position relative to 5 UTR (Kb)')
    ax1.set_xticks([0,1000,2000])
    ax1.set_xticklabels(['+1Kb','TSS','-1kb'])
  #  ax1.set_ylim(0.00022,0.00234)
    

sen1forstalls = plotall(flanktrnawte, flanktrnat2e, 'tRNA start position flank')
utrplot = plotall(flankUTR5wte, flankUTR5t2e, 'Protein coding genes 5 UTR')


["#1D5B79", "#468B97", "#EF6262", "#F3AA60"] 
trna1 = tRNA[(tRNA['Chroo'] == 'I')]
trna2 = tRNA[(tRNA['Chroo'] == 'II')]
trna3 = tRNA[(tRNA['Chroo'] == 'III')]




def check_plot (centro, genee, telo, ggr, frig):
    ff, (ax6) = plt.subplots(1,1, sharex=True)

    ax6.set_ylabel('Gene Annotations')
    ax6.set_xlabel('Chromosome position')

 
    for c in centro.itertuples(index=False, name=None):
            ax6.axvspan(c[0],c[1],0.1,0.15,color="black",alpha=0.8)
            
    for t in telo.itertuples(index=False, name=None):
            ax6.axvspan(t[0],t[1],0.1,0.15,color="grey",alpha=0.8)
            
    for ge in genee.itertuples(index=False, name=None):
        
        if ge[3] == 'protein coding':
            if ge[7] == 'reverse':
                ax6.axvspan(ge[4],ge[5],0.2,0.4,color="firebrick",alpha=0.3)
                #ax6.annotate(ge[0], xy = [ge[4],0.45]) 
            elif ge[7] == 'forward':
                ax6.axvspan(ge[4],ge[5],0.7,0.9,color="steelblue",alpha=0.3)
                #ax6.annotate(ge[0], xy = [ge[4],0.65]) 
                
    for gg in ggr.itertuples(index=False, name=None):

        ax6.axvspan(gg[1],gg[2],0.4,0.6,color="green",alpha=0.3)
    
    for fr in frig.itertuples(index=False, name=None):
        if fr[0] == 'I':

            ax6.axvspan(fr[1],fr[2],0.4,0.6,color="black",alpha=0.3)

    return ff

chromosome1 = check_plot(centro1, feat1, telo1, trna1, promoter)


common_ids = pd.merge(tRNA, df_within_centromeres, on='ID')['ID']

# Step 2: Drop rows from tRNA where 'ID' occurs in both tRNA and df_within_centromeres
tRNA = tRNA[~tRNA['ID'].isin(common_ids)]

# Print the updated tRNA DataFrame
print(tRNA)


#%%

def check_plot (centro, genee, telo, ggr):
    ff, (ax6) = plt.subplots(1,1, sharex=True)
   # ax6.set_ylim(0.2,0.8)

    ax6.set_ylabel('Gene Annotations')
    ax6.set_xlabel('Chromosome position')

 
    for c in centro.itertuples(index=False, name=None):
            ax6.axvspan(c[0],c[1],0.2,0.4,color="black",alpha=0.8)
            
    for t in telo.itertuples(index=False, name=None):
            ax6.axvspan(t[0],t[1],0.9,0.8,color="grey",alpha=0.8)
            
 #   for ge in genee.itertuples(index=False, name=None):
        
  #      if ge[3] == 'protein coding':
   #         if ge[7] == 'reverse':
    #            ax6.axvspan(ge[4],ge[5],0.2,0.4,color="firebrick",alpha=0.3)
                #ax6.annotate(ge[0], xy = [ge[4],0.45]) 
     #       elif ge[7] == 'forward':
      #          ax6.axvspan(ge[4],ge[5],0.7,0.9,color="steelblue",alpha=0.3)
                #ax6.annotate(ge[0], xy = [ge[4],0.65]) 
                
    for gg in ggr.itertuples(index=False, name=None):

        ax6.axvspan(gg[1],gg[2],0.3,0.8,color="green",alpha=0.3)

    return ff

chromosome3 = check_plot(centro, feat3,telo3, tRNA3)


tRNA1 = tRNA[(tRNA['Chroo'] == 'I')]
tRNA2 = tRNA[(tRNA['Chroo'] == 'II')]
tRNA3 = tRNA[(tRNA['Chroo'] == 'III')]
