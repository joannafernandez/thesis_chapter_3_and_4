#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 11:52:11 2023

@author: patricfernandez
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


import os


#%%
#cc loci extracted from FullMaps
#csvs renamed with Fullmap. prefix
#CC pipeline pombe seq bias 

#BASH for extraction 

root_directory="/Users/patricfernandez/Documents/911filter/sen1stallseqBias"

common_directory="/Users/patricfernandez/Documents/911filter/sen1stallseqBias"


# Iterate through the gene folders
for gene_folder in "$root_directory"/*; do
    if [ -d "$gene_folder" ]; then
        # Move all the text files from the gene folder to the common directory
        mv "$gene_folder"/*.txt "$common_directory/"
    fi
done

#move txt files to new folder, 
#replace '.' in txt file name with 'a' to get full spac ID 
#this new folder is the directory input to top3finder 

#%%

#try to make this into a function please joanna
def top3finder(directory):

    gene_status_df = pd.DataFrame(columns=['ID', 'Status'])
    
    file_dict = {}
    for filename in os.listdir(directory):
        if filename.endswith('.txt'):
            # Extract the gene name from the file name
            parts = filename.split('_')
            gene_name = parts[2]  # Assuming "aat1" is the third part
            file_path = os.path.join(directory, filename)
            
            # Add the gene name and file path to the dictionary
            with open(file_path, 'r') as file:
                #file_contents = file.read()
                file_contents = pd.read_csv(file, delimiter="\t")
                #C = file_contents.split("\t")
              # print(f"File Contents:\n{file_contents}\n")
    
                file_dict[gene_name] = file_contents
     
            
            for key, df in file_dict.items():
              #  print(key)
                interest = df.loc[15]
                interest = interest.drop(['Position'])
              #  print(interest.max())
                if interest.loc['G'] == interest.max():
                    is_largest = True
                elif  interest.loc['C'] == interest.max():
                    is_largest = True
                        
                else:
                #if interest.loc['G'] != interest.max() or interest.loc['C'] != interest.max(): 
                    is_largest = False
                
                        
            
                status = 'YES' if is_largest else 'NO'
                
                temp = pd.DataFrame({'ID': [key], 'Status': [status] })
    
               # print(temp)
            gene_status_df = pd.concat([gene_status_df, temp])
            gene_status_df['ID'] = gene_status_df['ID'].str.replace('a','.')
            
    return file_dict, gene_status_df
                
DScontrol, DScontrolTop3 = top3finder ("/Users/patricfernandez/Documents/seqbiasprocessed/DScontrolCORRECT")

WTcontrol, WTcontrolTop3 = top3finder("/Users/patricfernandez/Documents/seqbiasprocessed/WTcontrolCORRECT")

DSstall, DSstallTop3 = top3finder("/Users/patricfernandez/Documents/seqbiasprocessed/DSstallsCORRECT")

WTstall, WTstallTop3 = top3finder("/Users/patricfernandez/Documents/seqbiasprocessed/WTstallsCORRECT")



SEN1stall, SEN1stallTop3 = top3finder("/Users/patricfernandez/Documents/911filter/sen1seqbiasstalltxt")

DBL8stall, DBL8stallTop3 = top3finder("/Users/patricfernandez/Documents/911filter/dbl8seqbiasstalltxt")


SEN1control, SEN1controlTop3 = top3finder("/Users/patricfernandez/Documents/911filter/sen1seqbiascontroltxt")


DBL8control, DBL8controlTop3 = top3finder("/Users/patricfernandez/Documents/911filter/dbl8seqbiascontroltxt")


DSALLgenes, DSALLgenesTop3 = top3finder("/Users/patricfernandez/Documents/allgenestxt")

#%%

#attempt version with rank
def top3finder(directory):
    gene_status_df = pd.DataFrame(columns=['ID', 'Status', 'Rank'])
    file_dict = {}

    for filename in os.listdir(directory):
        if filename.endswith('.txt'):
            # Extract the gene name from the file name
            parts = filename.split('_')
            gene_name = parts[2]  # Assuming "aat1" is the third part
            file_path = os.path.join(directory, filename)

            # Add the gene name and file path to the dictionary
            with open(file_path, 'r') as file:
                file_contents = pd.read_csv(file, delimiter="\t")

            file_dict[gene_name] = file_contents

    for key, df in file_dict.items():
        interest = df.loc[15]
        interest = interest.drop(['Position'])

        # Create a DataFrame with the ID, Status, and Rank
        temp = pd.DataFrame({'ID': [key],
                             'Status': ['YES' if interest.loc['G'] == interest.max() or interest.loc['C'] == interest.max() else 'NO'],
                             'Rank': [interest.max()]})

        gene_status_df = pd.concat([gene_status_df, temp])

    # Rank the DataFrame based on the 'Rank' column in descending order
    gene_status_df = gene_status_df.sort_values(by='Rank', ascending=False).reset_index(drop=True)
    gene_status_df['ID'] = gene_status_df['ID'].str.replace('a', '.')

    return file_dict, gene_status_df

DSstall2, DSstallTop32 = top3finder("/Users/patricfernandez/Documents/seqbiasprocessed/DSstallsCORRECT")
DSALLgenes, DSALLgenesTop3 = top3finder("/Users/patricfernandez/Documents/allgenestxt")
WTALLgenes, WTALLgenesTop3 = top3finder('/Users/patricfernandez/Documents/wtalltxt')
#%%


def Find(file, delim):
    genes = pd.read_csv(file, delimiter=delim)
 #   print(genes)
   # genes['length'] = genes['end'] - genes['start']
  #  print(genes['length'])
    genesfor = genes[(genes['coding_strand'] == 'forward')]
    genesrev = genes[(genes['coding_strand'] == 'reverse')]
    
    genes.loc[genes['chro'] == "chr1", 'chro'] = 1
    genes.loc[genes['chro'] == "chr2", 'chro'] = 2
    genes.loc[genes['chro'] == "chr3", 'chro'] = 3


    return genesfor, genesrev, genes

ggenesfor, ggenesrev, gggenes = Find("dbl8_stall_sites_direction.txt", "\t")
controlfor, controlrev, ccontrol = Find('new_control.csv', ",")


#sen1stall = gggenes[(gggenes['genotype'] == 'sen1D')]
#dbl8stall = gggenes[(gggenes['genotype'] == 'dbl8D')]
#doublestall = gggenes[(gggenes['genotype'] == 'sen1dbl8DD_unique')]



#compiling df
gggenes = gggenes.merge(DSstallTop3, how='left', on='ID')

gggenes = gggenes.merge(WTstallTop3, how='left', on='ID')

gggenes = gggenes.rename(columns={'Status_x': 'DStop3', 'Status_y': 'WTtop3'})

gggenes = gggenes.merge(SEN1stallTop3, how='left', on='ID')

gggenes = gggenes.merge(DBL8stallTop3, how='left', on='ID')

gggenes = gggenes.rename(columns={'Status_x': 'SEN1top3', 'Status_y': 'DBL8top3'})


ccontrol = ccontrol.merge(DScontrolTop3, how='left', on='ID')

ccontrol = ccontrol.merge(WTcontrolTop3, how='left', on='ID')

ccontrol = ccontrol.rename(columns={'Status_x': 'DStop3', 'Status_y': 'WTtop3'})

ccontrol = ccontrol.merge(SEN1controlTop3, how='left', on='ID')

ccontrol = ccontrol.merge(DBL8controlTop3, how='left', on='ID')

ccontrol = ccontrol.rename(columns={'Status_x': 'SEN1top3', 'Status_y': 'DBL8top3'})

gggenes['length'] = (gggenes['end'] - gggenes['start'])


sen1stall = gggenes[(gggenes['genotype'] == 'sen1D')]
dbl8stall = gggenes[(gggenes['genotype'] == 'dbl8D')]
doublestall = gggenes[(gggenes['genotype'] == 'sen1dbl8DD_unique')]

condition = (
    (doublestall['genotype'] == 'sen1dbl8DD_unique') & 
    ((doublestall['coding_strand'] == 'reverse') & (doublestall['stalled_fork'] == 'leftward') | 
     (doublestall['coding_strand'] == 'forward') & (doublestall['stalled_fork'] == 'rightward'))
)

# Apply the condition to filter the DataFrame
HTstall = doublestall.loc[condition]


condition2 = (
    (doublestall['genotype'] == 'sen1dbl8DD_unique') & 
    ((doublestall['coding_strand'] == 'reverse') & (doublestall['stalled_fork'] == 'rightward') | 
     (doublestall['coding_strand'] == 'forward') & (doublestall['stalled_fork'] == 'leftward'))
)

# Apply the condition to filter the DataFrame
HHstall = doublestall.loc[condition2]



#proportion sen
sen1stallmelted_df = pd.melt(sen1stall, id_vars=['ID', 'gene', 'start', 'end'], value_vars=['DStop3', 'WTtop3', 'SEN1top3', 'DBL8top3'], var_name='Condition', value_name='Value')

sen1stallmelted_dfy = sen1stallmelted_df.copy()
sen1stallmelted_dfy["Condition"] = pd.Categorical(sen1stallmelted_dfy["Condition"],
                                                 ordered = True,
                                                 categories = ["WTtop3", "SEN1top3", "DBL8top3", "DStop3"])

custom_palette = ["#1D5B79", "#EF6262"]
# Plot the countplot
sns.set(style='whitegrid')
plt.figure(figsize=(10, 5))
sns.countplot(x='Condition', hue='Value', data=sen1stallmelted_dfy, palette = custom_palette, stat='count', dodge=True)
plt.title('Countplot of Conditions')
plt.show()



#proportion dbl8
sen1stallmelted_df = pd.melt(doublestall, id_vars=['ID', 'gene', 'start', 'end'], value_vars=['DStop3', 'WTtop3', 'SEN1top3', 'DBL8top3'], var_name='Condition', value_name='Value')

sen1stallmelted_dfy = sen1stallmelted_df.copy()
sen1stallmelted_dfy["Condition"] = pd.Categorical(sen1stallmelted_dfy["Condition"],
                                                 ordered = True,
                                                 categories = ["WTtop3", "SEN1top3", "DBL8top3", "DStop3"])

custom_palette = ["#1D5B79", "#EF6262"]
# Plot the countplot
sns.set(style='whitegrid')
plt.figure(figsize=(10, 5))
sns.countplot(x='Condition', hue='Value', data=sen1stallmelted_dfy, palette = custom_palette, stat='count', dodge=True)
plt.title('Countplot of Conditions')
plt.show()


#HTstall
sen1stallmelted_df = pd.melt(HTstall, id_vars=['ID', 'gene', 'start', 'end'], value_vars=['DStop3', 'WTtop3', 'SEN1top3', 'DBL8top3'], var_name='Condition', value_name='Value')

sen1stallmelted_dfy = sen1stallmelted_df.copy()
sen1stallmelted_dfy["Condition"] = pd.Categorical(sen1stallmelted_dfy["Condition"],
                                                 ordered = True,
                                                 categories = ["WTtop3", "SEN1top3", "DBL8top3", "DStop3"])

custom_palette = ["#1D5B79", "#EF6262"]
# Plot the countplot
sns.set(style='whitegrid')
plt.figure(figsize=(10, 5))
sns.countplot(x='Condition', hue='Value', data=sen1stallmelted_dfy, palette = custom_palette, stat='count', dodge=True)
plt.title('Countplot of Conditions')
plt.show()

#HHstall
sen1stallmelted_df = pd.melt(HHstall, id_vars=['ID', 'gene', 'start', 'end'], value_vars=['DStop3', 'WTtop3', 'SEN1top3', 'DBL8top3'], var_name='Condition', value_name='Value')

sen1stallmelted_dfy = sen1stallmelted_df.copy()
sen1stallmelted_dfy["Condition"] = pd.Categorical(sen1stallmelted_dfy["Condition"],
                                                 ordered = True,
                                                 categories = ["WTtop3", "SEN1top3", "DBL8top3", "DStop3"])

custom_palette = ["#1D5B79", "#EF6262"]
# Plot the countplot
sns.set(style='whitegrid')
plt.figure(figsize=(10, 5))
sns.countplot(x='Condition', hue='Value', data=sen1stallmelted_dfy, palette = custom_palette, stat='count', dodge=True)
plt.title('Countplot of Conditions')
plt.show()

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
ffeat = ffeat.merge(DSALLgenesTop3, how='left', on='ID')
ffeat = ffeat.merge(WTALLgenesTop3, how ='left', on='ID')
ffeat = ffeat.rename(columns={'Status_x': 'DStop3', 'Status_y': 'WTtop3'})
ffeat = ffeat.rename(columns={'Rank_x': 'Rank_DS', 'Rank_y': 'Rank_WT'})


def determine_category(row):
    if row['DStop3'] == 'YES' and row['WTtop3'] == 'YES':
        return 'bothYES'
    elif row['DStop3'] == 'YES' and row['WTtop3'] == 'NO':
        return 'DSonly'
    elif row['DStop3'] == 'NO' and row['WTtop3'] == 'YES':
        return 'WTonly'
    elif row['DStop3'] == 'NO' and row['WTtop3'] == 'NO':
        return 'bothNO'
    else:
        return 'Unknown'

ffeat['DET'] = ffeat.apply(determine_category, axis=1)
featwithrna['DET'] = featwithrna.apply(determine_category, axis=1)


if 'DSonly' in ffeat['DET'].values:
    print("DSonly is present in the 'DET' column.")
else:
    print("DSonly is not present in the 'DET' column.")

sns.set(style='whitegrid')
plt.figure(figsize=(10, 5))
sns.countplot(x='DET', data=ffeat, palette = custom_palette, stat='percent')
plt.title('Countplot of Conditions')
plt.show()


sns.set_style("ticks")



ffeat['Log2_length'] =np.log2(ffeat['end'] - ffeat['start'])
ffeat['Log2_Rank'] = np.log2(ffeat['Rank_DS'])

sns.lmplot(data=ffeat, x="Log2_length", y="Log2_Rank", hue="DStop3", scatter_kws={'s': 5})






xls = pd.ExcelFile('S_pombe_expression_data_PMID_23101633.xlsx')
df1 = pd.read_excel(xls, 'Table_S4')
df1.columns=df1.iloc[5]
df1 = df1.iloc[6: , :]
rob = df1[["Systematic.name", "MM.mRNA.cpc"]]
rob.rename({'Systematic.name': 'ID'}, axis=1, inplace=True)
rob.rename({"MM.mRNA.cpc": 'rna'}, axis=1, inplace=True)



featwithrna = rob.merge(ffeat, how='left', on='ID')
featwithrna = featwithrna.dropna(how='any', axis = 0)
featwithrna['rna'] = featwithrna['rna'].astype(int)


featwithrna['Log2_rna'] = np.log2(featwithrna['rna'])
featwithrna['Log2_Rank'] = np.log2(featwithrna['Rank'])
featwithrnay = featwithrna[featwithrna['rna'] != 0]
joinedd1 = featwithrnay[(featwithrnay['chro'] == 1)]
joinedd2 = featwithrnay[(featwithrnay['chro'] == 2)]
joinedd3 = featwithrnay[(featwithrnay['chro'] == 3)]

featwithrnay['Log2_length'] =np.log2(featwithrnay['end'] - featwithrnay['start'])
featwithrnay['length'] = featwithrnay['end'] - featwithrnay['start']

g = sns.JointGrid()
x, y = featwithrnay["Log2_rna"], featwithrnay["Log2_Rank"]
sns.regplot(x=x, y=y, scatter_kws={'s': 1}, line_kws={'color': 'black'}, ax=g.ax_joint)
sns.histplot(data=featwithrnay, x="Log2_rna", hue="Status",
    element="step",  common_norm = False, ax=g.ax_marg_x, legend = False)
sns.histplot(data=featwithrnay, y="Log2_Rank", hue="Status",
    element="step",  common_norm = False,  ax=g.ax_marg_y, legend = False)




custom_palette = ["#1D5B79", "#EF6262"]
sns.set_palette(custom_palette)

sns.lmplot(data=featwithrnay, x="length", y="Log2_Rank", hue="Status", scatter_kws={'s': 5})



#%%
#now let's try to do this but with my RNA-seq data.
#first write a function to read in significant DE 


def Collect_data(file1):
    full = pd.read_csv(file1)
    full.rename(columns={'Unnamed: 0': 'ID'}, inplace=True)
    full['ID'] = full['ID'].str.replace('a', '.')
      
    return full

wtvsds = Collect_data('results_table.csv')

#next we write a function to assign "Sig DE" to my ffeat df
#we want a column for pvalue and log2fc 

featwithrna = wtvsds.merge(ffeat, how='left', on='ID')
featwithrna = featwithrna.dropna(subset=['Rank_DS'])
featwithrna['Log2_Rank'] = np.log2(featwithrna['Rank_DS'])
sns.lmplot(data=featwithrna, x="Log2_rna", y="Log2_Rank", hue="DStop3", scatter_kws={'s': 5})


####HERE IS DIFFERENTIAL TOP3 AND RNA_SEQ
featwithrna = wtvsds.merge(ffeatdsonly, how='left', on='ID')
featwithrna = featwithrna.dropna(subset=['Rank_DS'])
featwithrna['Log2_Rank'] = np.log2(featwithrna['Rank_DS'])


#featwithrna = featwithrna.dropna(how='any', axis = 0)
#featwithrna['rna'] = featwithrna['rna'].astype(int)


featwithrna['Log2_length'] =np.log2(featwithrna['end'] - featwithrna['start'])
featwithrna['Log2_rna'] =np.log2(featwithrna['baseMean'])
featwithrna['length'] = featwithrna['end'] - featwithrna['start']

g = sns.JointGrid()
x, y = featwithrna["log2FoldChange"], featwithrna["Log2_Rank"]
sns.regplot(x=x, y=y, scatter_kws={'s': 1}, line_kws={'color': 'black'}, ax=g.ax_joint)
sns.histplot(data=featwithrna, x="log2FoldChange", hue="DET",
    element="step",  common_norm = False, ax=g.ax_marg_x, legend = False)
sns.histplot(data=featwithrna, y="Log2_Rank", hue="DET",
    element="step",  common_norm = False,  ax=g.ax_marg_y, legend = False)




custom_palette = ["#1D5B79", "#EF6262"]
sns.set_palette(custom_palette)

sns.lmplot(data=featwithrna, x="Log2_length", y="Log2_Rank", hue="Status", scatter_kws={'s': 5})


featwithrnasig = featwithrna[featwithrna['padj'] < 0.05]
sns.lmplot(data=featwithrnasig, x="Log2_rna", y="Log2_Rank", hue="Status", scatter_kws={'s': 5})

featwithrnasig['-logpadj'] = -np.log10(featwithrnasig['padj'])
featwithrna['-logpadj'] = -np.log10(featwithrna['padj'])


sns.scatterplot(data=featwithrnasig, x="log2FoldChange", y="-logpadj", hue = 'Status', size='Log2_Rank', sizes=(20, 200))
plt.show()

sns.scatterplot(data=featwithrna, x="log2FoldChange", y="-logpadj", hue = 'DET', size='Log2_Rank')
plt.show()


custom_palette = ['#1D5B79','#EF6262','#557C55', '#FF5B22']
custom_palette = ['#219C90','#E9B824','#EE9322', '#D83F31']

sns.set_palette(custom_palette)
sns.set_style("ticks")

# Your scatterplot code
sns.scatterplot(data=featwithrna, x="log2FoldChange", y="-logpadj", hue='DET', size='Log2_Rank')

# Add a horizontal dotted line at -log10(0.05)
plt.axhline(y=-np.log10(0.05), color='black', linestyle='--')

# Show the plot
plt.show()

custom_palette = ['#E9B824','#219C90','#EE9322', '#D83F31']


featwithrna['-logpadj'] = -np.log10(featwithrna['padj'])
sns.scatterplot(data=featwithrna, x="log2FoldChange", y="-logpadj", hue = 'Status', size='Log2_Rank', sizes=(20, 200))
plt.show()
sns.lmplot(data=featwithrna, x="Log2_rna", y="Log2_Rank",col ='DET', hue="DET", scatter_kws={'s': 5})


ffeatdsonly = featwithrna[(featwithrna['DET'] == 'DSonly')]
ffeatwtonly = featwithrna[(featwithrna['DET'] == 'WTonly')]
ffeatbothyes = featwithrna[(featwithrna['DET'] == 'bothYES')]
ffeatbothno = featwithrna[(featwithrna['DET'] == 'bothNO')]


import statsmodels.api as sm

# Assuming you have defined 'xx' and 'yy' as in your previous code
xx = np.array(ffeatbothno['Log2_rna']).reshape((-1, 1))
yy = np.array(ffeatbothno['Log2_Rank'])

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

featwithrnayYES = featwithrnay[(featwithrnay['Status'] == 'YES')]
featwithrnayNO = featwithrnay[(featwithrnay['Status'] == 'NO')]



plt.figure(figsize=(10, 4))
sns.kdeplot(
   data=featwithrnay, x="Log2_Rank",
   fill=True, common_norm=False,
   alpha=.5, linewidth=1,
)
#%%

featwithrnay = featwithrna[featwithrna['rna'] != 0]
joinedd1 = featwithrna[(featwithrna['chro'] == 1)]
joinedd2 = featwithrna[(featwithrna['chro'] == 2)]
joinedd3 = featwithrna[(featwithrna['chro'] == 3)]


joinedd1s = featwithrnasig[(featwithrnasig['chro'] == 1)]
joinedd2s = featwithrnasig[(featwithrnasig['chro'] == 2)]
joinedd3s = featwithrnasig[(featwithrnasig['chro'] == 3)]



#consider YES and NO on seperate graphs, and look at IGR type for hue
def last_resort (file):
    #data_output = pd.DataFrame()
    igr_start = []
    igr_end = []
    rnalist = []
    fclist = []
    statlist = []
    DETlist =[]
    ranklist = []

    igr_types = [] 
    name1 = []
    name2 = []
    name3 = []
    
    for i in range(len(file)-2):
        tempname = file.iloc[i]['ID']
        tempstart = file.iloc[i]["start"]
        tempend = file.iloc[i]["end"]
        tempstrand = file.iloc[i]['coding_strand']
        
        nextname = file.iloc[i+1]['ID']
        nextstart = file.iloc[i+1]["start"]
        nextend = file.iloc[i+1]["end"]
        nextstrand = file.iloc[i+1]['coding_strand']
        rnale = file.iloc[i+1]['Log2_rna']
        statale = file.iloc[i+1]['DStop3']
        defale = file.iloc[i+1]['DET']
        fcale = file.iloc[i+1]['log2FoldChange']
        rankale = file.iloc[i+1]['Log2_Rank']


        lastname = file.iloc[i+2]['ID']
        laststart = file.iloc[i+2]["start"]
        lastend = file.iloc[i+2]["end"]
        laststrand = file.iloc[i+2]['coding_strand']
        
        if nextstart - tempend < 0:
            #this is an overlap

            name2.append(nextname)
            rnalist.append(rnale)
            statlist.append(statale)
            DETlist.append(defale)
            fclist.append(fcale)
            ranklist.append(rankale)
            
            igr_types.append('overlap')
            igr_start.append(tempend)
            igr_end.append(nextstart)
            
                    
        elif nextstart > tempend and tempstrand == nextstrand:
            #this is a tandem igr 
            tandem_len = nextstart - tempend
            rnalist.append(rnale)
            statlist.append(statale)
            DETlist.append(defale)
            fclist.append(fcale)
            ranklist.append(rankale)
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
            statlist.append(statale)
            DETlist.append(defale)
            fclist.append(fcale)
            ranklist.append(rankale)
            divergentlen = nextstart - tempend
            igr_types.append('divergent')
            igr_start.append(tempend)
            igr_end.append(nextstart)

            
            
            
               # print(divergentlen)
                
        elif nextstart > tempend and tempstrand == 'forward' and nextstrand == 'reverse':
            #this is a convergent igr
            name2.append(nextname)
            rnalist.append(rnale)
            fclist.append(fcale)
            statlist.append(statale)
            DETlist.append(defale)
            ranklist.append(rankale)
            convergentlen = nextstart - tempend
            igr_types.append('convergent')
            igr_start.append(tempend)
            igr_end.append(nextstart)

            
            
            
    print(len(igr_start))
    print(len(igr_end))
    print(len(igr_types))
    print(len(nextname))
    print(len(rnalist))
    print(len(statlist))
    print(len(ranklist))
    print(len(fclist))
    igr_df = pd.DataFrame({
        'Start position': igr_start,
        'End position': igr_end,
        'IGR Type': igr_types,
        'IGR name': name2,
        'log2_rna': rnalist,
        'log2_Rank': ranklist,
        'Status': statlist,
        'log2FoldChange': fclist,
        'DET': DETlist
        
        
        # Add the IGR Type column
    })
        
    return igr_df

dye = last_resort(joinedd1)
dye2 = last_resort(joinedd2)     
dye3 = last_resort(joinedd3)     

dyes = last_resort(joinedd1s)
dye2s = last_resort(joinedd2s)     
dye3s = last_resort(joinedd3s)  
        
dye['Chr'] = 1
dye2['Chr'] = 2
dye3['Chr'] = 3

megeed_dye = pd.concat([dye, dye2, dye3], ignore_index=True)
megeed_dye = pd.concat([dyes, dye2s, dye3s], ignore_index=True)

megeed_dyey = megeed_dye.copy()
megeed_dyey["DET"] = pd.Categorical(sen1stallmelted_dfy["Condition"],
                                                 ordered = True,
                                                 categories = ["DSonly", "bothNO", "bothYES", "WTonly"])




custom_palette = ['#1D5B79','#EF6262','#557C55', '#FF5B22']
custom_palette = ['#219C90','#E9B824','#EE9322', '#D83F31']

sns.set_palette(custom_palette)
sns.set_style("ticks")
sns.lmplot(data=megeed_dyey, x="log2FoldChange", y="log2_Rank", col="IGR Type", hue = 'DET', scatter_kws={'s': 5}, facet_kws={'legend_out': True})



#####omg i think this is the better one
sns.lmplot(data=megeed_dyey, x="log2FoldChange", y="log2_Rank", col="DET", hue = 'IGR Type', scatter_kws={'s': 5}, facet_kws={'legend_out': True})
#try this
sns.catplot(data=megeed_dyey, x="DET", y="log2FoldChange", hue = 'IGR Type', kind="violin", inner=None)
sns.swarmplot(data=megeed_dyey, x="DET", y="log2FoldChange", hue = 'DET', size=1 )
sns.violinplot(data=megeed_dyey, x="DET", y="log2FoldChange", density_norm="count")
plt.show()
######I think what i need to do is only consider significant 

sns.lmplot(data=megeed_dye, x="log2_rna", y="log2_Rank", col="IGR Type", hue='DET', scatter_kws={'s': 5})


custom_palette = ["#EF6262", "#1D5B79"]
sns.set_palette(custom_palette)

sns.lmplot(data=megeed_dye, x="log2_rna", y="log2_Rank", col="IGR Type", hue='Status', scatter_kws={'s': 5})

# Show the plot
plt.show()
g = sns.JointGrid()
x, y = megeed_dye["log2_rna"], megeed_dye["log2_Rank"]
sns.regplot(x=x, y=y, scatter_kws={'s': 1}, line_kws={'color': 'black'}, ax=g.ax_joint)
sns.histplot(data=megeed_dye, x="log2_rna", hue="IGR Type",
    element="step",  common_norm = False, ax=g.ax_marg_x, legend = False)
sns.histplot(data=megeed_dye, y="log2_Rank", hue="IGR Type",
    element="step",  common_norm = False,  ax=g.ax_marg_y, legend = False)


megeed_dyetan = megeed_dye[(megeed_dye['IGR Type'] == 'tandem')]
megeed_dyecon = megeed_dye[(megeed_dye['IGR Type'] == 'convergent')]
megeed_dyediv = megeed_dye[(megeed_dye['IGR Type'] == 'divergent')]
megeed_dyeov = megeed_dye[(megeed_dye['IGR Type'] == 'overlap')]


def return_me_DET_per_IGR_and_model (df, the_class):
    for de in df.itertuples():
        if de[8] == the_class:
            dfsub = df[(df['DET'] == the_class)]
           # print(dfsub)
            
            
            xx = np.array(dfsub['log2_rna']).reshape((-1, 1))
            yy = np.array(dfsub['log2_Rank'])
            xx = sm.add_constant(xx)
            model = sm.OLS(yy, xx).fit()
            p_values = model.pvalues
            slope = model.params[1]
            intercept = model.params[0]
    print('P-values:', p_values)

    print(f'Equation of the regression line: y = {slope:.4f}x + {intercept:.4f}')
    print(model.summary())
    print('no of obs:', len(dfsub))
    return df


megeed_dyetan = return_me_DET_per_IGR_and_model(megeed_dyetan, 'DSonly')
megeed_dyetan = return_me_DET_per_IGR_and_model(megeed_dyetan, 'WTonly')
megeed_dyetan = return_me_DET_per_IGR_and_model(megeed_dyetan, 'bothNO')
megeed_dyetan = return_me_DET_per_IGR_and_model(megeed_dyetan, 'bothYES')

megeed_dyeov = return_me_DET_per_IGR_and_model(megeed_dyeov, 'DSonly')
megeed_dyeov = return_me_DET_per_IGR_and_model(megeed_dyeov, 'WTonly')
megeed_dyeov = return_me_DET_per_IGR_and_model(megeed_dyeov, 'bothNO')
megeed_dyeov = return_me_DET_per_IGR_and_model(megeed_dyeov, 'bothYES')


megeed_dyediv = return_me_DET_per_IGR_and_model(megeed_dyediv, 'DSonly')
megeed_dyediv = return_me_DET_per_IGR_and_model(megeed_dyediv, 'WTonly')
megeed_dyediv = return_me_DET_per_IGR_and_model(megeed_dyediv, 'bothNO')
megeed_dyediv = return_me_DET_per_IGR_and_model(megeed_dyediv, 'bothYES')

megeed_dyecon = return_me_DET_per_IGR_and_model(megeed_dyecon, 'DSonly')
megeed_dyecon = return_me_DET_per_IGR_and_model(megeed_dyecon, 'WTonly')
megeed_dyecon = return_me_DET_per_IGR_and_model(megeed_dyecon, 'bothNO')
megeed_dyecon = return_me_DET_per_IGR_and_model(megeed_dyecon, 'bothYES')

    
   

    
    

#%%


count = 0

for x in ffeat.itertuples():
    if x[15] == 'WTonly':
        count += 1
        print(count)
        

megeed_dyetan = megeed_dye[(megeed_dye['IGR Type'] == 'tandem')]
megeed_dyecon = megeed_dye[(megeed_dye['IGR Type'] == 'convergent')]
megeed_dyediv = megeed_dye[(megeed_dye['IGR Type'] == 'divergent')]
megeed_dyeov = megeed_dye[(megeed_dye['IGR Type'] == 'overlap')]

count = 0

for x in megeed_dye.itertuples():
    if x[7] == 'YES':
        count += 1
        print(count)
        
        
featwithrna  

count = 0

for x in featwithrna.itertuples():
    if x[23] == 'bothNO':
        count += 1
        print(count)      
        
        
featwithrnayYES = megeed_dyeov[(megeed_dyeov['Status'] == 'YES')]
featwithrnayNO = megeed_dyeov[(megeed_dyeov['Status'] == 'NO')]

# Assuming you have defined 'xx' and 'yy' as in your previous code
xx = np.array(featwithrnayNO['log2_rna']).reshape((-1, 1))
yy = np.array(featwithrnayNO['log2_Rank'])

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



####just for the big corre plots 
        
ffeatYES = ffeat[(ffeat['DStop3'] == 'YES')]
ffeatNO = ffeat[(ffeat['DStop3'] == 'NO')]
ffeatNO = ffeatNO.dropna(subset=['Rank_DS'])
ffeatNO

# Assuming you have defined 'xx' and 'yy' as in your previous code
xx = np.array(ffeatNO['Log2_length']).reshape((-1, 1))
yy = np.array(ffeatNO['Log2_Rank'])

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



ffeatYES = featwithrna[(featwithrna['DStop3'] == 'YES')]
ffeatNO = featwithrna[(featwithrna['DStop3'] == 'NO')]


xx = np.array(ffeatYES['Log2_length']).reshape((-1, 1))
yy = np.array(ffeatYES['Log2_Rank'])

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


######

# Data
groups = ['G1', 'G2', 'G3', 'G4']
values1 = [144, 13, 178, 51]
values2 = [3, 134, 11, 138]

####new plot, for stalls only
groups = ['WT_stall', 'SEN1_stall', 'DBL8_stall', 'DS_stall']
values1 = [13, 27, 80, 144]
values2 = [134, 120, 67, 3]

####new plot, for controls only
groups = ['WT_control', 'SEN1_control', 'DBL8_control', 'DS_control']
values1 = [51, 42, 95, 178]
values2 = [138, 147, 94, 11]


groups = ['WT_stall', 'SEN1_stall', 'DBL8_stall', 'DS_stall', 'WT_control', 'SEN1_control', 'DBL8_control', 'DS_control']
values1 = [13, 27, 80, 144, 51, 42, 95, 178]
values2 = [134, 120, 67, 3, 138, 147, 94, 11]

fig, ax = plt.subplots()

# Stacked bar chart
ax.bar(groups, values1, color = "#1D5B79", alpha = 0.8, label = "Top3 positive")
ax.bar(groups, values2, bottom = values1, color = "#EF6262", alpha = 0.8,  label = "Top3 negative")

for bar in ax.patches:
  ax.text(bar.get_x() + bar.get_width() / 2,
          bar.get_height() / 2 + bar.get_y(),
          round(bar.get_height()), ha = 'center',
          color = 'w', weight = 'bold', size = 10)

total_values = np.add(values1, values2)

# Total values labels
for i, total in enumerate(total_values):
  ax.text(i, total + 0.5, round(total),
          ha = 'center', weight = 'bold', color = 'black')


#ax.legend()
plt.show()

#%%

#%%

#%%

#%%



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


#new senataxin cc-ethanol
wtchr1, wtchr2, wtchr3, wt = Collect_data('FullMap.pombase220208-316_cc_23.txt', 10000)
sen1chr1, sen1chr2, sen1chr3, sen1 = Collect_data('FullMap.pombase220208-317_cc_23.txt', 10000)
dbl8chr1, dbl8chr2, dbl8chr3, dbl8 = Collect_data('FullMap.pombase220208-318_cc_23.txt', 10000)
dschr1, dschr2, dschr3, ds = Collect_data('FullMap.pombase220208-319_cc_23.txt', 10000)

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



gene1 = gggenes[(gggenes['chro'] == 1)]
gene2 = gggenes[(gggenes['chro'] == 2)]
gene3 = gggenes[(gggenes['chro'] == 3)]


aaaah = {'start':[3753687,1602264,1070904], 'end': [3789421,1644747,1137003], 'Chromosome': [1,2,3]}
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





#%%

def Chromosome_plot (cc1, cc2, cc3, cc4, featurex, centro, genee, telo):
    ff, (ax2, ax3, ax4, ax5, ax6) = plt.subplots(5,1, sharex=True)
   
    #ax1.set_title('WT')
    ax2.plot(cc1['Pos'], cc1['Wnorm'], color ='firebrick', alpha=0.8, linewidth = 1)
    ax2.plot(cc1['Pos'], cc1['Cnorm'], color ='steelblue', alpha=0.8, linewidth = 1)
    ax2.set_ylim(-0.2,0.2)
    ax3.set_ylim(-0.2,0.2)
    ax4.set_ylim(-0.2,0.2)
    ax5.set_ylim(-0.2,0.2)

    ax2.set_ylabel('WT (HpM)')
    #ax2.set_ylim(-5,5)
    ax2.set_ylabel('WT')

    ax3.plot(cc2['Pos'], cc2['Wnorm'], color ='firebrick', alpha=0.8, linewidth = 1)
    ax3.plot(cc2['Pos'], cc2['Cnorm'], color ='steelblue', alpha=0.8, linewidth = 1)
  #  ax3.set_ylim(-0.5,0.5)
   # ax3.set_ylabel('SEN1d (HpM)')
    ax3.set_ylabel('sen1D (HpM)')

    ax6.set_xlabel('Chromosome position')
    ax4.plot(cc3['Pos'], cc3['Wnorm'], color ='firebrick', alpha=0.8,linewidth = 1)
    ax4.plot(cc3['Pos'], cc3['Cnorm'], color ='steelblue', alpha=0.8, linewidth = 1)
  #  ax4.set_ylim(-0.5,0.5)
    ax4.set_ylabel('dbl8D (HpM)')

    ax5.plot(cc4['Pos'], cc4['Wnorm'], color ='firebrick', alpha=0.8, linewidth = 1)
    ax5.plot(cc4['Pos'], cc4['Cnorm'], color ='steelblue', alpha=0.8, linewidth = 1)
  #  ax5.set_ylim(-0.5,0.5)
    ax5.set_ylabel('sen1Ddbl8D (HpM)')
 #   ax4.set_xlabel('Chromosome position')

            
                  
    for fe in featurex.itertuples(index=False, name=None):
        if fe[5] == 'reverse':
            if fe[7] == 'sen1D':
                ax6.axvspan(fe[2],fe[3],0.3,0.5,color="teal",alpha=0.5)
               # ax6.annotate(fe[0],xy=(fe[2],0.3))
            if fe[7] == 'dbl8D':
                ax6.axvspan(fe[2],fe[3],0.3,0.5,color="orange",alpha=0.5)
                #ax6.annotate(fe[0],xy=(fe[2],0.3))
            if fe[7] == 'sen1dbl8DD_unique':
                ax6.axvspan(fe[2],fe[3],0.3,0.5,color="tomato",alpha=0.5)
                #ax6.annotate(fe[0],xy=(fe[2],0.3))
            
        elif fe[5] == 'forward':
            
            if fe[7] == 'sen1D':
                ax6.axvspan(fe[2],fe[3],0.5,0.8,color="teal",alpha=0.5)
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

    return ff


chromosome1 = Chromosome_plot(wtchr1, sen1chr1, dbl8chr1, dschr1, gene1, centro1, feat1, telo1)
chromosome2 = Chromosome_plot(wtchr2, sen1chr2, dbl8chr2, dschr2, gene2, centro2, feat2, telo2)
chromosome3 = Chromosome_plot(wtchr3, sen1chr3, dbl8chr3, dschr3, gene3, centro3, feat3, telo3)


#%%

def Chromosome_plot (cc1, cc2, cc3, cc4, featurex, centro, genee, telo):
    ff, (ax2, ax3, ax4, ax5, ax6) = plt.subplots(5,1, sharex=True)
   
    #ax1.set_title('WT')
    ax2.plot(cc1['Pos'], cc1['Wnorm'], color ='firebrick', alpha=0.8, linewidth = 1)
    ax2.plot(cc1['Pos'], cc1['Cnorm'], color ='steelblue', alpha=0.8, linewidth = 1)
    ax2.set_ylim(-0.2,0.2)
    ax3.set_ylim(-0.2,0.2)
    ax4.set_ylim(-0.2,0.2)
    ax5.set_ylim(-0.2,0.2)
    ax2.get_xaxis().set_ticks([])
    ax3.get_xaxis().set_ticks([])
    ax4.get_xaxis().set_ticks([])
    ax5.get_xaxis().set_ticks([])
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.spines['bottom'].set_visible(False)
    
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)
    ax4.spines['bottom'].set_visible(False)
    ax5.spines['top'].set_visible(False)
    ax5.spines['right'].set_visible(False)
    ax5.spines['bottom'].set_visible(False)

    ax2.set_ylabel('WT (HpM)')
    #ax2.set_ylim(-5,5)
    ax2.set_ylabel('WT')

    ax3.plot(cc2['Pos'], cc2['Wnorm'], color ='firebrick', alpha=0.8, linewidth = 1)
    ax3.plot(cc2['Pos'], cc2['Cnorm'], color ='steelblue', alpha=0.8, linewidth = 1)
  #  ax3.set_ylim(-0.5,0.5)
   # ax3.set_ylabel('SEN1d (HpM)')
    ax3.set_ylabel('sen1D (HpM)')

    ax6.set_xlabel('Chromosome position')
    ax4.plot(cc3['Pos'], cc3['Wnorm'], color ='firebrick', alpha=0.8,linewidth = 1)
    ax4.plot(cc3['Pos'], cc3['Cnorm'], color ='steelblue', alpha=0.8, linewidth = 1)
  #  ax4.set_ylim(-0.5,0.5)
    ax4.set_ylabel('dbl8D (HpM)')

    ax5.plot(cc4['Pos'], cc4['Wnorm'], color ='firebrick', alpha=0.8, linewidth = 1)
    ax5.plot(cc4['Pos'], cc4['Cnorm'], color ='steelblue', alpha=0.8, linewidth = 1)
  #  ax5.set_ylim(-0.5,0.5)
    ax5.set_ylabel('sen1Ddbl8D (HpM)')
 #   ax4.set_xlabel('Chromosome position')

            
                  
    for fe in featurex.itertuples(index=False, name=None):
        if fe[5] == 'reverse':
            if fe[7] == 'sen1D':
                ax6.axvspan(fe[2],fe[3],0.3,0.5,color="teal",alpha=0.5)
               # ax6.annotate(fe[0],xy=(fe[2],0.3))
            if fe[7] == 'dbl8D':
                ax6.axvspan(fe[2],fe[3],0.3,0.5,color="orange",alpha=0.5)
                #ax6.annotate(fe[0],xy=(fe[2],0.3))
            if fe[7] == 'sen1dbl8DD_unique':
                ax6.axvspan(fe[2],fe[3],0.3,0.5,color="tomato",alpha=0.5)
                #ax6.annotate(fe[0],xy=(fe[2],0.3))
            if fe[8] == 'NO':
                ax6.axvspan(fe[2],fe[3],0.4,0.6,color="black",alpha=0.5)
                ax6.annotate(fe[0],xy=(fe[2],0.3))
                    
            
        elif fe[5] == 'forward':
            
            if fe[7] == 'sen1D':
                ax6.axvspan(fe[2],fe[3],0.5,0.8,color="teal",alpha=0.5)
               # ax6.annotate(fe[0],xy=(fe[3],0.6))
            if fe[7] == 'dbl8D':
                ax6.axvspan(fe[2],fe[3],0.5,0.8,color="orange",alpha=0.5)
               # ax6.annotate(fe[0],xy=(fe[3],0.6))
            if fe[7] == 'sen1dbl8DD_unique':
                ax6.axvspan(fe[2],fe[3],0.5,0.8,color="tomato",alpha=0.5)
               # ax6.annotate(fe[0],xy=(fe[3],0.6))
                ax6.set_ylabel('Gene annotations')
            if fe[8] == 'NO':
                ax6.axvspan(fe[2],fe[3],0.4,0.7,color="black",alpha=0.5)
                ax6.annotate(fe[0],xy=(fe[3],0.6))
                
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


chromosome1 = Chromosome_plot(wtchr1, sen1chr1, dbl8chr1, dschr1, gene1, centro1, feat1, telo1)
chromosome2 = Chromosome_plot(wtchr2, sen1chr2, dbl8chr2, dschr2, gene2, centro2, feat2, telo2)


####I think i need to do a genome wide check, filter all, correlate + or - with transcription?


#%%


def Chromosome_plot (cc1, featurex, centro, genee, telo, join):
    ff, (ax2) = plt.subplots(1,1, sharex=True)
   
    #ax1.set_title('WT')
    ax2.plot(cc1['Pos'], cc1['Wnorm'], color ='black', alpha=0.8, linewidth = 1)
    ax2.plot(cc1['Pos'], cc1['Cnorm'], color ='grey', alpha=0.8, linewidth = 1)
    ax2.set_ylim(-0.2,0.2)

    ax2.get_xaxis().set_ticks([])

    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)

 
                
    for c in centro.itertuples(index=False, name=None):
            ax2.axvspan(c[0],c[1],0.1,0.15,color="black",alpha=0.8)
            
    for t in telo.itertuples(index=False, name=None):
            ax2.axvspan(t[0],t[1],0.1,0.15,color="grey",alpha=0.8)
            
    for ge in genee.itertuples(index=False, name=None):
        if ge[3] == 'protein coding':
            if ge[7] == 'reverse':
                ax2.axvspan(ge[4],ge[5],0.2,0.3,color="rosybrown",alpha=0.3)
            elif ge[7] == 'forward':
                ax2.axvspan(ge[4],ge[5],0.8,0.9,color="rosybrown",alpha=0.3)
                
    for j in join.itertuples(index=False, name=None):
        if j[11] == "YES":
            ax2.axvspan(j[5],j[6],0.1,0.9,color="steelblue",alpha=0.3)
        else:
            ax2.axvspan(j[5],j[6],0.1,0.9,color="firebrick",alpha=0.3)
            

    return ff


chromosome1 = Chromosome_plot(dschr1, gene1, centro1, feat1, telo1, joinedd1)
#%%



#BASH for extraction 
'''
directory="/Users/patricfernandez/Documents/seqbiasprocessed/DSstalls"



# Run fdupes to find and list duplicate files
duplicate_files=$(fdupes -r -1 "$directory")

# Iterate through the list of duplicate files and remove them
while read -r file; do
    echo "Removing duplicate file: $file"
    rm "$file"
done <<< "$duplicate_files"
'''
#%%

# Define the directory where your text files are located
#directory = '/Users/patricfernandez/Documents/seqbiasprocessed/DSstalls'

directory="/Users/patricfernandez/Documents/seqbiasprocessed/DScontrol"

directory="/Users/patricfernandez/Documents/seqbiasprocessed/WTcontrol"


directory="/Users/patricfernandez/Documents/seqbiasprocessed/WTstallsCORRECT"

directory="/Users/patricfernandez/Documents/seqbiasprocessed/DSstallsCORRECT"


#%%
#so this is the one that works
gene_status_df = pd.DataFrame(columns=['Gene Name', 'Status'])

file_dict = {}
for filename in os.listdir(directory):
    if filename.endswith('.txt'):
        # Extract the gene name from the file name
        parts = filename.split('_')
        gene_name = parts[2]  # Assuming "aat1" is the third part
        file_path = os.path.join(directory, filename)
        
        # Add the gene name and file path to the dictionary
        with open(file_path, 'r') as file:
            #file_contents = file.read()
            file_contents = pd.read_csv(file, delimiter="\t")
            #C = file_contents.split("\t")
          # print(f"File Contents:\n{file_contents}\n")

            file_dict[gene_name] = file_contents
 
        
        for key, df in file_dict.items():
          #  print(key)
            interest = df.loc[15]
            interest = interest.drop(['Position'])
          #  print(interest.max())
            if interest.loc['G'] == interest.max():
                is_largest = True
            if interest.loc['G'] != interest.max():
                    is_largest = False
                    
        
            status = 'YES' if is_largest else 'NO'
            
            temp = pd.DataFrame({'Gene Name': [key], 'Status': [status] })

           # print(temp)
        gene_status_df = pd.concat([gene_status_df, temp])
            

#%%
                
x = file_dict.get("apc2")           
feature5 = [file_dict.get('-5') for d in file_dict.dic]  
    
#%%

count = 0

for x in gene_status_df.itertuples():
    if x[2] == 'YES':
        count += 1
        print(count)
        
        
#%%
'''



# Create an empty dictionary to store the file paths
file_dict = {}

# List all the files in the directory
for filename in os.listdir(directory):
    if filename.endswith('.txt'):
        # Extract the gene name from the file name
        parts = filename.split('_')
        gene_name = parts[2]  # Assuming "aat1" is the third part
        file_path = os.path.join(directory, filename)
        
        # Add the gene name and file path to the dictionary
        file_dict[gene_name] = file_path

# Now, you have a dictionary where the keys are gene names and the values are file paths
# You can iterate through the dictionary based on gene name
for gene_name, file_path in file_dict.items():
    print(f"Gene Name: {gene_name}")
    print(f"File Path: {file_path}")
    with open(file_path, 'r') as file:
        file_contents = file.read()
        print(f"File Contents:\n{file_contents}\n")
        
        
        
        
        
        
        

        with open(file_path, 'r') as file:
            file_contents = file.readlines()
            # Split the header line (column names) and data rows
            header = file_contents[0].strip().split('\t')
            data = file_contents[1:]

            # Find the index of '-5' in the header
            column_index = header.index('-5')

            # Check if 'G' at '-5' is the largest value in the row
            is_largest = True
            for row in data:
                values = row.strip().split('\t')
                if float(values[column_index]) < max(float(val) for val in values):
                    is_largest = False
                    break

            # Add the gene name and status to the DataFrame
            status = 'YES' if is_largest else 'NO'
            gene_status_df = gene_status_df.append({'Gene Name': gene_name, 'Status': status}, ignore_index=True)

# Now, you have a dictionary where the keys are gene names and the values are file paths
# You also have a DataFrame with gene names and their 'YES' or 'NO' status
print(gene_status_df)

#%%

#new version to ,mine 


# Define the directory where your text files are located
directory = '/Users/patricfernandez/Documents/seqbiasprocessed/DSstalls'

# Create an empty dictionary to store the file paths
file_dict = {}

# Create an empty DataFrame to store the gene name and 'YES' or 'NO'
gene_status_df = pd.DataFrame(columns=['Gene Name', 'Status'])

# List all the files in the directory
for filename in os.listdir(directory):
    if filename.endswith('.txt'):
        # Extract the gene name from the file name
        parts = filename.split('_')
        gene_name = parts[2]  # Assuming "aat1" is the third part
        file_path = os.path.join(directory, filename)
        
        # Add the gene name and file path to the dictionary
        file_dict[gene_name] = file_path

        # Read and analyze the contents of the text file
        with open(file_path, 'r') as file:
            file_contents = file.readlines()
            # Split the header line (column names) and data rows
            header = file_contents[0].strip().split('\t')
            data = file_contents[1:]

            # Find the index of '-5' in the header
            column_index = header.index('-5')

            # Check if 'G' at '-5' is the largest value in the row
            is_largest = True
            for row in data:
                values = row.strip().split('\t')
                if float(values[column_index]) < max(float(val) for val in values):
                    is_largest = False
                    break

            # Add the gene name and status to the DataFrame
            status = 'YES' if is_largest else 'NO'
            gene_status_df = gene_status_df.append({'Gene Name': gene_name, 'Status': status}, ignore_index=True)

# Now, you have a dictionary where the keys are gene names and the values are file paths
# You also have a DataFrame with gene names and their 'YES' or 'NO' status
print(gene_status_df)

'''


