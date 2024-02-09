#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 16:01:16 2023

@author: patricfernandez
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns



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


#averaged
wtETchr1, wtETchr2, wtETchr3, wtET = Collect_data('FullMap.pombase220208-WT.txt', 10000)
t2Echr1, t2Echr2, t2Echr3, t2E = Collect_data('FullMap.pombase220208-T2.txt',10000)

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

#%%

#binned correlation data


normscore_df1 = wtET['score']
normscore_df2 = wtETB['score']

normscore_df3 = t2E['score']
normscore_df4 = t2EB['score']

# Create a new DataFrame with the selected 'normscore' columns
new_df = pd.DataFrame({'WT_A': normscore_df1, 'WT_B': normscore_df2, 'T2_A': normscore_df3, 'T2_B': normscore_df4})

iris_corr_matrix = new_df.corr()
print(iris_corr_matrix)

# Create the heatmap using the `heatmap` function of Seaborn
sns.heatmap(iris_corr_matrix, cmap='coolwarm', annot=True)

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

#need to bare in mind that in order to call IGRs i need to format a df of genes with expression data
#first take all features data


Genesy = pd.read_csv('Schizosaccharomyces_pombe_all_chromosomes.gff3',sep='\t', names=('1', '2','3', '4','5', '6','7', '8', '9'),skiprows=1)
unique_items_array = Genesy['3'].unique()
print(unique_items_array)

#make a happy mask to find only genes (original for all in CC_ethanol)
mask = Genesy['3'] == 'gene'
Genes = Genesy[mask]

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





#now retrieve cpc RNA data (published Margret)

xls = pd.ExcelFile('S_pombe_expression_data_PMID_23101633.xlsx')
df1 = pd.read_excel(xls, 'Table_S4')
df1.columns=df1.iloc[5]

df1 = df1.iloc[6: , :]

rob = df1[["Systematic.name", "MM.mRNA.cpc"]]


rob.rename({'Systematic.name': 'Systematic ID'}, axis=1, inplace=True)
rob.rename({"MM.mRNA.cpc": 'rna'}, axis=1, inplace=True)


GenesandRob = rob.merge(Genes, how='left', on='Systematic ID')
GenesandRoby = GenesandRob.dropna(how='any', axis = 0)
joinedd1 = GenesandRoby[(GenesandRoby['Chr'] == 1)]
joinedd2 = GenesandRoby[(GenesandRoby['Chr'] == 2)]
joinedd3 = GenesandRoby[(GenesandRoby['Chr'] == 3)]


#%%
#function for id of IGR
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



#%%

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





#%%








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

'''
sns.catplot(data=melted_df, x="IGR Type", y="log2_score", kind="box", hue= 'source')

sns.jointplot(
    data=dataIGR,
    x="log2_length", y="log2_epb", hue="IGR Type",
    kind = 'kde'
)


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

'''


#%%%%

#now using combined GenesandRob which is expression data + gene annotated 

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

#################
merged_feat['length_quartile'] = pd.qcut(merged_feat['log2_length'], q=4, labels=False)
merged_feat['epb_quartile'] = pd.qcut(merged_feat['log2_epb'], q=4, labels=False)
merged_feat['expression_quartile'] = pd.qcut(merged_feat['rna'], q=4, labels=False)

merged_feat.loc[merged_feat['expression_quartile'] == 0, 'expression_quartile'] = 'Q1'
merged_feat.loc[merged_feat['expression_quartile'] == 1, 'expression_quartile'] = 'Q2'
merged_feat.loc[merged_feat['expression_quartile'] == 2, 'expression_quartile'] = 'Q3'
merged_feat.loc[merged_feat['expression_quartile'] == 3, 'expression_quartile'] = 'Q4'







df_cleanedg = merged_feat.dropna(subset=['rna','log2_rna', 'log2_epb'])

df_cleanedg['rna'] = pd.to_numeric(df_cleanedg['rna'], errors='coerce')

# Check for NaN values in the 'rna' column
nan_count = df_cleanedg['rna'].isna().sum()

# If there are NaN values, drop them and proceed with qcut
if nan_count > 0:
    df_cleanedg = df_cleanedg.dropna(subset=['rna'])
    df_cleanedg['expression_quartiley'] = pd.qcut(df_cleanedg['rna'], q=4, labels=False)
else:
    # No NaN values, proceed with qcut directly
    df_cleanedg['expression_quartiley'] = pd.qcut(df_cleanedg['rna'], q=4, labels=False)




df_cleanedg.loc[df_cleanedg['expression_quartiley'] == 0, 'expression_quartiley'] = 'Q1'
df_cleanedg.loc[df_cleanedg['expression_quartiley'] == 1, 'expression_quartiley'] = 'Q2'
df_cleanedg.loc[df_cleanedg['expression_quartiley'] == 2, 'expression_quartiley'] = 'Q3'
df_cleanedg.loc[df_cleanedg['expression_quartiley'] == 3, 'expression_quartiley'] = 'Q4'


df_cleanedgy = df_cleanedg.copy()
df_cleanedgy["expression_quartiley"] = pd.Categorical(df_cleanedgy["expression_quartiley"],
                                                 ordered = True,
                                                 categories = ["Q1", "Q2", "Q3", "Q4"])
#the plot

g = sns.JointGrid()
x, y = df_cleanedgy["log2_rna"], df_cleanedgy["log2_epb"]
sns.regplot(x=x, y=y, scatter_kws={'s': 1}, line_kws={'color': 'black'}, ax=g.ax_joint)
sns.histplot(data=df_cleanedgy, x="log2_rna", hue="expression_quartiley",
    element="step",  common_norm = False, ax=g.ax_marg_x, legend = False)
sns.histplot(data=df_cleanedgy, y="log2_epb", hue="expression_quartiley",
    element="step",  common_norm = False,  ax=g.ax_marg_y, legend = False)

plt.ylim(-7.5,-4.5)

#%%

from sklearn.linear_model import LinearRegression



xx = np.array(df_cleanedgy['log2_rna']).reshape((-1, 1))
yy = np.array(df_cleanedgy['log2_epb'])
modell = LinearRegression().fit(xx, yy)
r_sqC = modell.score(xx, yy)
print('coefficient of determination:', r_sqC)



# Assuming you have defined 'xx' and 'yy' as in your previous code
modell = LinearRegression().fit(xx, yy)

# Get the slope (coefficient) and intercept from the model
slope = modell.coef_[0]
intercept = modell.intercept_

# Print the equation of the regression line
print(f'Equation of the regression line: y = {slope:.4f}x + {intercept:.4f}')






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
    data=df_cleanedgy,
    x="log2_rna", y="log2_epb", hue="expression_quartiley"
    
)

sns.jointplot(
    data=df_cleanedgy,
    x="log2_rna", y="log2_epb", hue="expression_quartiley",
    ylim= [-8,-5]
)

sns.lmplot(data=df_cleanedgy, x="log2_rna", y="log2_epb", hue="expression_quartiley")




sns.lmplot(
    data=df_cleanedgy, x="log2_rna", y="log2_epb",
    col="expression_quartiley", height=3,
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
custom_palette = ["#1D5B79", "#468B97", "#EF6262", "#F3AA60"]  # Replace these colors with your desired ones

# Set the custom color palette
sns.set_palette(custom_palette)

df_cleaned_tandem = df_cleaned[(df_cleaned['IGR Type'] == 'tandem')]
df_cleaned_convergent = df_cleaned[(df_cleaned['IGR Type'] == 'divergent')]
df_cleaned_divergent = df_cleaned[(df_cleaned['IGR Type'] == 'convergent')]

#divergent IGRs only

g = sns.JointGrid()
x, y = df_cleaned_divergent["log2_rna"], df_cleaned_divergent["log2_epb"]
sns.regplot(x=x, y=y, scatter_kws={'s': 1}, line_kws={'color': 'black'}, ax=g.ax_joint)
sns.histplot(data=df_cleaned_divergent, x="log2_rna", hue="expression_quartile",
    element="step",  common_norm = False, ax=g.ax_marg_x, legend = False)
sns.histplot(data=df_cleaned_divergent, y="log2_epb", hue="expression_quartile",
    element="step",  common_norm = False,  ax=g.ax_marg_y, legend = False)
plt.ylim(-11,-3)

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


from sklearn.linear_model import LinearRegression

xx = np.array(df_cleaned_divergent['log2_rna']).reshape((-1,1))
yy = np.array(df_cleaned_divergent['log2_epb'])
#model = LinearRegression()
#model.fit(x, y)
model = LinearRegression().fit(xx, yy)
r_sq = model.score(xx, yy)
print('coefficient of determination:', r_sq)

#####



#convergent IGRs only
g = sns.JointGrid()
x, y = df_cleaned_convergent["log2_rna"], df_cleaned_convergent["log2_epb"]
sns.regplot(x=x, y=y, scatter_kws={'s': 1}, line_kws={'color': 'black'}, ax=g.ax_joint)
sns.histplot(data=df_cleaned_convergent, x="log2_rna", hue="expression_quartile",
    element="step",  common_norm = False, ax=g.ax_marg_x, legend = False)
sns.histplot(data=df_cleaned_convergent, y="log2_epb", hue="expression_quartile",
    element="step",  common_norm = False,  ax=g.ax_marg_y, legend = False)
plt.ylim(-11,-3)

xx = np.array(df_cleaned_convergent['log2_rna']).reshape((-1, 1))
yy = np.array(df_cleaned_convergent['log2_epb'])

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


xx = np.array(df_cleaned_convergent['log2_rna']).reshape((-1,1))
yy = np.array(df_cleaned_convergent['log2_epb'])
#model = LinearRegression()
#model.fit(x, y)
model = LinearRegression().fit(xx, yy)
r_sq = model.score(xx, yy)
print('coefficient of determination:', r_sq)




#####
#tandem IGRs only
g = sns.JointGrid()
x, y = df_cleaned_tandem["log2_rna"], df_cleaned_tandem["log2_epb"]
sns.regplot(x=x, y=y, scatter_kws={'s': 1}, line_kws={'color': 'black'}, ax=g.ax_joint)
sns.histplot(data=df_cleaned_tandem, x="log2_rna", hue="expression_quartile",
    element="step",  common_norm = False, ax=g.ax_marg_x, legend = False)
sns.histplot(data=df_cleaned_tandem, y="log2_epb", hue="expression_quartile",
    element="step",  common_norm = False,  ax=g.ax_marg_y, legend = False)
plt.ylim(-11,-3)

xx = np.array(df_cleaned_tandem['log2_rna']).reshape((-1, 1))
yy = np.array(df_cleaned_tandem['log2_epb'])

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


xx = np.array(df_cleaned_tandem['log2_rna']).reshape((-1,1))
yy = np.array(df_cleaned_tandem['log2_epb'])
#model = LinearRegression()
#model.fit(x, y)
model = LinearRegression().fit(xx, yy)
r_sq = model.score(xx, yy)
print('coefficient of determination:', r_sq)









#all IGRs
g = sns.JointGrid()
x, y = df_cleaned["log2_rna"], df_cleaned["log2_epb"]
sns.regplot(x=x, y=y, scatter_kws={'s': 1}, line_kws={'color': 'black'}, ax=g.ax_joint)
sns.histplot(data=df_cleaned, x="log2_rna", hue="expression_quartile",
    element="step",  common_norm = False, ax=g.ax_marg_x, legend = False)
sns.histplot(data=df_cleaned, y="log2_epb", hue="expression_quartile",
    element="step",  common_norm = False,  ax=g.ax_marg_y, legend = False)

plt.ylim(-7.5,-4.5)


xx = np.array(df_cleaned['log2_rna']).reshape((-1, 1))
yy = np.array(df_cleaned['log2_epb'])

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


xx = np.array(df_cleaned['log2_rna']).reshape((-1,1))
yy = np.array(df_cleaned['log2_epb'])
#model = LinearRegression()
#model.fit(x, y)
model = LinearRegression().fit(xx, yy)
r_sq = model.score(xx, yy)
print('coefficient of determination:', r_sq)


from sklearn.linear_model import LinearRegression

x = np.array(df_cleaned['log2_rna']).reshape((-1,1))
y = np.array(df_cleaned['log2_epb'])
#model = LinearRegression()
#model.fit(x, y)
model = LinearRegression().fit(x, y)
r_sq = model.score(x, y)
print('coefficient of determination:', r_sq)



#all protein 
g = sns.JointGrid()
x, y = df_cleanedgy["log2_rna"], df_cleanedgy["log2_epb"]
sns.regplot(x=x, y=y, scatter_kws={'s': 1}, line_kws={'color': 'black'}, ax=g.ax_joint)
sns.histplot(data=df_cleanedgy, x="log2_rna", hue="expression_quartiley",
    element="step",  common_norm = False, ax=g.ax_marg_x, legend = False)
sns.histplot(data=df_cleanedgy, y="log2_epb", hue="expression_quartiley",
    element="step",  common_norm = False,  ax=g.ax_marg_y, legend = False)

#plt.ylim(-7.5,-4.5)

xx = np.array(df_cleanedgy['log2_rna']).reshape((-1, 1))
yy = np.array(df_cleanedgy['log2_epb'])

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


xx = np.array(df_cleaned_tandem['log2_rna']).reshape((-1,1))
yy = np.array(df_cleaned_tandem['log2_epb'])
#model = LinearRegression()
#model.fit(x, y)
model = LinearRegression().fit(xx, yy)
r_sq = model.score(xx, yy)
print('coefficient of determination:', r_sq)


###################
####################
#####################
####################
##################

#so df_cleaned (IGR) is corrrect and matches plot
#and df_cleanedgy (proteins) also correct and matched 
#therefore do ECDF on these 
#%%

#protein coding by gene only 
sns.lmplot(
    data=df_cleanedgy, x="log2_rna", y="log2_epb",
    col="expression_quartiley", height=3, hue ="expression_quartiley",
    facet_kws=dict(sharex=False, sharey=False),scatter_kws={'s': 1}, line_kws={'color': 'black'}
)


proteinq1 = df_cleanedgy[(df_cleanedgy['expression_quartiley'] == 'Q1')]
proteinq2 = df_cleanedgy[(df_cleanedgy['expression_quartiley'] == 'Q2')]
proteinq3 = df_cleanedgy[(df_cleanedgy['expression_quartiley'] == 'Q3')]
proteinq4 = df_cleanedgy[(df_cleanedgy['expression_quartiley'] == 'Q4')]

xx = np.array(proteinq4['log2_rna']).reshape((-1, 1))
yy = np.array(proteinq4['log2_epb'])

xx = sm.add_constant(xx)
model = sm.OLS(yy, xx).fit()
p_values = model.pvalues
print('P-values:', p_values)
slope = model.params[1]
intercept = model.params[0]
print(f'Equation of the regression line: y = {slope:.4f}x + {intercept:.4f}')





#IGRs only 
sns.lmplot(
    data=df_cleaned, x="log2_rna", y="log2_epb",
    col="expression_quartile", height=3, hue ="expression_quartile",
    facet_kws=dict(sharex=False, sharey=False),scatter_kws={'s': 1}, line_kws={'color': 'black'}
)


igrq1 = df_cleaned[(df_cleaned['expression_quartile'] == 'Q1')]
igrq2 = df_cleaned[(df_cleaned['expression_quartile'] == 'Q2')]
igrq3 = df_cleaned[(df_cleaned['expression_quartile'] == 'Q3')]
igrq4 = df_cleaned[(df_cleaned['expression_quartile'] == 'Q4')]

xx = np.array(igrq4['log2_rna']).reshape((-1, 1))
yy = np.array(igrq4['log2_epb'])

xx = sm.add_constant(xx)
model = sm.OLS(yy, xx).fit()
p_values = model.pvalues
print('P-values:', p_values)
slope = model.params[1]
intercept = model.params[0]
print(f'Equation of the regression line: y = {slope:.4f}x + {intercept:.4f}')






#%%



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




sns.catplot(data=intrainterr, x="IGR Type", y="log2_score", kind="box", hue= 'source', showfliers=False)

dataIGR.rename(columns={'IGR name': 'Systematic ID'}, inplace=True)









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

'''
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
'''

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

#%%
from scipy.stats import ks_2samp




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

####

merged_featy = merged_feat.copy()
merged_featy["expression_quartile"] = pd.Categorical(merged_featy["expression_quartile"],
                                                 ordered = True,
                                                 categories = ["Q1", "Q2", "Q3", "Q4"])
custom_palette = ["#1D5B79", "#468B97", "#EF6262", "#F3AA60"]  # Replace these colors with your desired ones
sns.set_palette(custom_palette)


