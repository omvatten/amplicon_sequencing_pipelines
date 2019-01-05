#Python3 script to join frequency tables from samples processed separately using USEARCH
import os
import pandas as pd

#Returns list of names and sequences of SVs from fasta
def listofseqs(f):
    with open(f) as r:
        read = r.readlines()
    namesvs = [[],[]]
    
    for i in range(len(read)):
        if i==0:
            name = read[i].strip()[1:] 
            seq=''
        elif i==len(read)-1:
            seq=seq+read[i].strip()
            namesvs[0].append(name)
            namesvs[1].append(seq)
        elif read[i][0] != '>':
            seq = seq+read[i].strip()
        elif read[i][0] == '>':
            namesvs[0].append(name)
            namesvs[1].append(seq)
            name = read[i].strip()[1:] 
            seq=''
    return namesvs

#Make list of files names for tables and fasta
path = '../out/'
tabfiles = [f for f in os.listdir(path) if f.startswith('tab')]
fastafiles = [f for f in os.listdir(path) if f.startswith('otu')]

#Make dictionary holding fasta lists and tables dataframes
seqdicts = {}
tabdicts = {}
for t in tabfiles:
    tabname = t.split('_')[0].split('.')[1]
    tabdicts[tabname] = pd.read_csv(path+t, sep='\t', index_col=0)
for s in fastafiles:
    seqname = s.split('_')[0].split('.')[1]
    seqdicts[seqname] = listofseqs(path+s)

#Associate each Otu name in a table with a sequence
for smp in tabdicts.keys():
    df = tabdicts[smp]
    slist = seqdicts[smp]

    df['seq']=0
    for i in df.index:
        seqtoadd = slist[1][slist[0].index(i)]
        df.loc[i,'seq']=seqtoadd #Add the actual sequence to the df

    tabdicts[smp]=df

#Merge all tables, index is sequence
listofsamples = list(tabdicts.keys())
for snr in range(1,len(listofsamples)):
    if snr==1:
        df1 = tabdicts[listofsamples[snr-1]]
        df2 = tabdicts[listofsamples[snr]]
        df = pd.merge(df1, df2, on='seq', how='outer')
    else:
        df = pd.merge(df, tabdicts[listofsamples[snr]], on='seq', how='outer')
df = df.set_index('seq')
df = df.fillna(0)

#Give new names to sequences and return as fasta file and table
listofseqs=list(df.index)
svlist=[]
for i in range(len(listofseqs)):
    svlist.append('SV'+str(i+1))

df['sv']=svlist
df = df.set_index('sv')
df = df.applymap(int)

with open(path+'All_seqs_nonpooled.fa', 'w') as f:
    for i in range(len(listofseqs)):
        f.write('>'+svlist[i]+'\n')
        f.write(listofseqs[i]+'\n')
        
df.to_csv(path+'All_tab_nonpooled.txt', sep='\t')

print('Finished joining tabs')
