import csv
import numpy as np

pathDtab='../path_to_folder_with_DADA2.R_output/'
Dfile='DADA2_frequency_table.txt'
outtabname='path/converted_frequency_table.txt'
outfastaname='path/converted_fasta_file.fa'

dd2table=[]
with open(pathDtab+Dfile, 'r', newline='') as f:
    rd=csv.reader(f,delimiter='\t')
    for line in rd:
        dd2table.append(line)

samplenames=[]
for i in range(1,len(dd2table)):
    samplenames.append(dd2table[i][0])

ASVseqs=[]
for i in range(len(dd2table[0])):
    ASVseqs.append(dd2table[0][i])

valRaw=np.zeros((len(ASVseqs),len(samplenames)))
for smp in range(len(samplenames)):
    for asv in range(len(ASVseqs)):
        valRaw[asv,smp]=float(dd2table[smp+1][asv+1])

dtable_out=[['ASV_ID']+samplenames]
asvfile_out=[]

for asv in range(len(ASVseqs)):
    dtable_out.append(['ASV'+str(asv+1)]+valRaw[asv].tolist())
    asvfile_out.append('>ASV'+str(asv+1))
    asvfile_out.append(ASVseqs[asv])

for i in range(1,len(dtable_out)):
    for j in range(1,len(dtable_out[0])):
        dtable_out[i][j]=int(dtable_out[i][j])

with open(outfastaname,'w') as f:
    for line in asvfile_out:
        f.write('{}\n'.format(line))

with open(outtabname,'w',newline='') as f:
    wr=csv.writer(f,delimiter='\t')
    wr.writerows(dtable_out)

print("Finished")
