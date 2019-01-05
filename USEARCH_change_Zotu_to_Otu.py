#Python3 script to change the name of the sequences in zotus.fa files from zotu to otu
import os

files=os.listdir('../out')

for fil in files:
	if fil[0]=='z':
		
		with open('../out/'+fil,'r') as f:
			readfile=f.readlines()

		count=0
		for i in range(len(readfile)):
			if readfile[i][0]=='>':
				count+=1
				readfile[i]='>Otu'+str(count)+'\n'
		print(fil,'  ',count)

		with open('../out/'+fil,'w',newline='') as f:
			for i in readfile:
				f.write(i)
print('Finished')
