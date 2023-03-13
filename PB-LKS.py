import os
import sys
from PBLKS_package.PBLKS_pred import countkmer,pred

phage=input('please enter a file name of phage:')
phage_item=phage.split('/')[0:-1]
phage_dir=''
for i in range(0,len(phage_item)):
    phage_dir += phage_item[i] + '/'
phage_name=phage.split('/')[-1]

if phage_name not in os.listdir(phage_dir):
    print('There is no corresponding file in the folder')
    exit(0)
    
bac=input('please enter a file name of bacteria:')
bac_item=bac.split('/')[0:-1]
bac_dir=''
for j in range(0,len(bac_item)):
    bac_dir += bac_item[j] + '/'
bac_name=bac.split('/')[-1]

if bac_name not in os.listdir(bac_dir):
    print('There is no corresponding file in the folder')
    exit(0)



dic_phage=countkmer(phage_dir,phage_name)
dic_bac=countkmer(bac_dir,bac_name)

pred_result=pred(dic_phage,dic_bac)

print('The prediction is over.')

