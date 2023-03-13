import numpy as np
import os
from Bio import SeqIO
import pickle
def countkmer(my_dir,name):
    
    a_mer=''
    m=['A','T','G','C']
    my_mer=[]
    for o in range(len(m)):
        for p in range(len(m)):
            for q in range(len(m)):
                for r in range(len(m)):
                    my_item = str(m[o])+str(m[p])+str(m[q]) +str(m[r])
                    my_mer.append(my_item)
                    a_mer += str(my_item)+','

    kmer_dic={}
    
    
    contig_n = 0

    for record in SeqIO.parse(my_dir + name, "fasta") :


        seq=record.seq

        if len(seq) >= 9000:
            contig_n += 1
            num=[]
            for n  in range(0,len(seq)):
                num.append(n)
                
            every_1800=num[0:len(seq):1800]

            for i in every_1800:
                my_num=[]
                dic_kmer={}
                for item in my_mer:
                    dic_kmer[item]=0
                if i+9000 <= len(seq):
                    for j in range(i,i+9000):
                        if j+4 <= i+9000:
                            kmer = seq[j:j+4]
                            if kmer in dic_kmer:
                                dic_kmer[kmer] +=1
                    for item_n in my_mer:
                        my_num.append(int(dic_kmer[item_n]) )
                    k_num = str(int(i/1800))
                    kmer_dic[str(contig_n) + '+' + k_num] = my_num

            my_num=[]
            dic_kmer={}
            for itemk in my_mer:
                dic_kmer[itemk]=0
            for k in range(len(seq)-9000,len(seq)):
                if  k+4 <= len(seq):
                    kmerk = seq[k:k+4]
                    if kmerk in dic_kmer.keys():
                        dic_kmer[kmerk] +=1 
            for item_nk in my_mer:
                my_num.append(int(dic_kmer[item_nk]))
            lne= int(len(seq)-9000)
            k_num = str(int(lne//1800 + 1))
            kmer_dic[str(contig_n) + '+'  + k_num] = my_num

    
    if contig_n == 0:
        dic_kmer={}
        for item3 in my_mer:
            dic_kmer[item3]=0
        my_num=[]
        max_contig=''
        for record in SeqIO.parse(my_dir + name, "fasta") :
            seq=record.seq
            if len(seq) > len(max_contig):
                max_contig = seq
        for n in range(len(seq)):
            if n+4 <= len(seq):
                kmer3 = seq[n:n+4]
                if kmer3 in dic_kmer.keys():
                    dic_kmer[kmer3] +=1
        for item_n3 in my_mer:
            my_num .append(int(dic_kmer[item_n3]))
        kmer_dic['1'+ '+' + '0']=my_num

    return(kmer_dic)






def pred(dic_phage,dic_bac):
    
  
    corr_dic={}
    p_lst=[]
    b_lst=[]
    
    for p_key in dic_phage.keys():
        p_lst.append(p_key)

    for b_key in dic_bac.keys():
        b_lst.append(b_key)

    my_corr={}

    for m in range(0,len(p_lst)):
        for n in range(0,len(b_lst)):

            x = dic_phage[p_lst[m]]
            y = dic_bac[b_lst[n]]


            my_corr[p_lst[m] + '_' + b_lst[n]]=np.corrcoef(x,y)[0][-1]


    max_cor=-1
    for test_key in my_corr.keys():
        if my_corr[test_key]>max_cor:
            max_cor=my_corr[test_key]
            max_lne=test_key       

    #print(max_lne)

    phage_lne=str(max_lne.split('_')[0])
    bac_lne=str(max_lne.split('_')[1])
 
    pb_sub=[]
    for r in range(0,256):
        every_sub= dic_phage[phage_lne][r] - dic_bac[bac_lne][r]
        pb_sub.append(every_sub)

    test_x=[]
    test_x.append(pb_sub)
    input_model = open('./PBLKS_package/PBLKS_model.pkl','rb')
    clf = pickle.load(input_model)

    pred_y=clf.predict(test_x)

    if str(pred_y) == '[1]':
        print('The prediction result of PB-LKS is: 1.(The queried bacteria is the host of the queried phage.)')
    if str(pred_y) == '[0]':
        print('The prediction result of PB-LKS is: 0.(The queried bacteria is not the host of the queried phage.)')



