#!-*-coding=utf-8-*-
##---------------------------------
#		DATE:208-12-07
#		AUTHOR:Lijie Zhang
#		Function: Convert the Ensemble ID to Gene name using RNA-seq data in TCGA, depend on gene2ensembl.gz file and Homo_sapiens.gene_info.gz file in NCBI 
#
##---------------------------------

import sys,os
import requests
import json
import glob
import gzip

#gene2ensembl
def get_gene2ensembl_dict(gene2ensembl_file):
	f=gzip.GzipFile(gene2ensembl_file, "r")
	lines=f.readlines()
	f.close()
	result_dict={}
	for line in lines[1:]:
		line_data = line.split("\t")
		if "ENSG" in line_data[2]:
			if line_data[2].strip() not in result_dict.keys():
				result_dict[line_data[2].strip()]=line_data[1].strip() 
	return result_dict

#Homo_sapiens.gene_info.gz
def get_gene_info_dict(homo_sapiens_file,gene2ensembl_dict):
	f=gzip.GzipFile(homo_sapiens_file,"r")
	lines=f.readlines()
	f.close()
	result_dict={}
	for line in lines[1:]:
		line_data=line.split("\t")
		result_dict[line_data[1].strip()]=line_data[2].strip()

	for key,value in gene2ensembl_dict.items():
		if value in result_dict.keys():
			gene2ensembl_dict[key]=result_dict[value] 
	return gene2ensembl_dict


#recuve all htseq count files and store in a directory
def deal_all_htseq_files(input_dir,gene_ensembl_dict,file_suffix="*.htseq.counts.gz",store_dir="htseq"):
	path=os.path.join(input_dir,file_suffix)
	if not os.path.exists(store_dir):
		os.mkdirs(store_dir)
	all_files = glob.glob(path)
	for htseq_file in all_files:
		filename=os.path.basename(htseq_file)
		temp_path=os.path.join(store_dir,filename[:-3]+".txt")
		s=open(temp_path,"w")
		f=gzip.GzipFile(homo_sapiens_file,"r")
		line=f.readline()
		while line.startswith("ENSG"):
			line_data=line.split("\t")
			ensg = line_data[0][0:15]
			count = line_data[1].strip()
			if ensg in gene_ensembl_dict.keys():
				ensg = gene_ensembl_dict[ensg]
			s.write("%s\t%s\n"%(ensg,count))
			line=f.readline() 
		s.close()
	return all_files

	


  


if __name__=="__main__":
	#list_all_files("/home/zhanglj/skin_cancer/RNA-seq/*") 
	dict1=get_gene2ensembl_dict("/home/zhanglj/skin_cancer/gene2ensembl.gz")         
	gene_ensembl_dict=get_gene_info_dict("/home/zhanglj/skin_cancer/Homo_sapiens.gene_info.gz",dict1)  
	deal_all_htseq_files(input_dir,gene_ensembl_dict,file_suffix="*.htseq.counts.gz",store_dir="/home/zhanglj/skin_cancer/htseq")   
    
