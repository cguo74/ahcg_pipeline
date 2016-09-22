#!/usr/bin/env python
#coding:utf-8
need = []
flist = open("list.txt", 'r')
for fl in flist:
	a = fl.split('.')[0].split(' ')[1]
	need.append(a)
	print need

ff = open ('hg19_refGene.txt','r')
for f in ff:
	for n in need:
		if n in f:
			#print f
			#print f.split('\t')[9:11]
			aa = f.split('\t')
			c4 = aa[1]
			c5 = aa[3]
			c2 = aa[9].strip(',')
			c3 = aa[10].strip(',')
			c1 = aa[2]
			len_c2 = len(c2.split(','))

			#len_c3 = len(c3.split(','))
			
			c22 = c2.split(',')
			c33 = c3.split(',')
			
			for p in range(len_c2):
				print c1, c22[p], c33[p], c4, c5 
				

