#preprocessing

import sys
sys.path.insert(0, './PyUtils')
sys.path.insert(0, './PyGEO')

import os
import numpy
import getopt
from bcolors import *

def writeMoses(Dataobj, filename, delim):
		outfile = open(filename,'w')
		outfile.write(delim.join(list(Dataobj[0]))+'\n')
		for x in Dataobj[1:]:
			outfile.write(str(int(x[0]))+delim+delim.join(x[1:])+'\n')
		outfile.close()
		print 'PREPROCESSING COMPLETE, CREATED MOSES FILE '+filename+'\n'

def writeSVM(Dataobj, filename):
		outfile = open(filename,'w')
		#FEATUREINDEX = range(1:len(Dataobj[0]))
		DATA = Dataobj[1:]
		for x in DATA:
			outfile.write(x[0])
			for i in range(1,len(x[1:])):
				outfile.write(' '+repr(i)+':'+x[i])
			outfile.write('\n')
		outfile.close()
		print 'PREPROCESSING COMPLETE, CREATED SVM FILE '+filename+'\n'

def preprocess(argv):
	print bcolors.BOLD + bcolors.HEADER +'PREPROCESSING' + bcolors.ENDC +'\n'
	outdir = 'out'
	itype = ''
	binarize = 1
	ifile = ''
	ofile = ''
	otype = 'moses'
	p_cutoff = 0
	targetCategory = ''
	targetClass=''
	idf = 'IDENTIFIER'
	try:
		opts, args = getopt.getopt(argv,"hi:o:d:p:b:t:c:y:",["ifile=","ofile=","identifier=","pval=","binarize=","targetCategory=","targetClass=","outputType="])
	except getopt.GetoptError:
		print 'test.py -i <inputfile> -o <outputfile>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'test.py -i <inputfile> -o <outputfile>'
			sys.exit()
		elif opt in ("-i", "--ifile"):
			supported_itype = {	'dat'	:'Generic flatfile',
						'data'	:'Generic flatfile',
						'txt'	:'Generic flatfile',
						'soft'	:'NCBI microarray geo soft format',
						'moses'	:'MOSES learning set format',
						'svm'	:'libsvm specified svm format'
					}
			ifile = arg
			itype = arg.split('.')[-1]
			print 'INPUT FILE NAME: '
			print '\t' + ifile.split('/')[-1]+'\n'
			if not os.path.isfile(ifile):
				print bcolors.FAIL + bcolors.BOLD + '[ERROR]' + bcolors.ENDC + ' SPECIFIED INPUT FILE \'.'+ifile+'\' CANNOT BE READ'
				exit()
			if itype not in supported_itype:
				print bcolors.FAIL + bcolors.BOLD + '[ERROR]' + bcolors.ENDC + ' SPECIFIED INPUT TYPE \'.'+itype+'\' NOT SUPPORTED'
				exit()
			else:
				print 'INFERRED INPUT TYPE:'
				print '\t' + itype + ': '+supported_itype[itype]+'\n'		
		elif opt in ("-o", "--ofile"):
			ofile = arg
		elif opt in ("-d", "--identifier"):
			idf = arg
		elif opt in ("-p", "--pval"):
			p_cutoff = arg
		elif opt in ("-t", "targetCategory="):
			targetCategory = arg
		elif opt in ("-c", "targetClass="):
			targetClass = arg
		elif opt in ("-b", "--binarize"):
			binarize = arg
		elif opt in ("-y", "--outputType"):
			supported_otype = ['moses','svm']
			otype = arg
			if otype not in supported_otype:
				print bcolors.FAIL + bcolors.BOLD + '[ERROR]' + bcolors.ENDC +' SPECIFIED OUTPUT TYPE \''+otype+'\' NOT SUPPORTED, EXECUTION HALTED'
				exit()
	fileStem = '.'.join(ifile.split('/')[-1].split('.')[0:-1])
	if ofile == '':
		ofile = fileStem+'.'+otype
	D_Obj = open(ifile).readlines()
	############################################################################################START OF CASE SCRIPTS
	if itype == 'soft':
		import georead
		OB = georead.ob_transform(D_Obj, identifier=idf, enum=True, targetCategory=targetCategory, targetClass=targetClass)
		print 'CLASS georead RETURNS ' + str(len(OB)-1) + ' features\n'
		print bcolors.BOLD + bcolors.HEADER +'ANALYZING DATASET' + bcolors.ENDC +'\n'
		print 'FEATURE NOMENCLATURE:'
		if idf != 'IDENTIFIER':
			print '\t'+idf+'\n'
		else:
			print '\t\'IDENTIFIER\' (DEFAULT)\n'
		if p_cutoff != 0:
			from geo2moses import pvalFilter
			OB = pvalFilter(OB,p_cutoff)
		M_ob = georead.ob2moses(OB)
		if not os.path.exists(outdir):
			os.makedirs(outdir)
			print 'OUTPUT DIRECTORY ./'+outdir+' NOT FOUND, CREATING DIRECTORY\n'
		if binarize == (True or 1):
			print 'BINARIZING DATASET'
			M_ob = georead.binarizeMoses(M_ob)
		##########################################################################################################
		if otype == 'moses':
			writeMoses(M_ob, outdir+'/'+ofile, '\t') 
		elif otype == 'svm':
			writeSVM(M_ob, outdir+'/'+ofile) 
		print str(M_ob.shape[0]-1) + ' Samples, ' + str(M_ob.shape[1]) + ' Features'
	elif itype == 'moses':
		print itype
	elif itype == 'svm':
		print itype
	elif itype in ['txt','dat','text']:
		#INFER INTERNAL DATA TYPE (SNP/MICROARRAY)
		print 'INPUT TYPE \''+itype+'\' (Generic flatfile) IS AMBIBUOUS, PLEASE SPECIFY DATA TYPE:'
		dtypes = {1:'microarray', 2:'snp'}
		for x in dtypes:
			print '\t'+str(x)+': '+dtypes[x]
		itype = dtypes[int(raw_input("\nCHOOSE DATA TYPE: "))]
		if itype == 'snp':
			print 'ENTERING SNP PROCESSING PORTOCOOL:'
			from snptransform import *
		print itype

if __name__ == "__main__":
	preprocess(sys.argv[1:])
