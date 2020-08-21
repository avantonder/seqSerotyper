#!/usr/bin/env python

import os, os.path, argparse, glob, sys, subprocess, csv, operator, numpy
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

#from Bio import SearchIO

#Check command line inputs exist

def check_file_exists(filepath, file_description):

	if not os.path.exists(filepath):
		print("The " + file_description + " " + filepath + " does not exist")
		sys.exit(1)

#Define blaster function

def blaster(blast_ref, blast_query, blastall_out):

	blastall_cmd = ["blastall", "-p", "blastn", "-d", blast_ref, "-i", blast_query, "-m 8", "-o", blastall_out]
	process = subprocess.call(blastall_cmd)

#Define blatter function - normal blast output

def blatter(blat_ref, blat_query, blat_output):

	blatter_cmd = ["blat", "-t=dna", blat_query, "-q=dna", blat_ref, blat_output, "out=blast"]
	process = subprocess.call(blatter_cmd)

#Define blatter8 function - blast8 output

def blatter8(blat8_ref, blat8_query, blat8_output):

	blatter8_cmd = ["blat", "-t=dna", blat8_query, "-q=dna", blat8_ref, blat8_output, "out=blast8"]
        process = subprocess.call(blatter8_cmd)
	
#Define mviewer function

def mviewer(mviewIn, mviewOut):

        mview_cmd = ["/nfs/users/nfs_a/avt/Software/mview/bin/mview", "-in", "blast", mviewIn, "-out", "pearson", "-top", "1"]
        process = subprocess.Popen(mview_cmd, stdout=subprocess.PIPE)

        mview_out = process.stdout.read()
        mview_outfile = open(mviewOut, "w")
        mview_outfile.write(mview_out)
        mview_outfile.close()

#Parse command line options

def parse_args(args):

	global _parser

	_parser = argparse.ArgumentParser(description = 'seqSerotyper: a tool for assigning serotypes to pneumococcal WGS assemblies.')
	_parser.add_argument('--ref_dir', '-r', help='Path to the directory containing the serotype reference files [REQUIRED]')
	_parser.add_argument('--input_dir', '-i', help='Path to the directory containing the query fasta file [REQUIRED]')
	_parser.add_argument('--fasta', '-f', help='Query fasta file [REQUIRED]')
	_parser.add_argument('--output_dir', '-o', help='Output directory [OPTIONAL]')

	opts = _parser.parse_args(args)

	return opts

#Main script

def main(opts):
		
	fasta_file = []
	ids = None	

	#If no output directory specified then create output directory in input directory

	if not opts.output_dir:
		if opts.fasta:
			opts.output_dir = os.path.join(os.path.dirname(opts.fasta), 'seqSerotyper')
			if not os.path.isdir(opts.output_dir): os.makedirs(opts.output_dir)
		else:
			opts.output_dir = os.path.join(opts.input_dir, 'seqSerotyper')
			if not os.path.isdir(opts.output_dir): os.makedirs(opts.output_dir)
	else:
		if not os.path.exists(opts.output_dir):os.makedirs(opts.output_dir)

	#Check command line files exist

	if opts.input_dir:
		check_file_exists(opts.input_dir, 'input directory')
		fasta_file = glob.glob(os.path.join(opts.input_dir, opts.fasta))

		if len(fasta_file) != 1:
			print "Incorrect or no fasta file specified. Please specify the correct fasta file"
			sys.exit(1)

		(seqDir,seqFileName) = os.path.split(fasta_file[0])
		(id,suffix) = seqFileName.split(".", 1)

	if opts.ref_dir:
		check_file_exists(opts.ref_dir, 'reference directory')

	#Blast fasta file against main reference file to assign initial serotype

	reference_file = glob.glob(os.path.join(opts.ref_dir, "Pneumo_serotypes.fasta"))
	
	blat_out = os.path.join(opts.output_dir, id + '.tab')
	blatter8(reference_file[0], fasta_file[0], blat_out)
	
    	tab = open(blat_out)
    	tabread = csv.reader(tab, delimiter = "\t")
    	tabsort = sorted(tabread, key = lambda x: float(x[11]), reverse=True)
    	tophit = tabsort[0]
    
    	top_serotype = tophit[0]

	hit_length = int(tophit[3])
    	perc_id = tophit[2]
 	tophit.extend([top_serotype])

    	sero_out = os.path.join(opts.output_dir, id + '_serotype.result')
    	sero_write = open(sero_out, 'w')
    	sero_write.write('Isolate\tSerotype\t%id\n')   	

    	#Serogroup 6

    	if tophit[12] in ('6A', '6B', '6C', '6D', '6Bii(6E)', '6F', '6G'):
    	
    		sero6CD_ref = glob.glob(os.path.join(opts.ref_dir, "6D_wciN_beta_nuc.fas"))
    		sero6CD_out = os.path.join(opts.output_dir, id + ".6all.tab")
    		blaster(sero6CD_ref[0], fasta_file[0], sero6CD_out)
		
		sero6CDtab = open(sero6CD_out)
		sero6CDtabread = csv.reader(sero6CDtab, delimiter = "\t")
		sero6CDtabsort = sorted(sero6CDtabread, key = lambda x: float(x[11]), reverse=True)
		sero6CDhit = sero6CDtabsort[0]
		sero6CDaln = int(sero6CDhit[3])
		
		if sero6CDaln < 500:
			sero6CDhit.extend(['Not6CD'])
		elif sero6CDaln > 500:
			sero6CDhit.extend(['6CD'])

		#Assign serotypes 6C and 6D

		if sero6CDhit[12] == '6CD':
			check6CD_ref = glob.glob(os.path.join(opts.ref_dir, "6A_6B_wzy_nuc.fas"))		
			check6CD_out = os.path.join(opts.output_dir, id + ".6CD.tab")
			blatter(check6CD_ref[0], fasta_file[0], check6CD_out)
			mview6CD_out = os.path.join(opts.output_dir, id + ".6CD.fasta")
			mviewer(check6CD_out, mview6CD_out)			
			
			for index, CDseq_record in enumerate(SeqIO.parse(mview6CD_out, "fasta")):
				print("fasta read")

			query6CD_seq = str(CDseq_record.seq)
			query6CD_dna = Seq(query6CD_seq, generic_dna)
			wzy_458 = query6CD_dna[457]

			if wzy_458 == "a":
				sero_write.write(id+'\t'+"6D"+'\t'+"wzy_checked"+'\n')
			elif wzy_458 == "-":
				sero_write.write(id+'\t'+"6C"+'\t'+"wzy_checked"+'\n')

		#Assign serotype 6F	

		elif sero6CDhit[12] == 'Not6CD':

			sero6ABF_ref = glob.glob(os.path.join(opts.ref_dir, "6A_wciN_alpha_nuc.fas"))
                        sero6ABF_out = os.path.join(opts.output_dir, id + ".6ABF.tab")
                        blatter(sero6ABF_ref[0], fasta_file[0], sero6ABF_out)
                        mview6ABF_out = os.path.join(opts.output_dir, id + ".6ABF.fasta")
                        mviewer(sero6ABF_out, mview6ABF_out)

                        for index, ABFseq_record in enumerate(SeqIO.parse(mview6ABF_out, "fasta")):
	                	print("fasta read")

                        query6ABF_seq = str(ABFseq_record.seq)
                        query6ABF_dna = Seq(query6ABF_seq, generic_dna)
                        query6ABF_aa = query6ABF_dna.translate()
                        wciN_150 = query6ABF_aa[149]

                        if wciN_150 == "T":
	                	sero_write.write(id+'\t'+"6F"+'\t'+"wciN_checked"+'\n')
				
			#Assign 6A, 6B and 6E
			
			else:

				not6CD_ref = glob.glob(os.path.join(opts.ref_dir, "6A_6B_wzy_nuc.fas"))
				not6CD_out = os.path.join(opts.output_dir, id + ".6EG.tab")
				blatter(not6CD_ref[0], fasta_file[0], not6CD_out)
				mview6EG_out = os.path.join(opts.output_dir, id + ".6EG.fasta")
				mviewer(not6CD_out, mview6EG_out)

				for index, EGseq_record in enumerate(SeqIO.parse(mview6EG_out, "fasta")):
					print("fasta read")
	
				query6EG_seq = str(EGseq_record.seq)
                        	query6EG_dna = Seq(query6EG_seq, generic_dna)
				wzy_658 = query6EG_dna[657]

				sero6AB_ref = glob.glob(os.path.join(opts.ref_dir, "6A_6B_wciP_nuc.fas"))
                                sero6AB_out = os.path.join(opts.output_dir, id + ".6AB.tab")
                                blatter(sero6AB_ref[0], fasta_file[0], sero6AB_out)
                                mview6AB_out = os.path.join(opts.output_dir, id + ".6AB.fasta")
                                mviewer(sero6AB_out, mview6AB_out)

                                for index, ABseq_record in enumerate(SeqIO.parse(mview6AB_out, "fasta")):
                                	print("fasta read")

                                query6AB_seq = str(ABseq_record.seq)
                                query6AB_dna = Seq(query6AB_seq, generic_dna)
                                query6AB_aa = query6AB_dna.translate()
                                wciP_195 = query6AB_aa[194]

				if wzy_658 == 'a' and wciP_195 == "S":
					sero_write.write(id+'\t'+"6E(6A)"+'\t'+"wciP_wzy_checked"+'\n')
				elif wzy_658 != 'a' and wciP_195 == "S":
					sero_write.write(id+'\t'+"6A"+'\t'+"wciP_wzy_checked"+'\n')
				elif wzy_658 == 'a' and wciP_195 == "N":				
					sero_write.write(id+'\t'+"6E(6B)"+'\t'+"wciP_wzy_checked"+'\n')
				elif wzy_658 != 'a' and wciP_195 == "N":
					sero_write.write(id+'\t'+"6B"+'\t'+"wciP_wzy_checked"+'\n')
				else:
                                        sero_write.write(id+'\t'+"6variant"+'\t'+"wciP_wzy_check_fail"+'\n')
					
	#Assign serotypes 7A and 7F

	elif tophit[12] in ('7A', '7F'):

		sero7AF_ref = glob.glob(os.path.join(opts.ref_dir, "7A_7F_wcwD_nuc.fas"))
		sero7AF_out = os.path.join(opts.output_dir, id + ".7AF.tab")
		blatter8(sero7AF_ref[0], fasta_file[0], sero7AF_out)
		sero7AFtab = open(sero7AF_out)
                sero7AFtabread = csv.reader(sero7AFtab, delimiter = "\t")
                sero7AFtabsort = sorted(sero7AFtabread, key = lambda x: float(x[11]), reverse=True)
                sero7AFhit = sero7AFtabsort[0]
                sero7AFaln = int(sero7AFhit[5])

                if sero7AFaln == 1:
			sero_write.write(id+'\t'+"7A"+'\t'+"wcwD_checked"+'\n')
		else:
			sero_write.write(id+'\t'+"7F"+'\t'+"wcwD_checked"+'\n')

	#Assign serotypes 7B, 7C and 40

	elif tophit[12] in ('7B', '7C', '40'):

		sero7BC_ref = glob.glob(os.path.join(opts.ref_dir, "7B_wcwK_nuc.fas"))
                sero7BC_out = os.path.join(opts.output_dir, id + ".7BC.tab")
                blatter(sero7BC_ref[0], fasta_file[0], sero7BC_out)
                mview7BC_out = os.path.join(opts.output_dir, id + ".7BC.fasta")
                mviewer(sero7BC_out, mview7BC_out)

                for index, s7BCseq_record in enumerate(SeqIO.parse(mview7BC_out, "fasta")):
         	       print("fasta read")

                query7BC_seq = str(s7BCseq_record.seq)
                query7BC_dna = Seq(query7BC_seq, generic_dna)
                query7BC_aa = query7BC_dna.translate()
                wcwK_16 = query7BC_aa[15]

                if wcwK_16 == "D":
 	        	sero_write.write(id+'\t'+"7B"+'\t'+"wcwK_checked"+'\n')
                elif wcwK_16 == "G":
                	sero_write.write(id+'\t'+"7C"+'\t'+"wcwK_checked"+'\n')
                elif wcwK_16 == "N":
                	sero_write.write(id+'\t'+"40"+'\t'+"wcwK_checked"+'\n')
		else:
			sero_write.write(id+'\t'+"7B/7C/40"+'\t'+"wcwK_check_fail"+'\n')

	#Assign serotypes 9A and 9V

	elif tophit[12] in ('9A', '9V'):

		sero9AV_ref = glob.glob(os.path.join(opts.ref_dir, "9A_9V_wcjE_nuc.fas"))
                sero9AV_out = os.path.join(opts.output_dir, id + ".9AV.tab")
		blatter8(sero9AV_ref[0], fasta_file[0], sero9AV_out)
		sero9AVtab = open(sero9AV_out)
                sero9AVtabread = csv.reader(sero9AVtab, delimiter = "\t")
                sero9AVtabsort = sorted(sero9AVtabread, key = lambda x: float(x[11]), reverse=True)
                sero9AVhit = sero9AVtabsort[0]
                sero9AValn = int(sero9AVhit[5])

		if sero9AValn == 1:
			sero_write.write(id+'\t'+"9A"+'\t'+"wcjE_checked"+'\n')
		else:
			sero_write.write(id+'\t'+"9V"+'\t'+"wcjE_checked"+'\n')

	#Assign serotypes 9L and 9N

	elif tophit[12] in ('9L', '9N'):

                sero9LN_ref = glob.glob(os.path.join(opts.ref_dir, "9L_wcjA_nuc.fas"))
                sero9LN_out = os.path.join(opts.output_dir, id + ".9LN.tab")
                blatter(sero9LN_ref[0], fasta_file[0], sero9LN_out)
                mview9LN_out = os.path.join(opts.output_dir, id + ".9LN.fasta")
                mviewer(sero9LN_out, mview9LN_out)

                for index, s9LNseq_record in enumerate(SeqIO.parse(mview9LN_out, "fasta")):
                       print("fasta read")

                query9LN_seq = str(s9LNseq_record.seq)
                query9LN_dna = Seq(query9LN_seq, generic_dna)
                query9LN_aa = query9LN_dna.translate()
                wcjA_139 = query9LN_aa[138]

                if wcjA_139 == "H":
                        sero_write.write(id+'\t'+"9L"+'\t'+"wcjA_checked"+'\n')
                elif wcjA_139 == "Y":
                        sero_write.write(id+'\t'+"9N"+'\t'+"wcjA_checked"+'\n')
		else:
			sero_write.write(id+'\t'+"9L/9N"+'\t'+"wcjA_check_fail"+'\n')

	#Serogroup 11

	elif tophit[12] in ('11A', '11D', '11F'):

		#Assign serotype 11D, 11F and 11X

		sero11X_ref = glob.glob(os.path.join(opts.ref_dir, "11X_wcrL_nuc.fas"))
		sero11X_out = os.path.join(opts.output_dir, id + ".11X.tab")
		blatter8(sero11X_ref[0], fasta_file[0], sero11X_out)
		sero11Xtab = open(sero11X_out)
		sero11Xtabread = csv.reader(sero11Xtab, delimiter = "\t")
		sero11Xtabsort = sorted(sero11Xtabread, key = lambda x: float(x[11]), reverse=True)		
		sero11Xhit = sero11Xtabsort[0]
		sero11Xaln = float(sero11Xhit[2])

		if sero11Xaln > 90:
			sero_write.write(id+'\t'+"11X"+'\t'+"wcrL_checked"+'\n')
		else:			
			sero11AD_ref = glob.glob(os.path.join(opts.ref_dir, "11A_wcrL_nuc.fas"))
			sero11AD_out = os.path.join(opts.output_dir, id + ".11AD.tab")
			blatter(sero11AD_ref[0], fasta_file[0], sero11AD_out) 	
			mview11AD_out = os.path.join(opts.output_dir, id + ".11AD.fasta")
			mviewer(sero11AD_out, mview11AD_out)

                	for index, ADseq_record in enumerate(SeqIO.parse(mview11AD_out, "fasta")):
	                	print("fasta read")

	                query11AD_seq = str(ADseq_record.seq)
		        query11AD_dna = Seq(query11AD_seq, generic_dna)
	                query11AD_aa = query11AD_dna.translate()
        	        wcrL_112 = query11AD_aa[111]

                	if wcrL_112 == "S":
        	        	sero_write.write(id+'\t'+"11D"+'\t'+"wcrL_checked"+'\n')
			elif wcrL_112 == "A":
				sero_write.write(id+'\t'+"11F"+'\t'+"wcrL_checked"+'\n')
			else:
				sero_write.write(id+'\t'+"11A"+'\t'+"wcrL_checked"+'\n')

	elif tophit[12] in ('11B', '11C'):

		sero11BC_ref = glob.glob(os.path.join(opts.ref_dir, "11C_gct_nuc.fas"))
		sero11BC_out = os.path.join(opts.output_dir, id + ".11BC.tab")
		blatter8(sero11BC_ref[0], fasta_file[0], sero11BC_out)
		sero11BCtab = open(sero11BC_out)
                sero11BCtabread = csv.reader(sero11BCtab, delimiter = "\t")
                sero11BCtabsort = sorted(sero11BCtabread, key = lambda x: float(x[11]), reverse=True)
                sero11BChit = sero11BCtabsort[0]
                sero11BCaln = float(sero11BChit[5])

		if sero11BCaln == 1:
                        sero_write.write(id+'\t'+"11B"+'\t'+"gct_checked"+'\n')
		else:
			sero_write.write(id+'\t'+"11C"+'\t'+"gct_checked"+'\n')

	#Assign serogroup 12 and serotypes 44 and 46

	elif tophit[12] in ('12A', '12B', '12F', '44', '46'):

		sero12_wcxF = glob.glob(os.path.join(opts.ref_dir, "12F_wcxF_nuc.fas"))
                sero12_wcxD = glob.glob(os.path.join(opts.ref_dir, "12F_wcxD_nuc.fas"))
		sero12_wcxF_out = os.path.join(opts.output_dir, id + ".12wcxF.tab")
		sero12_wcxD_out = os.path.join(opts.output_dir, id + ".12wcxD.tab")
                
		blatter(sero12_wcxF[0], fasta_file[0], sero12_wcxF_out)
                blatter(sero12_wcxD[0], fasta_file[0], sero12_wcxD_out)
		
		mview12_wcxF_out = os.path.join(opts.output_dir, id + ".12wcxF.fasta")
                mview12_wcxD_out = os.path.join(opts.output_dir, id + ".12wcxD.fasta")
		
		mviewer(sero12_wcxF_out, mview12_wcxF_out)
		mviewer(sero12_wcxD_out, mview12_wcxD_out)
		
                for index1, s12wcxF_record in enumerate(SeqIO.parse(mview12_wcxF_out, "fasta")):
                        print("fasta read")

		for index2, s12wcxD_record in enumerate(SeqIO.parse(mview12_wcxD_out, "fasta")):
			print("fasta read")

                query12wcxF_seq = str(s12wcxF_record.seq)
                query12wcxF_dna = Seq(query12wcxF_seq, generic_dna)
                query12wcxD_seq = str(s12wcxD_record.seq)
		query12wcxD_dna = Seq(query12wcxD_seq, generic_dna)

		wcxF_256 = query12wcxF_dna[255]
		wcxF_512 = query12wcxF_dna[511]
		wcxF_560 = query12wcxF_dna[559]
		wcxF_703 = query12wcxF_dna[702]
		wcxD_265 = query12wcxD_dna[264]
		wcxD_376 = query12wcxD_dna[375]
			
                if wcxF_256 == "g" and wcxF_512 == "c" and wcxF_560 == "t" and wcxF_703 == "a" and wcxD_265 == "c" and wcxD_376 == "c":
                        sero_write.write(id+'\t'+"12A"+'\t'+"wcxF_checked"+'\n')
		elif wcxF_256 == "g" and wcxF_512 == "c" and wcxF_560 == "t" and wcxF_703 == "c" and wcxD_265 == "c" and wcxD_376 == "c":
			sero_write.write(id+'\t'+"12B"+'\t'+"wcxF_checked"+'\n')
		elif wcxF_256 == "g" and wcxF_512 == "t" and wcxF_560 == "t" and wcxF_703 == "c" and wcxD_265 == "t" and wcxD_376 == "t":
			sero_write.write(id+'\t'+"12F"+'\t'+"wcxF_checked"+'\n')
		elif wcxF_256 == "g" and wcxF_512 == "c" and wcxF_560 == "c" and wcxF_703 == "c" and wcxD_265 == "t" and wcxD_376 == "t":
			sero_write.write(id+'\t'+"44"+'\t'+"wcxF_checked"+'\n')
		elif wcxF_256 == "a" and wcxF_512 == "c" and wcxF_560 == "t" and wcxF_703 == "c" and wcxD_265 == "c" and wcxD_376 == "c":
	        	sero_write.write(id+'\t'+"46"+'\t'+"wcxF_checked"+'\n')
		elif wcxF_256 == "g" and wcxF_512 == "c" and wcxF_560 == "t" and wcxF_703 == "c" and wcxD_265 == "t" and wcxD_376 == "c":
			sero_write.write(id+'\t'+"12F"+'\t'+"wcxF_checked"+'\n')
		else:
			sero_write.write(id+'\t'+"12/44/46"+'\t'+"wcxF_check_fail"+'\n')

	#Assign serotypes 15A, 15B/C, 15F

	elif tophit[12] in ('15A', '15B/15C', '15F'):

		#Assign serotype 15F
				
		sero15F_ref = glob.glob(os.path.join(opts.ref_dir, "15F_wcjE_nuc.fas"))
		sero15F_out = os.path.join(opts.output_dir, id + ".15F.tab")
		blaster(sero15F_ref[0], fasta_file[0], sero15F_out)
		
		sero15Ftab = open(sero15F_out)
                sero15Ftabread = csv.reader(sero15Ftab, delimiter = "\t")
                sero15Ftabsort = sorted(sero15Ftabread, key = lambda x: float(x[11]), reverse=True)
                sero15Fhit = sero15Ftabsort[0]
                sero15Faln = int(sero15Fhit[3])

                if sero15Faln < 1000:
                        sero15Fhit.extend(['Not15F'])
                elif sero15Faln > 1000:
			sero15Fhit.extend(['15F'])
			sero_write.write(id+'\t'+"15F"+'\t'+"wcjE_checked"+'\n')
                        
                #Assign serotypes 15A and 15B/C

                if sero15Fhit[12] == 'Not15F':
			
			sero15A_ref = glob.glob(os.path.join(opts.ref_dir, "15A_wchL_nuc.fas"))
	                sero15A_out = os.path.join(opts.output_dir, id + ".15A.tab")
        	        blatter8(sero15A_ref[0], fasta_file[0], sero15A_out)
                	sero15Atab = open(sero15A_out)
                	sero15Atabread = csv.reader(sero15Atab, delimiter = "\t")
                	sero15Atabsort = sorted(sero15Atabread, key = lambda x: float(x[11]), reverse=True)
                	sero15Ahit = sero15Atabsort[0]
                	sero15Aaln = float(sero15Ahit[2])

                	if sero15Aaln > 99.8:
                        	sero_write.write(id+'\t'+"15A"+'\t'+"wchL_checked"+'\n')
                	elif sero15Aaln < 99.8:
                        	sero_write.write(id+'\t'+"15B/15C"+'\t'+"wchL_checked"+'\n')			
		else:
			print("serogroup 15 analysed")

	#Assign serotypes 18B and 18C

	elif tophit[12] in ('18B', '18C'):

		sero18BC_ref = glob.glob(os.path.join(opts.ref_dir, "18B_18C_wciX_nuc.fas"))
		sero18BC_out = os.path.join(opts.output_dir, id + ".18BC.tab")
                blatter(sero18BC_ref[0], fasta_file[0], sero18BC_out)
  		mview18BC_out = os.path.join(opts.output_dir, id + ".18BC.fasta")
                mviewer(sero18BC_out, mview18BC_out)

                for index, s18BCseq_record in enumerate(SeqIO.parse(mview18BC_out, "fasta")):
                        print("fasta read")

                query18BC_seq = str(s18BCseq_record.seq)
                query18BC_dna = Seq(query18BC_seq, generic_dna)
                wciX_169 = query18BC_dna[168]             

                if wciX_169 == "t":
	                sero_write.write(id+'\t'+"18B"+'\t'+"wciX_checked"+'\n')
                elif wciX_169 == "g":
                        sero_write.write(id+'\t'+"18C"+'\t'+"wciX_checked"+'\n')
		else:
			sero_write.write(id+'\t'+"18B/18C"+'\t'+"wciX_check_fail"+'\n')

	#Assign serotypes 19A and 19A(19F)

	elif tophit[12] == '19A':

		sero19AF_ref = glob.glob(os.path.join(opts.ref_dir, "19AF_wzy_nuc.fas"))
                sero19AF_out = os.path.join(opts.output_dir, id + ".19AF.tab")

                blaster(sero19AF_ref[0], fasta_file[0], sero19AF_out)
                sero19AFtab = open(sero19AF_out)
                sero19AFtabread = csv.reader(sero19AFtab, delimiter = "\t")
                sero19AFtabsort = sorted(sero19AFtabread, key = lambda x: float(x[11]), reverse=True)
                sero19AFhit = sero19AFtabsort[0]
                sero19AFaln = float(sero19AFhit[2])

                if sero19AFaln < 99:
                        sero_write.write(id+'\t'+"19A"+'\t'+"wzy_checked"+'\n')
                elif sero19AFaln > 99:
                        sero_write.write(id+'\t'+"19A(19F)"+'\t'+"wzy_checked"+'\n')

	#Assign serotypes 19B and 19C

	elif tophit[12] in ('19B', '19C'):

		sero19BC_ref = glob.glob(os.path.join(opts.ref_dir, "19C_wchU_nuc.fas"))
                sero19BC_out = os.path.join(opts.output_dir, id + ".19BC.tab")
		
		blaster(sero19BC_ref[0], fasta_file[0], sero19BC_out)
                sero19BCtab = open(sero19BC_out)
                sero19BCtabread = csv.reader(sero19BCtab, delimiter = "\t")
                sero19BCtabsort = sorted(sero19BCtabread, key = lambda x: float(x[11]), reverse=True)
                sero19BChit = sero19BCtabsort[0]
                sero19BCaln = int(sero19BChit[3])

                if sero19BCaln < 500:
                        sero_write.write(id+'\t'+"19B"+'\t'+"wchU_checked"+'\n')
                elif sero19BCaln > 500:
                        sero_write.write(id+'\t'+"19C"+'\t'+"wchU_checked"+'\n')

	#Assign serotypes 22A and 22F

	elif tophit[12] in ('22A', '22F'):

		sero22AF_ref = glob.glob(os.path.join(opts.ref_dir, "22F_wcwA_nuc.fas"))
		sero22AF_out = os.path.join(opts.output_dir, id + ".22AF.tab")
		blaster(sero22AF_ref[0], fasta_file[0], sero22AF_out)
		sero22AFtab = open(sero22AF_out)
                sero22AFtabread = csv.reader(sero22AFtab, delimiter = "\t")
                sero22AFtabsort = sorted(sero22AFtabread, key = lambda x: float(x[11]), reverse=True)
                sero22AFhit = sero22AFtabsort[0]
                sero22AFaln = int(sero22AFhit[3])

                if sero22AFaln < 500:
                	sero_write.write(id+'\t'+"22A"+'\t'+"wcwA_checked"+'\n')
                elif sero22AFaln > 500:
                        sero_write.write(id+'\t'+"22F"+'\t'+"wcwA_checked"+'\n')

	#Assign serotypes 25A and 25F

	elif tophit[12] in ('25A', '25F'):

		sero25AF_ref = glob.glob(os.path.join(opts.ref_dir, "25A_wcyC_nuc.fas"))
                sero25AF_out = os.path.join(opts.output_dir, id + ".25AF.tab")
		blatter(sero25AF_ref[0], fasta_file[0], sero25AF_out)
		mview25AF_out = os.path.join(opts.output_dir, id + ".25AF.fasta")
                mviewer(sero25AF_out, mview25AF_out)
 
                for index, s25AFseq_record in enumerate(SeqIO.parse(mview25AF_out, "fasta")):
 	               print("fasta read")

                query25AF_seq = str(s25AFseq_record.seq)
                query25AF_dna = Seq(query25AF_seq, generic_dna)
                wcyC_430 = query25AF_dna[429]
				
                if wcyC_430 == "g":
        	        sero_write.write(id+'\t'+"25A"+'\t'+"wcyC_checked"+'\n')
		else:
			sero_write.write(id+'\t'+"25F"+'\t'+"wcyC_checked"+'\n')

	#Assign serotypes 28A and 28F

	elif tophit[12] in ('28A', '28F'):

		sero28AF_ref = glob.glob(os.path.join(opts.ref_dir, "28A_wciU_nuc.fas"))
                sero28AF_out = os.path.join(opts.output_dir, id + ".28AF.tab")
                blatter(sero28AF_ref[0], fasta_file[0], sero28AF_out)
                mview28AF_out = os.path.join(opts.output_dir, id + ".28AF.fasta")
                mviewer(sero28AF_out, mview28AF_out)

                for index, s28seq_record in enumerate(SeqIO.parse(mview28AF_out, "fasta")):
 	               print("fasta read")

                query28AF_seq = str(s28seq_record.seq)
                query28AF_dna = Seq(query28AF_seq, generic_dna)
                query28AF_aa = query28AF_dna.translate()
                wciU_40 = query28AF_aa[39]
		wciU_72 = query28AF_aa[71]

                if wciU_40 == "C" and wciU_72 == "I":
        	        sero_write.write(id+'\t'+"28A"+'\t'+"wciU_checked"+'\n')
                elif wciU_40 == "Y" and wciU_72 == "M":
                	sero_write.write(id+'\t'+"28F"+'\t'+"wciU_checked"+'\n')
		else:
			sero_write.write(id+'\t'+"28var"+'\t'+"wciU_checked"+'\n')
                        
	#Assign serotypes 33A, 33F, 33X and 37
	
	#Assign serotype 37
		
	elif tophit[12] in ('33A', '33F'):

		sero37_ref = glob.glob(os.path.join(opts.ref_dir, "37_tts_nuc.fas"))
		sero37_out = os.path.join(opts.output_dir, id + ".37.tab")
		blaster(sero37_ref[0], fasta_file[0], sero37_out)
		sero37tab = open(sero37_out)
                sero37tabread = csv.reader(sero37tab, delimiter = "\t")
                sero37tabsort = sorted(sero37tabread, key = lambda x: float(x[11]), reverse=True)
                sero37hit = sero37tabsort[0]
                sero37aln = int(sero37hit[3])

                if sero37aln > 500:
                        sero_write.write(id+'\t'+"37"+'\t'+"tts present"+'\n')
                elif sero37aln < 500:
		
			#Assign serotype 33X

			sero33X_ref = glob.glob(os.path.join(opts.ref_dir, "33A_33F_wcjE_nuc.fas"))
			sero33X_out = os.path.join(opts.output_dir, id + ".33X.tab")
                	blaster(sero33X_ref[0], fasta_file[0], sero33X_out)
                	sero33Xtab = open(sero33X_out)
                	sero33Xtabread = csv.reader(sero33Xtab, delimiter = "\t")
                	sero33Xtabsort = sorted(sero33Xtabread, key = lambda x: float(x[11]), reverse=True)
                	sero33Xhit = sero33Xtabsort[0]
                	sero33Xaln = int(sero33Xhit[3])

                	if sero33Xaln < 300:
                        	sero_write.write(id+'\t'+"33X"+'\t'+"no_wcjE_present"+'\n')
                	elif sero33Xaln > 300:
		
				#Assign serotypes 33A and 33F

				sero33AF_ref = glob.glob(os.path.join(opts.ref_dir, "33A_33F_wcjE_nuc.fas"))
        	        	sero33AF_out = os.path.join(opts.output_dir, id + ".33AF.tab")			
				blatter8(sero33AF_ref[0], fasta_file[0], sero33AF_out)
				sero33AFtab = open(sero33AF_out)
                		sero33AFtabread = csv.reader(sero33AFtab, delimiter = "\t")
                		sero33AFtabsort = sorted(sero33AFtabread, key = lambda x: float(x[11]), reverse=True)
                		sero33AFhit = sero33AFtabsort[0]
                		sero33AFaln = int(sero33AFhit[3])

                		if sero33AFaln > 1075:
                        		sero_write.write(id+'\t'+"33A"+'\t'+"wcjE_checked"+'\n')
                		else:
                        		sero_write.write(id+'\t'+"33F"+'\t'+"wcjE_checked"+'\n')

	#Assign serotypes 35A, 35C and 42

	elif tophit[12] in ('35A', '35C', '42'):

		#Assign serotype 35A

		sero35AC_ref = glob.glob(os.path.join(opts.ref_dir, "35A_mnp1_nuc.fas"))
                sero35AC_out = os.path.join(opts.output_dir, id + ".35AC.tab")
                blatter(sero35AC_ref[0], fasta_file[0], sero35AC_out)
                mview35AC_out = os.path.join(opts.output_dir, id + ".35AC.fasta")
                mviewer(sero35AC_out, mview35AC_out)

                for index, s35seq_record in enumerate(SeqIO.parse(mview35AC_out, "fasta")):
                       print("fasta read")

                query35AC_seq = str(s35seq_record.seq)
                query35AC_dna = Seq(query35AC_seq, generic_dna)
                query35AC_aa = query35AC_dna.translate()
                mnp_67 = query35AC_aa[66]
                mnp_97 = query35AC_aa[96]

                if mnp_67 == "R" and mnp_97 == "V":
                        sero_write.write(id+'\t'+"35A"+'\t'+"mnp_checked"+'\n')
                else:

			#Assign serotypes 35C and 42

			sero42_ref = glob.glob(os.path.join(opts.ref_dir, "35C_wzh_nuc.fas"))
                	sero42_out = os.path.join(opts.output_dir, id + ".42.tab")
                	blatter(sero42_ref[0], fasta_file[0], sero42_out)
                	mview42_out = os.path.join(opts.output_dir, id + ".42.fasta")
                	mviewer(sero42_out, mview42_out)

                	for index, s42seq_record in enumerate(SeqIO.parse(mview42_out, "fasta")):
                        	print("fasta read")

                	query42_seq = str(s42seq_record.seq)
                	query42_dna = Seq(query42_seq, generic_dna)
                	query42_aa = query42_dna.translate()
                	wzh_156 = query42_aa[155]
                	
			if wzh_156 == "R":
                        	sero_write.write(id+'\t'+"35C"+'\t'+"wzh_checked"+'\n')
                	elif wzh_156 == "M":
                        	sero_write.write(id+'\t'+"42"+'\t'+"wzh_checked"+'\n')
			else:
				sero_write.write(id+'\t'+"35A/35C/42"+'\t'+"wzh_check_fail"+'\n')

	#Assign serotypes 35B and 35D

	elif tophit[12] == '35B':

		sero35BD_ref = glob.glob(os.path.join(opts.ref_dir, "35B_wciG_nuc.fas"))
                sero35BD_out = os.path.join(opts.output_dir, id + ".35BD.tab")
		blatter8(sero35BD_ref[0], fasta_file[0], sero35BD_out)
                sero35BDtab = open(sero35BD_out)
                sero35BDtabread = csv.reader(sero35BDtab, delimiter = "\t")
                sero35BDtabsort = sorted(sero35BDtabread, key = lambda x: float(x[11]), reverse=True)
                sero35BDhit = sero35BDtabsort[0]
                sero35BDaln = int(sero35BDhit[5])

                if sero35BDaln > 0:
	                sero_write.write(id+'\t'+"35D"+'\t'+"wciG_checked"+'\n')
                else:
        	        sero_write.write(id+'\t'+"35B"+'\t'+"wciG_checked"+'\n')

	#Assign other serotypes that don't require further classification	 
	    	
   	else:
    		sero_write .write(id+'\t'+tophit[12]+'\t'+str(perc_id)+'\n')
    
    	sero_write.close()
    
if __name__ == "__main__":
  opts= parse_args(sys.argv[1:])
  main(opts)

		
