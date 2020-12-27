#!/usr/bin/python3
from collections import namedtuple, OrderedDict
import argparse

''' --- FUNCTIONS BLOCK --- '''
def gene_dict(gene_id, record):
	''' Pass gene id, add to dict'''
	gene_id_records[gene_id]=record
def gene_coords(gene_id, chromosome, start, end):
	''' Passing gene_id, save mapping positions '''
	gene_id_coords.setdefault(gene_id, []).extend([chromosome, int(start), int(end)])
def transcript_dict(gene_id, record):
	''' Pass gene id, create dict of transcripts'''
	transcript_id_records.setdefault(gene_id, []).append(record)
def cds_dict(transcript_id, record):
	''' Pass transcript id, create list of CDSs'''
	cdss_records.setdefault(transcript_id, []).append(record)
def exons_dict(transcript_id, record):
	''' Pass transcript id, create list of exons'''
	exons_records.setdefault(transcript_id, []).append(record)
def five_utrs_dict(transcript_id, record):
	''' Pass transcript id, create list of 5' utr '''
	five_prime_utrs_records.setdefault(transcript_id, []).append(record)
def three_utrs_dict(transcript_id, record):
	''' Pass transcript id, create list of 3' utr '''
	three_prime_utrs_records.setdefault(transcript_id, []).append(record)	
def codons_dict(transcript_id, record):
	''' Populate dict of start & stop codons'''
	start_stop_codons.setdefault(transcript_id, []).append(record)
def misc_records(gene_id, record):
	''' Populate misc dict '''
	misc_features.setdefault(gene_id, []).append(record)
def parse_attributes(kvps):
    ''' Parse the 9th column of GTF into a dictionary of attributes.'''
    f_dict = {}
    kvp = kvps.split(';')[:-1]
    for i in kvp:
        ri = i.replace(' "', '="')
        try:
            k, v = ri.split('=')
            rk = k.replace(' ', '')
            f_dict[rk] = v
        except ValueError:
            print('Malformed file.\nThis line isn\'t good. %s' % kvp)
            exit()
    return f_dict
def fill_dict(gtf_file):
	''' Fill dictionaries using provided features '''
	with open(gtf_file, 'r') as gtf_file:
		for line in gtf_file:
			if line.startswith('#'):
				hashtag_lines.append(line)
			else:
				record = GTF_record(*line.rstrip().split('\t'))
				record_attributes=parse_attributes(record.attributes)
				if record.feature_type == 'gene':
					gene_dict(record_attributes['gene_id'], line.rstrip())
					#populate also coords dict
					gene_coords(record_attributes['gene_id'], record.chromosome, record.start, record.end)
				elif record.feature_type == 'transcript':
					transcript_dict(record_attributes['gene_id'], line.rstrip())
				elif record.feature_type == 'CDS':
					cds_dict(record_attributes['transcript_id'], line.rstrip())
				elif record.feature_type == 'exon':
					exons_dict(record_attributes['transcript_id'], line.rstrip())
				elif record.feature_type == 'three_prime_utr':
					three_utrs_dict(record_attributes['transcript_id'], line.rstrip())
				elif record.feature_type == 'five_prime_utr':
					five_utrs_dict(record_attributes['transcript_id'], line.rstrip())
				elif record.feature_type == 'stop_codon' or record.feature_type == 'start_codon':
					codons_dict(record_attributes['transcript_id'], line.rstrip())
				else:
					misc_records(record_attributes['gene_id'], line.rstrip())
def create_three_prime_utr_fw(last_cds_end,transcript_record,transcript_attributes):
	''' Return a 3'UTR string to be fitted into utrs_dictionary '''
	attr_as_string='; '.join(['{} {}'.format(k,v) for k,v in transcript_attributes.items()])
	return f'{transcript_record.chromosome}\t{transcript_record.source}\t{"three_prime_utr"}\t{last_cds_end}\t{transcript_record.end}\t{transcript_record.score}\t{transcript_record.strand}\t{transcript_record.phase}\t{attr_as_string}'
def create_three_prime_utr_rev(first_cds_start, transcript_record, transcript_attributes):
	''' Return a 3'UTR string to be fitted into utrs_dict'''
	attr_as_string='; '.join(['{} {}'.format(k,v) for k,v in transcript_attributes.items()])
	return f'{transcript_record.chromosome}	{transcript_record.source}	{"three_prime_utr"}	'\
		f'{transcript_record.start}	{first_cds_start}	{transcript_record.score}	{transcript_record.strand}	'\
		f'{transcript_record.phase}	{attr_as_string}'
def fill_missing_utrs(gene_id_records):
	global added_three_prime_utr_fw, added_three_prime_utr_rev
	for gene_id in gene_id_records.keys():
		#check if gene is coding 
		gene_record = GTF_record(*gene_id_records[gene_id].rstrip().split('\t'))
		gene_attributes = parse_attributes(gene_record.attributes)
		if gene_attributes['gene_biotype'] == '"protein_coding"':
			#get transcripts: iterate through them
			transcripts = transcript_id_records[gene_id]
			for transcript in transcripts:
				#get the transcript id
				transcript_record = GTF_record(*transcript.rstrip().split('\t'))
				transcript_attributes=parse_attributes(transcript_record.attributes)
				transcript_id=transcript_attributes['transcript_id']
				# check existance of 3' UTR record
				# first check strand
				if gene_record.strand == '+':
					try:
						utr_three=three_prime_utrs_records[transcript_id]
					except KeyError:
						# no 3'UTR for this transcript: add one 
						#get last CDS entry for this transcript
						last_cds= GTF_record(*cdss_records[transcript_id][-1].split('\t')) 
						#try to make a novel 3'utr record
						novel_utr=create_three_prime_utr_fw(last_cds.end,transcript_record,transcript_attributes)
						#add to utrs_dict
						three_prime_utrs_records[transcript_id]=[]
						three_prime_utrs_records[transcript_id].append(novel_utr)
						added_three_prime_utr_fw+=1
				elif gene_record.strand == '-':
					try:
						utr_three=three_prime_utrs_records[transcript_id]
					except KeyError:
						#no 3'UTR for this transcript
						#get first CDS for this transcript 
						first_cds=GTF_record(*cdss_records[transcript_id][0].split('\t'))
						novel_utr=create_three_prime_utr_rev(first_cds.start,transcript_record,transcript_attributes)
						added_three_prime_utr_rev+=1
						# add to utrs_dict
						three_prime_utrs_records[transcript_id]= []
						three_prime_utrs_records[transcript_id].append(novel_utr)
				else:
					raise Exception("ERROR! This gene contain wrong information about strand: %s" %gene_id_records[gene_id])
def return_features_from_dict(dict_of_features, transcript_id,output_file):
	'''Passing transcript_id, return features: continue if none is present'''
	try:
		for feature in dict_of_features[transcript_id]:
			output_file.write(feature+'\n')
	except KeyError:
		pass
def gtf_rebuild(output_gtf):
	''' Rebuild the GTF file '''
	with open(output_gtf, 'w') as oput_gtf:
	#first write hashtag lines
		for hashtag_line in hashtag_lines:
			oput_gtf.write(hashtag_line)
		for gene_id in gene_id_records.keys():
			oput_gtf.write(gene_id_records[gene_id]+'\n')
			#now transcripts
			for transcript in transcript_id_records[gene_id]:
				oput_gtf.write(transcript+'\n')
				transcript_record=GTF_record(*transcript.rstrip().split('\t'))
				transcript_attributes=parse_attributes(transcript_record.attributes)
				transcript_id=transcript_attributes['transcript_id']
				# now print exons
				for exon in exons_records[transcript_id]:
					oput_gtf.write(exon+'\n')
				#now print CDS: they're present only if the gene is coding
				if transcript_attributes['transcript_biotype'] == '"protein_coding"':
					list_of_dicts=[cdss_records,three_prime_utrs_records,five_prime_utrs_records,start_stop_codons]
					for k in list_of_dicts:
						return_features_from_dict(k,transcript_id,oput_gtf)
				else:
					pass			
			try:
				for misc in misc_features[gene_id]:
					oput_gtf.write(misc+'\n')
			except KeyError:
				pass	
''' ---- DATA STRUCTURES ------ '''
gtf_fields=['chromosome','source','feature_type','start','end','score','strand','phase','attributes']
GTF_record=namedtuple('GTFrecord',gtf_fields)
hashtag_lines=[]
gene_id_records=OrderedDict()
gene_id_coords=OrderedDict()
transcript_id_records=OrderedDict()
cdss_records=OrderedDict()
exons_records=OrderedDict()
five_prime_utrs_records=OrderedDict()
three_prime_utrs_records=OrderedDict()
start_stop_codons=OrderedDict()
misc_features=OrderedDict()

'''---- LOG VARIABLES --- '''
added_three_prime_utr_fw=0
added_three_prime_utr_rev=0

''' --- CODE --- '''
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Giving an Ensembl GTF file, add explicit 3'UTR when missing to protein coding genes.")
	parser.add_argument("--input", help="Input file. Ensembl GTF only, other versions could lead to errors.")
	parser.add_argument("--output", help="Output GTF file")
	args = parser.parse_args()

	if (args.input == None or args.output == None):
		parser.print_help()
		exit()
	#populate those dictionaries. 
	fill_dict(gtf_file=args.input)
	#add missing three prime UTR
	fill_missing_utrs(gene_id_records)
	#rebuild the output GTF
	gtf_rebuild(output_gtf=args.output)
	#cheers
	print("We finished here. Here's some info:")
	print("Processed genes: %s" %len(gene_id_records.keys()))
	print("Processed transcripts: %s" %len(cdss_records.keys()))
	print("Added 3'UTR onto forward strand: %s" %added_three_prime_utr_fw)
	print("Added 3'UTR onto reverse strand: %s" %added_three_prime_utr_rev)
	print('Report any error to @danilotat [ github.com/danilotat ]'+'\n')
	print('Enjoy!')

