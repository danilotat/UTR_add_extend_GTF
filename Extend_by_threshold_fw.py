#!/usr/bin/env python3

from collections import OrderedDict,namedtuple
from gtf_advanced_parser import *
import argparse

class GTF_element(object):
	def __init__(self,chromosome,source,feature_type,start,end,score,strand,phase,attributes):
		self.chromosome=chromosome
		self.source=source
		self.feature_type=feature_type
		self.start=int(start) 
		self.end=int(end) 
		self.score=score
		self.strand=strand
		self.phase=phase
		self.attributes=attributes 

	def gtf_line(self):
		'''Passing a GTF_element, return a line as string'''
		return f'{self.chromosome}\t{self.source}\t{self.feature_type}\t{self.start}\t{self.end}\t{self.score}\t{self.strand}\t{self.phase}\t{self.attributes}'

def extend_return(dict_type,threshold,transcript_attributes):
	try:
		for item in dict_type[transcript_attributes['transcript_id']]:
			item_as_gtf=GTF_element(*item.split('\t'))
			if item_as_gtf.feature_type == 'stop_codon' or item_as_gtf.feature_type == 'CDS':
				print(item_as_gtf.gtf_line())
			else:
				if item_as_gtf.end == gene_id_coords[prev_gene_id][2]:
					item_as_gtf.end+=threshold
					print(item_as_gtf.gtf_line())
				else:
					print(item_as_gtf.gtf_line())
	except KeyError:
		pass
def iterate_and_extend(gene_id,threshold):
	gene_record = GTF_element(*gene_id_records[gene_id].split('\t'))
	gene_record.end+=threshold
	gene_gtf_line=gene_record.gtf_line()
	print(gene_gtf_line)
	for transcript in transcript_id_records[gene_id]:
		transcript_as_gtf=GTF_element(*transcript.split('\t'))
		transcript_attributes=parse_attributes(transcript_as_gtf.attributes)
		if transcript_as_gtf.end == gene_id_coords[gene_id][2]:
			transcript_as_gtf.end+=threshold
			transcript_gtf_line=transcript_as_gtf.gtf_line()
			print(transcript_gtf_line)
			# iterate through dict to return features
			for k in list_of_dicts:
				extend_return(k,threshold,transcript_attributes)
		else:
			print(transcript)
			#internal transcript
			for dict_type in list_of_dicts:
				try:
					for record in dict_type[transcript_attributes['transcript_id']]:
						print(record)
				except KeyError:
					continue 
	try:
		for misc in misc_features[gene_id]:
			print(misc)
	except KeyError:
		pass

def iterate_not_extend(gene_id):
	gene_record = GTF_element(*gene_id_records[gene_id].split('\t'))
	gene_gtf_line=gene_record.gtf_line()
	print(gene_gtf_line)
	for transcript in transcript_id_records[gene_id]:
		transcript_as_gtf=GTF_element(*transcript.split('\t'))
		transcript_attributes=parse_attributes(transcript_as_gtf.attributes)
		print(transcript)
		for k in list_of_dicts:
			try:
				for record in k[transcript_attributes['transcript_id']]:
					print(record)
			except KeyError:
				continue
	try:
		for misc in misc_features[gene_id]:
			print(misc)
	except KeyError:
		pass
		
# ---- PARSE ARGUMENTS ----

parser = argparse.ArgumentParser(description="Extend 3' terminus of a GTF forward")
parser.add_argument("--input", help="Input file. Must be a fw only version of a valid GTF")
parser.add_argument("--threshold", type=int, default=1000, help="Threshold to use for 3'-terminal extension. Default 1000")
args=parser.parse_args()

# check if no argument was passed 
if args.input == None:
	print('No input provided. Exiting..')
	exit()



fill_dict(gtf_file=args.input)
threshold=args.threshold



#first make a list of keys
list_keys=[]
for k in gene_id_coords.keys():
	list_keys.append(k)

#iterate through records
list_of_dicts=[exons_records,cdss_records,five_prime_utrs_records,three_prime_utrs_records,start_stop_codons]
for idx,key in enumerate(list_keys):
	if idx == 0:
		continue
	elif 0 < idx < len(list_keys) - 1:
		cur_gene_id=list_keys[idx]
		prev_gene_id=list_keys[idx-1]
		#check if prev_gene is coding
		prev_gene_as_gtf=GTF_element(*gene_id_records[prev_gene_id].split('\t'))
		prev_gene_attributes=parse_attributes(prev_gene_as_gtf.attributes)
		if prev_gene_attributes['gene_biotype'] == '"protein_coding"':
			#check chromosome if same
			if gene_id_coords[cur_gene_id][0] == gene_id_coords[prev_gene_id][0]:
				residue=int(gene_id_coords[cur_gene_id][1]) - int(gene_id_coords[prev_gene_id][2])
				#calc residue if higher than threshold
				if residue > threshold:
					iterate_and_extend(prev_gene_id,threshold)
				elif 0 < residue < threshold:
					# here extension could be done by a value lower than threshold
					corrected_threshold=int(residue-1)
					iterate_and_extend(prev_gene_id,corrected_threshold)
				elif int(gene_id_coords[cur_gene_id][1]) < int(gene_id_coords[prev_gene_id][2]):
					iterate_not_extend(prev_gene_id)
			elif gene_id_coords[cur_gene_id][0] != gene_id_coords[prev_gene_id][0]:
				#chromosome is different. Extend if possible
				iterate_and_extend(prev_gene_id,threshold)
		else:
			#not coding gene
			iterate_not_extend(prev_gene_id)
	elif idx == len(list_keys) - 1:
		cur_gene_id=list_keys[idx]
		prev_gene_id=list_keys[idx-1]
		#check if prev_gene is coding
		prev_gene_as_gtf=GTF_element(*gene_id_records[prev_gene_id].split('\t'))
		prev_gene_attributes=parse_attributes(prev_gene_as_gtf.attributes)
		if prev_gene_attributes['gene_biotype'] == '"protein_coding"':
			#check chromosome if same
			if gene_id_coords[cur_gene_id][0] == gene_id_coords[prev_gene_id][0]:
				residue=int(gene_id_coords[cur_gene_id][1]) - int(gene_id_coords[prev_gene_id][2])
				#calc residue if higher than threshold
				if residue > threshold:
					iterate_and_extend(prev_gene_id,threshold)
				elif 0 < residue < threshold:
					# here extension could be done by a value lower than threshold
					corrected_threshold=int(residue-1)
					iterate_and_extend(prev_gene_id,corrected_threshold)
				elif int(gene_id_coords[cur_gene_id][1]) < int(gene_id_coords[prev_gene_id][2]):
					iterate_not_extend(prev_gene_id)
			elif gene_id_coords[cur_gene_id][0] != gene_id_coords[prev_gene_id][0]:
				#chromosome is different. Extend if possible
				iterate_and_extend(prev_gene_id,threshold)
		else:
			#not coding gene
			iterate_not_extend(prev_gene_id)
		#last entry, check only if coding
		cur_gene_as_gtf=GTF_element(*gene_id_records[cur_gene_id].split('\t'))
		cur_gene_attributes=parse_attributes(cur_gene_as_gtf.attributes)
		if cur_gene_attributes['gene_biotype'] == '"protein_coding"':
			iterate_and_extend(cur_gene_id,threshold)
		else:
			iterate_not_extend(cur_gene_id)

