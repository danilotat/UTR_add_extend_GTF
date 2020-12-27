from collections import OrderedDict,namedtuple
from gtf_advanced_parser import *

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

def extend_return_rev(dict_type, threshold, transcript_attributes):
	try:
		for item in dict_type[transcript_attributes['transcript_id']]:
			item_as_gtf=GTF_element(*item.split('\t'))
			if item_as_gtf.feature_type == 'stop_codon' or item_as_gtf.feature_type == 'CDS':
				print(item_as_gtf.gtf_line())
			else:
				if item_as_gtf.start == gene_id_coords[cur_gene_id][1]:
					item_as_gtf.start-=threshold
					print(item_as_gtf.gtf_line())
				else:
					print(item_as_gtf.gtf_line())
	except KeyError:
		pass
def iterate_and_extend_rev(gene_id, threshold):
	gene_record = GTF_element(*gene_id_records[gene_id].split('\t'))
	gene_record.start-=threshold
	gene_gtf_line=gene_record.gtf_line()
	print(gene_gtf_line)
	for transcript in transcript_id_records[gene_id]:
		transcript_as_gtf=GTF_element(*transcript.split('\t'))
		transcript_attributes=parse_attributes(transcript_as_gtf.attributes)
		if transcript_as_gtf.start == gene_id_coords[gene_id][1]:
			transcript_as_gtf.start-=threshold
			transcript_gtf_line=transcript_as_gtf.gtf_line()
			print(transcript_gtf_line)
			# iterate through dict to return features
			for k in list_of_dicts:
				extend_return_rev(k,threshold,transcript_attributes)
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

def iterate_not_extend_rev(gene_id):
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
#first make a list of keys

parser = argparse.ArgumentParser(description="Extend 3' terminus of a GTF reverse")
parser.add_argument("--input", help="Input file. Must be a rev only version of a valid GTF")
parser.add_argument("--threshold", type=int, default=1000, help="Threshold to use for 3'-terminal extension. Default 1000")
args=parser.parse_args()

# check if no argument was passed 
if args.input == None:
	print('No input provided. Exiting..')
	exit()



fill_dict(gtf_file=args.input)
threshold=args.threshold

list_keys=[]
for k in gene_id_coords.keys():
	list_keys.append(k)

#iterate through records
list_of_dicts=[exons_records,cdss_records,five_prime_utrs_records,three_prime_utrs_records,start_stop_codons]
for idx, key in enumerate(list_keys):
	if idx==0:
		cur_gene_id=list_keys[idx]
		cur_gene_as_gtf=GTF_element(*gene_id_records[cur_gene_id].split('\t'))
		cur_gene_attributes=parse_attributes(cur_gene_as_gtf.attributes)
		#now check position: remember that edits must be done against cur gene
		#is it coding?
		if cur_gene_attributes['gene_biotype']=='"protein_coding"':
			if int(gene_id_coords[cur_gene_id][1]) > threshold:
				iterate_and_extend_rev(cur_gene_id,threshold)
			elif 1 < int(gene_id_coords[cur_gene_id][1]) < threshold:
				adj_threshold= int(gene_id_coords[cur_gene_id][1]) - 1
				iterate_and_extend_rev(cur_gene_id,threshold=threshold)
			elif int(gene_id_coords[cur_gene_id][1]) == 1:
				iterate_not_extend_rev(cur_gene_id)
		#this is the first record: just extend if possible
	if 0 < idx < len(list_keys) - 1:
		cur_gene_id=list_keys[idx]
		prev_gene_id=list_keys[idx-1]
		cur_gene_as_gtf=GTF_element(*gene_id_records[cur_gene_id].split('\t'))
		cur_gene_attributes=parse_attributes(cur_gene_as_gtf.attributes)
		#now check position: remember that edits must be done against cur gene
		#is it coding?
		if cur_gene_attributes['gene_biotype'] == '"protein_coding"':
			#first check if records are under the same chromosome
			if gene_id_coords[cur_gene_id][0] == gene_id_coords[prev_gene_id][0]:
				#same chromosome
				residue=int(gene_id_coords[cur_gene_id][1]) - int(gene_id_coords[prev_gene_id][2])
				if residue > threshold:
					iterate_and_extend_rev(cur_gene_id,threshold=threshold)
				elif 0 < residue <= threshold:
					adj_threshold=int(residue-1)
					iterate_and_extend_rev(cur_gene_id,threshold=adj_threshold)
				elif int(gene_id_coords[prev_gene_id][2]) > int(gene_id_coords[cur_gene_id][1]):
					iterate_not_extend_rev(cur_gene_id)
			else:
				#different chromosome: extend if possibile
				if int(gene_id_coords[cur_gene_id][1]) > threshold:
					iterate_and_extend_rev(cur_gene_id,threshold=threshold)
				elif 1 < int(gene_id_coords[cur_gene_id][1]) < threshold:
					adj_threshold= int(gene_id_coords[cur_gene_id][1]) - 1
					iterate_and_extend_rev(cur_gene_id,threshold=adj_threshold)
				elif int(gene_id_coords[cur_gene_id][1]) == 1:
					iterate_not_extend_rev(cur_gene_id)

		else:
			#non coding: just iterate
			iterate_not_extend_rev(cur_gene_id)
	elif idx == len(list_keys) - 1:
		#last record 
		cur_gene_id=list_keys[idx]
		cur_gene_as_gtf=GTF_element(*gene_id_records[cur_gene_id].split('\t'))
		cur_gene_attributes=parse_attributes(cur_gene_as_gtf.attributes)
		if cur_gene_attributes['gene_biotype']=='"protein_coding"':
			if int(gene_id_coords[cur_gene_id][1]) > threshold:
				iterate_and_extend_rev(cur_gene_id,threshold=threshold)
			elif int(gene_id_coords[cur_gene_id][1]) < threshold:
				adj_threshold= int(gene_id_coords[cur_gene_id][1]) - 1
				iterate_and_extend_rev(cur_gene_id,threshold=adj_threshold)
			elif int(gene_id_coords[cur_gene_id][1]) == 1:
				iterate_not_extend_rev(cur_gene_id)