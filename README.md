# UTR_add_extend_GTF

Provided tool will add explicit 3'UTR and extends it into a valid Ensembl GTF file. 
This seems to be useful while working with 3'RNAseq.


The repository is composed by 3 script file:

**gtf_advanced_parser**

Yet another gtf parser. This is intended to be more close to the biological meaning of a GTF file, so it works by creating OrderedDict of features by exploiting gene_id or transcript_id as indexes. So *genes*, *transcripts* and miscellaneous will be stored under the *gene_id* index; any other feature, like *exon*, *CDS*, *UTR* will be stored under the *transcript_id* index. If this script is executed as a standalone python script, it will add explicit 3'UTR and reformat the provided GTF by this scheme:
For any gene:

- gene
- transcript
- exons
- CDS 
- 5' UTR
- 3' UTR
- start & stop codons
- miscellaneous


**Extend_by_threshold_{fw-rev}.py**

Using the scheme above, this script use a {fw-rev} only version of a GTF* to extend the 3' terminus of genes by a given threshold when possible. By possible it means that no overlap with exisisting genes will be added and, if an overlap already exists, it won't be extended. Obviously no coding sequence will be extended, so this could be done only into the untranslated region with 3' coordinates.
It works by extending a gene, its transcript(s) with same end coordinates, its exon(s) with same end coordinates and its 3'UTR(s).


*(please refer to the bash script to see how to automatically split the GTF and pass results to these two scripts)*

**process_them**


Extremely sample bash script provided as a POC to see how to extend multiple times a valid GTF file using multiple threshold. Here a valid GTF is considered a file already processed with **gtf_advanced_parser**, like:

	$ python3 gtf_advanced_parser.py --input Input.gtf --output Output.gtf





