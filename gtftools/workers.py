#!/usr/bin/python3

import itertools
from .objects import *
from .transcript_graph import BST, TreeNode


def other_entries(entry: GTF_record, transcript_collection: dict):
    """
    Append different types of entries to the corresponding transcript in the collection.

    This function takes a GTF record and a collection of transcripts. Depending on the feature type of the GTF record, it appends the record to the corresponding transcript in the collection. The function handles the following feature types: 'exon', 'CDS', 'five_prime_utr', 'three_prime_utr', 'start_codon', 'stop_codon'.

    Parameters
    ----------
    entry : GTF_record
        A GTF record object representing an entry in a GTF file.
    transcript_collection : dict
        A dictionary where the key is a transcript ID and the value is a transcript object.

    Returns
    -------
    None

    Raises
    ------
    KeyError
        If the 'transcript_id' attribute is not found in the GTF record.
    """
    try:
        transcript_id = entry.attributes["transcript_id"]
    except KeyError:
        print(entry.asStr())
        exit()
    if entry.feature_type == "exon":
        transcript_collection[transcript_id]._add_exon(entry)
    if entry.feature_type == "CDS":
        transcript_collection[transcript_id]._add_cds(entry)
    if entry.feature_type == "five_prime_utr":
        transcript_collection[transcript_id]._add_fiveprimeutr(entry)
    if entry.feature_type == "three_prime_utr":
        transcript_collection[transcript_id]._add_threeprimeutr(entry)
    if entry.feature_type == "start_codon":
        transcript_collection[transcript_id]._add_startcodon(entry)
    if entry.feature_type == "stop_codon":
        transcript_collection[transcript_id]._add_stopcodon(entry)


def first_iteration(gtf: str) -> dict:
    """
    The first iteration is used just to create a dependancy graph of the GTF file, using the IDs of the transcripts.
    This is used to know which record comes before and after each item, and which transcript has missing UTR

    Parameters
    ----------
    gtf : str
        The path to the GTF file to be parsed.

    Returns
    -------
    tuple
        A tuple containing four elements:
        - headers: A list of header lines from the GTF file.
        - trees: A dictionary where the keys are chromosome IDs and the values are dictionaries mapping strand IDs to BST objects.
        - genes: A dictionary where the keys are gene IDs and the values are gene objects.
        - transcripts: A dictionary where the keys are transcript IDs and the values are transcript objects.

    Raises
    ------
    IndexError
        If a transcript record is encountered before its corresponding gene record in the GTF file.
    """
    headers = []
    trees = {}
    genes = {}
    transcripts = {}
    with open(gtf, "r") as gtf_file:
        for line in gtf_file:
            if line.startswith("#"):
                headers.append(line)
            else:
                entry = GTF_record(*line.rstrip().split("\t"))
                if entry.feature_type == "gene":
                    if entry.chromosome in trees.keys():
                        if entry.strand in trees[entry.chromosome].keys():
                            # the tree already exists
                            trees[entry.chromosome][entry.strand].branch._insert_record(
                                TreeNode(
                                    entry.attributes["gene_id"], entry.start, entry.end
                                )
                            )
                        else:
                            # strand is seen for the first time
                            trees[entry.chromosome][entry.strand] = BST(
                                entry.chromosome, entry.strand
                            )
                            trees[entry.chromosome][entry.strand].branch._insert_record(
                                TreeNode(
                                    entry.attributes["gene_id"], entry.start, entry.end
                                )
                            )
                    else:
                        # first time seen this chromosome
                        trees[entry.chromosome] = {}
                        trees[entry.chromosome][entry.strand] = BST(
                            entry.chromosome, entry.strand
                        )
                        trees[entry.chromosome][entry.strand].branch._insert_record(
                            TreeNode(
                                entry.attributes["gene_id"], entry.start, entry.end
                            )
                        )
                    if entry.attributes["gene_id"] not in genes:
                        genes[entry.attributes["gene_id"]] = gene(entry)
                        continue
                if entry.feature_type == "transcript":
                    transcripts[entry.attributes["transcript_id"]] = transcript(entry)
                    # populate the genes dictionary
                    if entry.attributes["gene_id"] in genes.keys():
                        genes[entry.attributes["gene_id"]]._add_transcript(
                            entry.attributes["transcript_id"]
                        )
                    else:
                        raise IndexError("Malformed gtf")
                else:
                    other_entries(entry, transcripts)
    # sort the trees
    for chromosome in trees:
        for strand in trees[chromosome]:
            trees[chromosome][strand].branch._sort_entries()
    return headers, trees, genes, transcripts


def add_utr(genes: dict, transcripts: dict):
    """
    This function does the second iteration.
    Add just UTRs to the transcripts

    Parameters
    ----------
    genes : dict
        A dictionary where the key is a gene ID and the value is a gene object.
    transcripts : dict
        A dictionary where the key is a transcript ID and the value is a transcript object.

    Returns
    -------
    dict
        A dictionary where the key is a transcript ID and the value is the length of the added 3' UTR. Transcripts that did not have a 3' UTR added are not included in the dictionary.
    """
    added_utr = {}
    for gene in genes:
        # the idea is that:
        # is coding?
        geneEntry = genes[gene]
        if geneEntry.iscoding:
            # cool: let's see if there's the 3'UTR
            # iterate through transcripts
            for tid in geneEntry.transcripts:
                transcript = transcripts[tid]
                if len(transcript.threeprimeutrs) == 0:
                    # no UTR
                    if transcript.strand == "+":
                        # is forward
                        if len(transcript.stopcodon) > 0:
                            # there's a stop codon
                            if transcript.end > transcript.stopcodon[0].end:
                                # the last exon goes beyond the stop codon
                                # let's add the UTR
                                transcript._add_threeprimeutr(
                                    GTF_record(
                                        transcript.exons[-1].chromosome,
                                        transcript.source,
                                        "three_prime_utr",
                                        transcript.stopcodon[0].end,
                                        transcript.end,
                                        transcript.score,
                                        transcript.strand,
                                        transcript.phase,
                                        transcript.attributes,
                                    )
                                )
                                added_utr[tid] = transcript.threeprimeutrs[-1].length
                    else:
                        # is reverse
                        if len(transcript.stopcodon) > 0:
                            # there's a stop codon
                            if transcript.start < transcript.stopcodon[0].start:
                                # the last exon goes beyond the stop codon
                                # let's add the UTR
                                transcript._add_threeprimeutr(
                                    GTF_record(
                                        transcript.exons[-1].chromosome,
                                        transcript.source,
                                        "three_prime_utr",
                                        transcript.start,
                                        transcript.stopcodon[0].start,
                                        transcript.score,
                                        transcript.strand,
                                        transcript.phase,
                                        transcript.attributes,
                                    )
                                )
                                added_utr[tid] = transcript.threeprimeutrs[-1].length
    return added_utr


def extend_utr(genes, transcripts, trees, extension_mode, min_dist):
    """
    This adds extended UTRs to the given transcript and adjusts the coords
    of the source gene and transcript.

    Parameters
    ----------
    genes : dict
        A dictionary where the key is a gene ID and the value is a gene object.
    transcripts : dict
        A dictionary where the key is a transcript ID and the value is a transcript object.
    trees : dict
        A dictionary where the key is a chromosome ID and the value is a tree object.
    extension_mode : int or str
        The mode of extension. If an integer is provided, it is used as the length of the extension. If 'max' is provided, the maximum possible extension is used.
    min_dist : int
        The minimum distance between the end of the transcript and the start of the next gene.

    Returns
    -------
    dict
        A dictionary where the key is a transcript ID and the value is the length of the extension. Transcripts that were not extended have a length of 0.

    Raises
    ------
    AttributeError
        If `extension_mode` is neither an integer nor 'max'.
    """
    if str(extension_mode) != "max":
        try:
            extension_mode = int(extension_mode)
        except ValueError:
            raise AttributeError('extension mode accepts only an integer or "max"')
    # iterate through genes.
    logged_infos = {}
    for gene in genes:
        # UTR only on protein coding genes
        geneEntry = genes[gene]
        if geneEntry.iscoding:
            for tid in geneEntry.transcripts:
                # IMPORTANT:
                # we want to extend only the transcript which ends at the same position of the gene.
                transcriptEntry = transcripts[tid]
                if (
                    transcriptEntry.end == geneEntry.end
                    and not tid in logged_infos.keys()
                ):
                    # get the next record. we need the right tree to select from
                    nextRecord = trees[transcriptEntry.chromosome][
                        transcriptEntry.strand
                    ].branch._query_next(gene)
                    if nextRecord is not None:
                        # try the extension
                        txid, ext_len = transcriptEntry._extend_three_prime_utr(
                            geneEntry, nextRecord, extension_mode, min_dist
                        )
                        if txid is not None:
                            logged_infos[txid] = ext_len
                    else:
                        continue
                else:
                    # the transcript ends before the gene end.
                    # just keep iterating.
                    continue
    return logged_infos


def recompose_gtf(genes, transcripts, out):
    """
    This function recompose the final gtf.

    Parameters
    ----------
    genes : dict
        A dictionary where the key is a gene ID and the value is a gene object. Generated
        with the `first_iteration` function.
    transcripts : dict
        A dictionary where the key is a transcript ID and the value is a transcript object.
        Generated with the `first_iteration` function.
    outfile : '_io.TextIOWrapper'
        The output file where the GTF will be written.

    Returns
    -------
    None
    """
    for gene in genes:
        geneEntry = genes[gene]
        out.write(geneEntry.asStr())
        for tid in geneEntry.transcripts:
            transcriptEntry = transcripts[tid]
            out.write(transcriptEntry.asStr())
            for utr in transcriptEntry.fiveprimeutrs:
                out.write(utr.asStr())
            for startcodon in transcriptEntry.startcodon:
                out.write(startcodon.asStr())
            for exon in transcriptEntry.exons:
                out.write(exon.asStr())
            for cds in transcriptEntry.cdss:
                out.write(cds.asStr())
            for utr in transcriptEntry.threeprimeutrs:
                out.write(utr.asStr())
            for stopcodon in transcriptEntry.stopcodon:
                out.write(stopcodon.asStr())


def chromosome_wise(gtf, outfile, extension_mode, min_distance, logfile=None):
    """
    This is a worker designed to parse the GTF file and create objects chromosome-wise, to reduce
    the memory footprint.

    Parameters
    ----------
    gtf : str
        The path to the GTF file to be parsed.
    outfile : str
        The path to the output file where the GTF file will be written.
    extension_mode : str or int
        The mode of extension. Can be 'max', 'min', or an integer representing the fixed length of extension.
    min_distance : int, optional
        The minimum distance between the end of the 3' UTR and the start of the next record. Default is 10.
    logfile : str, optional
        The path to the log file where the extension statistics will be written.
    """
    whole_logging = {}
    outfile = open(outfile, "w")
    print("--" * 20)
    print(f"Starting to parse {gtf}")
    print("--" * 20)
    with open(gtf, "r") as gtf_file:
        get_line = ((x.split("\t")[0], x) for x in gtf_file)
        for chromosome, lines in itertools.groupby(get_line, lambda x: x[0]):
            if chromosome.startswith("#"):
                outfile.writelines([x[1] for x in lines])
            else:
                print(f"Processing chromosome {chromosome}")
                trees = {}
                trees[chromosome] = {}
                genes = {}
                transcripts = {}
                for _, line in lines:
                    entry = GTF_record(*line.rstrip().split("\t"))
                    if entry.feature_type == "gene":
                        if entry.strand in trees[entry.chromosome].keys():
                            # the tree already exists
                            trees[entry.chromosome][entry.strand].branch._insert_record(
                                TreeNode(
                                    entry.attributes["gene_id"], entry.start, entry.end
                                )
                            )
                        else:
                            # strand is seen for the first time
                            trees[entry.chromosome][entry.strand] = BST(
                                entry.chromosome, entry.strand
                            )
                            trees[entry.chromosome][entry.strand].branch._insert_record(
                                TreeNode(
                                    entry.attributes["gene_id"], entry.start, entry.end
                                )
                            )
                        if entry.attributes["gene_id"] not in genes:
                            genes[entry.attributes["gene_id"]] = gene(entry)
                            continue
                    if entry.feature_type == "transcript":
                        transcripts[entry.attributes["transcript_id"]] = transcript(
                            entry
                        )
                        # populate the genes dictionary
                        if entry.attributes["gene_id"] in genes.keys():
                            genes[entry.attributes["gene_id"]]._add_transcript(
                                entry.attributes["transcript_id"]
                            )
                        else:
                            raise IndexError("Malformed gtf")
                    else:
                        other_entries(entry, transcripts)
                # sort the trees
                for strand in trees[chromosome]:
                    trees[chromosome][strand].branch._sort_entries()
                # now the relevant data structures are ready
                # let's do the extension
                added_utr = add_utr(genes, transcripts)
                print(f"Added {len(added_utr)} 3' UTRs")
                print("--" * 20)
                logged_infos = extend_utr(
                    genes, transcripts, trees, extension_mode, min_distance
                )
                # add the logged infos to the whole logging
                for geneid in logged_infos:
                    whole_logging[geneid] = logged_infos[geneid]
                recompose_gtf(genes, transcripts, outfile)
    print(f"Done.")
    print(
        f"Successfully extended {len([x for x in whole_logging if whole_logging[x] > 0])} transcripts out of {len(whole_logging)} transcripts."
    )
    outfile.close()
    if logfile:
        with open(logfile, "w") as f:
            for geneid in whole_logging:
                f.write(f"{geneid}\t{whole_logging[geneid]}\n")
