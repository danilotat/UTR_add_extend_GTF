#!/usr/bin/python3

def write_in_red(text: str):
    """
    Print the text in red color.

    Parameters
    ----------
    text : str
        The text to be printed in red color.
    """
    return f"\033[91m{text}\033[0m"


class GTF_record(object):
    """
    A GTF record is the first building block of the parser.
    The attribute field is parsed resulting in a dict.

    Attributes
    ----------
    chromosome : str
        The chromosome that the GTF record belongs to.
    source : str
        The source of the GTF record.
    feature_type : str
        The type of feature that the GTF record represents.
    start : int
        The start position of the feature in the chromosome.
    end : int
        The end position of the feature in the chromosome.
    score : str
        The score of the GTF record.
    strand : str
        The strand that the GTF record belongs to.
    phase : str
        The phase of the GTF record.
    length : int
        The length of the feature.
    attributes : dict
        A dictionary of attributes from the GTF record.

    Methods
    -------
    __init__(chromosome, source, feature_type, start, end, score, strand, phase, attributes)
        Initialize a GTF record object.
    parse_attributes(attributes)
        Parse the attributes of a GTF record.
    is_coding(feat_dict)
        Check if a GTF record is coding.
    """

    def __init__(
        self,
        chromosome,
        source,
        feature_type,
        start,
        end,
        score,
        strand,
        phase,
        attributes,
    ):
        self.chromosome = str(chromosome)
        self.source = source
        self.feature_type = feature_type
        self.start = int(start)
        self.end = int(end)
        self.score = score
        self.strand = strand
        self.phase = phase
        self.length = abs(self.end - self.start)
        self.attributes = GTF_record.parse_attributes(attributes)

    @staticmethod
    def parse_attributes(attributes):
        if isinstance(attributes, dict):
            return attributes
        else:
            feat_dict = {}
            keyVal = attributes.split(";")[:-1]
            for item in keyVal:
                replItem = item.replace(' "', '="')
                # populate the dict
                try:
                    k, v = replItem.split("=")
                    rk = k.replace(" ", "")
                    vk = v.replace('"', "")
                    feat_dict[rk] = vk
                except ValueError:
                    print(f"Unable to parse this attribute field\n{attributes}")
                    exit()
            return feat_dict

    @staticmethod
    def is_coding(feat_dict: dict) -> bool:
        """
        Simple check if the given record is protein coding.

        Parameters
        ----------
        feat_dict : dict
            A dictionary of attributes from the GTF record.

        Returns
        -------
        bool
            True if the record is protein coding, False otherwise.
        """
        if feat_dict["gene_biotype"] == "protein_coding":
            return True
        else:
            return False

    def __repr__(self):
        return f"{write_in_red("Feature type")}: {self.feature_type}\t{write_in_red("chr")}:{self.chromosome}\t{write_in_red("start")}:{self.start}\t{write_in_red("end")}:{self.end}\t{write_in_red("strand")}:({self.strand})\t{write_in_red("Length")}:{self.length}\t" \
            f"{write_in_red('Attributes')}: {self.attributes}\n" 

    def asStr(self):
        """
        Return the record as a string for writing out.

        Returns
        -------
        str
            The record as a string, separated by tabs.
        """
        attributesAsStr = "; ".join([f'{k} "{v}"' for k, v in self.attributes.items()])
        return (
            "\t".join(
                [
                    str(x)
                    for x in [
                        self.chromosome,
                        self.source,
                        self.feature_type,
                        self.start,
                        self.end,
                        self.score,
                        self.strand,
                        self.phase,
                        attributesAsStr,
                    ]
                ]
            )
            + "\n"
        )

    def _update_length(self, strand, value):
        """
        Update the length of the transcript based on the 3'UTR extension

        Parameters
        ----------
        strand : str
            The strand that the transcript belongs to.
        value : int
            The value to update the terminal with. It can be either the start or the
            end, it'll be handled accordingly to the strand.

        Raises
        ------
        ValueError
            If the strand is not valid.
        """
        if strand == "+":
            self.end = value
            self.length = abs(self.end - self.start)
        elif strand == "-":
            self.start = value
            self.length = abs(self.end - self.start)
        else:
            raise ValueError("The strand is not valid. Please use either '+' or '-'")


class gene(GTF_record):
    """
    A gene is a GTF_record with a list of transcripts associated

    Attributes
    ----------
    attributes : dict
        A dictionary of attributes from the GTF record.
    iscoding : bool
        A boolean indicating whether the gene is coding or not.
    transcripts : list
        A list of transcript IDs associated with the gene.

    Methods
    -------
    _add_transcript(transcript)
        Add a transcript ID to the gene.
    _longest_transcript(transcripts)
        Return the ID of the longest transcript associated with the gene.
    """

    def __init__(self, gtf_record):
        super().__init__(
            gtf_record.chromosome,
            gtf_record.source,
            gtf_record.feature_type,
            gtf_record.start,
            gtf_record.end,
            gtf_record.score,
            gtf_record.strand,
            gtf_record.phase,
            gtf_record.attributes,
        )
        self.attributes = gtf_record.attributes
        self.iscoding = GTF_record.is_coding(self.attributes)
        self.transcripts = []

    def _add_transcript(self, transcript: str):
        self.transcripts.append(transcript)

    def _longest_transcript(self, transcripts: dict) -> str:
        """
        Return the longest transcript for this gene.

        Parameters
        ----------
        transcripts : dict
            A dictionary where the key is a transcript ID and the value is a transcript object.

        Returns
        -------
        str
            The transcript_ID of the longest transcript.
        """
        # get lengths from the treeBranch
        lengths = {}
        for tid in self.transcripts:
            transcript = transcripts[tid]
            lengths[tid] = transcript.length
        # sort by length
        lengths = {
            k: v for k, v in sorted(lengths.items(), key=lambda x: x[1], reverse=True)
        }
        return list(lengths.keys())[0]


class transcript(GTF_record):
    """
    A transcript is a key class of the GTF parser. It has multiple elements associated with it including exons, CDSs, 3'UTR, 5'UTR, start-codon, stop-codon.
    In our case it acts as a collector, under which you could collect all the associated features.

    This class inherits from the GTF_record class and adds additional attributes and methods specific to transcripts. It acts as a collector for all associated features of a transcript, including exons, CDSs, 3'UTRs, 5'UTRs, start codons, and stop codons.

    Attributes
    ----------
    gene_id : str
        The ID of the gene that the transcript belongs to.
    transcript_id : str
        The ID of the transcript.
    exons : list
        A list of exons associated with the transcript.
    cdss : list
        A list of CDSs associated with the transcript.
    threeprimeutrs : list
        A list of 3' UTRs associated with the transcript.
    fiveprimeutrs : list
        A list of 5' UTRs associated with the transcript.
    startcodon : list
        A list of start codons associated with the transcript.
    stopcodon : list
        A list of stop codons associated with the transcript.

    Methods
    -------
    _add_exon(exon)
        Add an exon to the transcript.
    _add_cds(cds)
        Add a CDS to the transcript.
    _add_threeprimeutr(threeprimeutr)
        Add a 3' UTR to the transcript.
    _add_fiveprimeutr(fiveprimeutr)
        Add a 5' UTR to the transcript.
    _add_startcodon(startcodon)
        Add a start codon to the transcript.
    _add_stopcodon(stopcodon)
        Add a stop codon to the transcript.
    _extend_three_prime_utr(geneEntry, next_record, extension_mode, min_distance=10)
        Extend the 3' UTR of the transcript based on the provided extension mode.
    """

    def __init__(self, gtf_record):
        super().__init__(
            gtf_record.chromosome,
            gtf_record.source,
            gtf_record.feature_type,
            gtf_record.start,
            gtf_record.end,
            gtf_record.score,
            gtf_record.strand,
            gtf_record.phase,
            gtf_record.attributes,
        )
        self.gene_id = self.attributes["gene_id"]
        self.transcript_id = self.attributes["transcript_id"]
        self.exons = []
        self.cdss = []
        self.threeprimeutrs = []
        self.fiveprimeutrs = []
        self.startcodon = []
        self.stopcodon = []

    def _add_exon(self, exon: str):
        self.exons.append(exon)

    def _add_cds(self, cds: str):
        self.cdss.append(cds)

    def _add_threeprimeutr(self, threeprimeutr: GTF_record):
        self.threeprimeutrs.append(threeprimeutr)

    def _add_fiveprimeutr(self, fiveprimeutr: str):
        self.fiveprimeutrs.append(fiveprimeutr)

    def _add_startcodon(self, startcodon: str):
        self.startcodon.append(startcodon)

    def _add_stopcodon(self, stopcodon: str):
        self.stopcodon.append(stopcodon)

    def _extend_three_prime_utr(
        self, geneEntry: gene, next_record, extension_mode, min_distance=10
    ):
        """
        Knowing the downstream record, the 3' UTR could be then extended following
        one of the parsable rules, passed via the `extension_mode` params.

        The extension mode could be:

        - "max": the 3' UTR is generated starting from the end of the last CDS of the actual transcript till the start of the next record,
        - "min": the 3' UTR is generated starting from the end of the last CDS of the actual transcript till the end of the transcript.
        - int: passing an integer will extend the 3'UTR with a fixed length starting from the end of the last CDS. If the passed integer is higher than the distance with the next transcript, it'll instead extended using the distance between the records minus the min_distance.

        REMIND: the extension mode is different for each strand.

        Parameters
        ----------
        geneEntry : gene
            A gene object representing the gene that the transcript belongs to.
        next_record : GTF_record
            A GTF record object representing the next record in the GTF file.
        extension_mode : str or int
            The mode of extension. Can be 'max', 'min', or an integer representing the fixed length of extension.
        min_distance : int, optional
            The minimum distance between the end of the 3' UTR and the start of the next record. Default is 10.

        Returns
        -------
        tuple
            A tuple where the first element is the transcript ID and the second element is the length of the extension.

        Raises
        ------
        ValueError
            If the extension mode is not 'max', 'min', or an integer.
        """
        if not isinstance(extension_mode, int):
            if not extension_mode in ["max", "min"]:
                raise ValueError(
                    "The extension mode for the 3' UTR is not valid. Please use one of the following: 'max', 'min' or an integer."
                )
        # check that there's at least a CDS
        edited_len = 0
        if len(self.cdss) != 0:
            if self.strand == "+":
                # given that we're already checking the transcript whose end is the same of the gene,
                # we could safely initialize both at the same value. It'll be adjusted later, if needed.
                maximum_space = next_record.start - (self.cdss[-1].end + 4) # as we need to account for the stop codon.
                if maximum_space < min_distance:
                    return (self.attributes["transcript_id"], 0)
                else:
                    if extension_mode == "max":
                        utr_start = self.cdss[-1].end + 4
                        utr_end = self.cdss[-1].end + 4 + maximum_space - min_distance
                    elif isinstance(extension_mode, int):
                        if extension_mode > maximum_space:
                            # the requested distance is higher than the available "space" on the chromosome.
                            # add just the minimum distance as requested.
                            # get the last cds
                            utr_start = self.cdss[-1].end + 4
                            utr_end = next_record.start - min_distance
                        else:
                            # we have enough space.
                            utr_start = self.cdss[-1].end + 4
                            utr_end = self.cdss[-1].end + 4 + extension_mode
                    # adjust also the transcript length based on that.
                    edited_len = abs(utr_end - utr_start)
                    self._update_length(self.strand, utr_end)
                    geneEntry.end = utr_end
                    # initialize a GTF record and assign the novel 3'UTR only if the 3'UTR is not already present.
                    if len(self.threeprimeutrs) == 0:
                        extension_length = abs(utr_end - utr_start)
                        ghostFeat = GTF_record(
                            self.chromosome,
                            self.source,
                            "three_prime_utr",
                            utr_start,
                            utr_end,
                            self.score,
                            self.strand,
                            self.phase,
                            self.attributes,
                        )
                        self._add_threeprimeutr(ghostFeat)
                    else:
                        extension_length = abs(
                            self.threeprimeutrs[-1].end - utr_end
                        )
                        self.threeprimeutrs[-1].start = utr_start
                        self.threeprimeutrs[-1].end = utr_end
                        self.threeprimeutrs[-1].length = abs(
                            self.threeprimeutrs[-1].end
                            - self.threeprimeutrs[-1].start
                        )
                    return (self.attributes["transcript_id"], extension_length)
            elif self.strand == "-":
                maximum_space = self.cdss[-1].start - next_record.end
                if maximum_space < min_distance:
                    return (self.attributes["transcript_id"], 0)
                else:
                    if extension_mode == "max":
                        utr_start = next_record.end + min_distance
                        utr_end = self.cdss[-1].start
                    if isinstance(extension_mode, int):
                        if extension_mode > (
                            (self.cdss[-1].start + min_distance) - next_record.end
                        ):
                            utr_start = next_record.end + min_distance
                            utr_end = self.cdss[-1].start
                        else:
                            utr_start = self.cdss[-1].start - extension_mode
                            utr_end = self.cdss[-1].start
                    edited_len = utr_start - self.cdss[-1].start
                    # adjust the transcript coord
                    self._update_length(self.strand, utr_start)
                    geneEntry.start = utr_start
                    ghostFeat = GTF_record(
                        self.chromosome,
                        self.source,
                        "three_prime_utr",
                        utr_start,
                        utr_end,
                        self.score,
                        self.strand,
                        self.phase,
                        self.attributes,
                    )
                    self._add_threeprimeutr(ghostFeat)
                    return (self.attributes["transcript_id"], edited_len)
        else:
            return (self.attributes["transcript_id"], edited_len)


class exon(GTF_record):
    """
    An exon is a simpler entry.

    Attributes
    ----------
    gene_id : str
        The ID of the gene that the exon belongs to.
    transcript_id : str
        The ID of the transcript that the exon belongs to.
    exon_id : str
        The ID of the exon. This attribute is optional and its value may be None if the 'exon_id' attribute is not present in the GTF record.
    """

    def __init__(self, gtf_record):
        super().__init__(
            gtf_record.chromosome,
            gtf_record.source,
            gtf_record.feature_type,
            gtf_record.start,
            gtf_record.end,
            gtf_record.score,
            gtf_record.strand,
            gtf_record.phase,
            gtf_record.attributes,
        )
        self.gene_id = self.attributes["gene_id"]
        self.transcript_id = self.attributes["transcript_id"]
        self.exon_id = self.attributes.get("exon_id")


class cds(GTF_record):
    """
    A coding sequence.

    Attributes
    ----------
    gene_id : str
        The ID of the gene that the exon belongs to.
    transcript_id : str
        The ID of the transcript that the exon belongs to.
    exon_id : str
        The ID of the exon. This attribute is optional and its value may be None if the 'exon_id' attribute is not present in the GTF record.
    """

    def __init__(self, gtf_record):
        super().__init__(
            gtf_record.chromosome,
            gtf_record.source,
            gtf_record.feature_type,
            gtf_record.start,
            gtf_record.end,
            gtf_record.score,
            gtf_record.strand,
            gtf_record.phase,
            gtf_record.attributes,
        )
        self.gene_id = self.attributes["gene_id"]
        self.transcript_id = self.attributes["transcript_id"]
        self.exon_id = self.attributes["exon_id"]


class threeprimeutr(GTF_record):
    """
    The three-prime UTR is the main focus here.

    Attributes
    ----------
    gene_id : str
        The ID of the gene that the exon belongs to.
    transcript_id : str
        The ID of the transcript that the exon belongs to.
    exon_id : str
        The ID of the exon. This attribute is optional and its value may be None if the 'exon_id' attribute is not present in the GTF record.
    """

    def __init__(self, gtf_record):
        super().__init__(
            gtf_record.chromosome,
            gtf_record.source,
            gtf_record.feature_type,
            gtf_record.start,
            gtf_record.end,
            gtf_record.score,
            gtf_record.strand,
            gtf_record.phase,
            gtf_record.attributes,
        )
        self.gene_id = self.attributes["gene_id"]
        self.transcript_id = self.attributes["transcript_id"]
        self.exon_id = self.attributes.get("exon_id")


class fiveprimeutr(GTF_record):
    """
    Simple class for 5' UTR.

    Attributes
    ----------
    gene_id : str
        The ID of the gene that the exon belongs to.
    transcript_id : str
        The ID of the transcript that the exon belongs to.
    exon_id : str
        The ID of the exon. This attribute is optional and its value may be None if the 'exon_id' attribute is not present in the GTF record.
    """

    def __init__(self, gtf_record):
        super().__init__(
            gtf_record.chromosome,
            gtf_record.source,
            gtf_record.feature_type,
            gtf_record.start,
            gtf_record.end,
            gtf_record.score,
            gtf_record.strand,
            gtf_record.phase,
            gtf_record.attributes,
        )
        self.gene_id = self.attributes["gene_id"]
        self.transcript_id = self.attributes["transcript_id"]
        self.exon_id = self.attributes["exon_id"]
