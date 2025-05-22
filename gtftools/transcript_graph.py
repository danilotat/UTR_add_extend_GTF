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


class TreeNode(object):
    """
    Basic entry for each gene ids

    Attributes
    ----------
    geneID : str
        The gene ID of the entry
    start : int
        The start position of the entry
    end : int
        The end position of the entry
    length : int
        The length of the entry

    Methods
    -------
    __init__(geneID, start, end)
        Initialize a TreeNode object.
    """

    def __init__(self, geneID, start, end) -> None:
        self.geneID = geneID
        self.start = start
        self.end = end
        self.length = abs(self.end - self.start)
    
    def __repr__(self):
        return f"{write_in_red('geneID')}: {self.geneID}\t{write_in_red('start')}: {self.start}\t{write_in_red('end')}: {self.end}\t \
               {write_in_red("Length")}: {self.length}"


class TreeBranch(object):
    """
    This class represents a branch of the tree, specific for a strand of a chromosome.
    Stores record in the form of TreeNode objects.

    Attributes
    ----------
    strand : str
        The strand that the tree branch belongs to.
    entries : dict
        A dictionary of entries in the tree branch.

    Methods
    -------
    __init__(strand)
        Initialize a TreeBranch object.
    _insert_record(record)
        Insert a TreeNode record into the entries.
    _sort_entries()
        Sort the entries based on the strand.
    _query_previous(geneID)
        Query the previous TreeNode based on the gene ID.
    _query_next(geneID)
        Query the next TreeNode based on the gene ID.
    """

    def __init__(self, strand) -> None:
        self.strand = strand
        self.entries = {}

    def _insert_record(self, record: TreeNode):
        self.entries[record.geneID] = (record.start, record.end)

    def _sort_entries(self):
        """
        If strand is fw, sort by start
        Else, sort by end
        """
        if self.strand == "+":
            self.entries = {
                k: v
                for k, v in sorted(self.entries.items(), key=lambda item: item[1][0])
            }
        if self.strand == "-":
            self.entries = {
                k: v
                for k, v in sorted(self.entries.items(), key=lambda item: item[1][1], reverse=True)
            }

    def _query_previous(self, geneID: str) -> TreeNode:
        """
        Given that the sorting is handled in different ways from forward and reverse,
        this method has the same behavior in both the strands

        Parameters
        ----------
        geneID : str
            The gene ID of the current TreeNode

        Returns
        -------
        TreeNode
            The previous TreeNode of the current TreeNode

        Raises
        ------
        ValueError
            If the gene ID is not found in the branch
        IndexError
            If the gene ID is the first TreeNode in the branch
        """
        ids = list(self.entries.keys())
        try:
            tid_idx = ids.index(geneID)
            try:
                prev_tid = ids[tid_idx - 1]
                return TreeNode(
                    prev_tid, self.entries[prev_tid][0], self.entries[prev_tid][1]
                )
            except IndexError:
                return None
        except ValueError:
            raise ValueError(f"Gene ID {geneID} not found in the branch")

    def _query_next(self, geneID: str) -> TreeNode:
        """
        Given that the sorting is handled in different ways from forward and reverse,
        this method has the same behavior in both the strands

        Parameters
        ----------
        geneID : str
            The gene ID of the current TreeNode

        Returns
        -------
        TreeNode
            The next TreeNode of the current TreeNode

        Raises
        ------
        ValueError
            If the gene ID is not found in the branch
        IndexError
            If the gene ID is the last TreeNode in the branch
        """
        ids = list(self.entries.keys())
        try:
            tid_idx = ids.index(geneID)
            try:
                next_tid = ids[tid_idx + 1]
                return TreeNode(
                    next_tid, self.entries[next_tid][0], self.entries[next_tid][1]
                )
            except IndexError:
                return None
        except ValueError:
            raise ValueError(f"Gene ID {geneID} not found in the branch")


class BST(object):
    """
    Class for the dependency graph of the GTF file.
    While the name was originally chosen because it'd use a binary search to retrieve previous and next item across the branch, now is just a collector.

    It stores the gene_IDs and tuples of (start, end)

    Attributes
    ----------
    root : TreeNode
        The root of the tree
    chromosome : str
        The chromosome that the tree belongs to
    strand : str
        The strand that the tree belongs to
    branch : TreeBranch
        The branch of the tree

    Methods
    -------
    __init__(chromosome, strand)
        Initialize a BST object.
    """

    def __init__(self, chromosome, strand) -> None:
        self.root = None
        self.chromosome = chromosome
        self.strand = strand
        self.branch = TreeBranch(self.strand)
