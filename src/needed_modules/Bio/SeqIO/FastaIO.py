# Copyright 2006-2017,2020 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
#
# This module is for reading and writing FASTA format files as SeqRecord
# objects.  The code is partly inspired  by earlier Biopython modules,
# Bio.Fasta.* and the now deprecated Bio.SeqIO.FASTA

"""Bio.SeqIO support for the "fasta" (aka FastA or Pearson) file format.

You are expected to use this module via the Bio.SeqIO functions.
"""


from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from .Interfaces import SequenceIterator


def SimpleFastaParser(handle):

    # Skip any text before the first record (e.g. blank lines, comments)
    for line in handle:
        if line[0] == ">":
            title = line[1:].rstrip()
            break
    else:
        # no break encountered - probably an empty file
        return

    # Main logic
    # Note, remove trailing whitespace, and any internal spaces
    # (and any embedded \r which are possible in mangled files
    # when not opened in universal read lines mode)
    lines = []
    for line in handle:
        if line[0] == ">":
            yield title, "".join(lines).replace(" ", "").replace("\r", "")
            lines = []
            title = line[1:].rstrip()
            continue
        lines.append(line.rstrip())

    yield title, "".join(lines).replace(" ", "").replace("\r", "")





class FastaIterator(SequenceIterator):
    """Parser for Fasta files."""

    def __init__(self, source, alphabet=None, title2ids=None):
        if alphabet is not None:
            raise ValueError("The alphabet argument is no longer supported")
        self.title2ids = title2ids
        super().__init__(source, mode="t", fmt="Fasta")

    def parse(self, handle):
        """Start parsing the file, and return a SeqRecord generator."""
        records = self.iterate(handle)
        return records

    def iterate(self, handle):
        """Parse the file and generate SeqRecord objects."""
        title2ids = self.title2ids
        if title2ids:
            for title, sequence in SimpleFastaParser(handle):
                id, name, descr = title2ids(title)
                yield SeqRecord(Seq(sequence), id=id, name=name, description=descr)
        else:
            for title, sequence in SimpleFastaParser(handle):
                try:
                    first_word = title.split(None, 1)[0]
                except IndexError:
                    assert not title, repr(title)
                    # Should we use SeqRecord default for no ID?
                    first_word = ""
                yield SeqRecord(
                    Seq(sequence), id=first_word, name=first_word, description=title,
                )

