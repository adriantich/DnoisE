# Copyright 2006-2013 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.SeqIO support module (not for general use).

Unless you are writing a new parser or writer for Bio.SeqIO, you should not
use this module.  It provides base classes to try and simplify things.
"""

from abc import ABC, abstractmethod
from Bio import StreamModeError


class SequenceIterator(ABC):
    """Base class for building SeqRecord iterators.

    You should write a parse method that returns a SeqRecord generator.  You
    may wish to redefine the __init__ method as well.
    """

    def __init__(self, source, alphabet=None, mode="t", fmt=None):
        """Create a SequenceIterator object.

        Arguments:
        - source - input file stream, or path to input file
        - alphabet - no longer used, should be None

        This method MAY be overridden by any subclass.

        Note when subclassing:
        - there should be a single non-optional argument, the source.
        - you do not have to require an alphabet.
        - you can add additional optional arguments.
        """
        if alphabet is not None:
            raise ValueError("The alphabet argument is no longer supported")
        try:
            self.stream = open(source, "r" + mode)
            self.should_close_stream = True
        except TypeError:  # not a path, assume we received a stream
            if mode == "t":
                if source.read(0) != "":
                    raise StreamModeError(
                        "%s files must be opened in text mode." % fmt
                    ) from None
            elif mode == "b":
                if source.read(0) != b"":
                    raise StreamModeError(
                        "%s files must be opened in binary mode." % fmt
                    ) from None
            else:
                raise ValueError("Unknown mode '%s'" % mode) from None
            self.stream = source
            self.should_close_stream = False
        try:
            self.records = self.parse(self.stream)
        except Exception:
            if self.should_close_stream:
                self.stream.close()
            raise

    def __next__(self):
        try:
            return next(self.records)
        except Exception:
            if self.should_close_stream:
                self.stream.close()
            raise

    def __iter__(self):
        """Iterate over the entries as a SeqRecord objects.

        Example usage for Fasta files::

            with open("example.fasta","r") as myFile:
                myFastaReader = FastaIterator(myFile)
                for record in myFastaReader:
                    print(record.id)
                    print(record.seq)

        This method SHOULD NOT be overridden by any subclass. It should be
        left as is, which will call the subclass implementation of __next__
        to actually parse the file.
        """
        return self

    @abstractmethod
    def parse(self, handle):
        """Start parsing the file, and return a SeqRecord iterator."""


