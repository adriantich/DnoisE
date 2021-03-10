
from Bio.SeqIO import FastaIO


# Convention for format names is "mainname-subtype" in lower case.
# Please use the same names as BioPerl or EMBOSS where possible.
#
# Note that this simple system copes with defining
# multiple possible iterators for a given format/extension
# with the -subtype suffix
#
# Most alignment file formats will be handled via Bio.AlignIO

_FormatToIterator = {
    "fasta": FastaIO.FastaIterator,
}






def parse(handle, format, alphabet=None):


    # Try and give helpful error messages:
    if not isinstance(format, str):
        raise TypeError("Need a string for the file format (lower case)")
    if not format:
        raise ValueError("Format required (lower case string)")
    if not format.islower():
        raise ValueError("Format string '%s' should be lower case" % format)
    if alphabet is not None:
        raise ValueError("The alphabet argument is no longer supported")

    iterator_generator = _FormatToIterator.get(format)
    if iterator_generator:
        return iterator_generator(handle)
    raise ValueError("Unknown format '%s'" % format)

