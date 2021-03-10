# Copyright 2000-2002 Andrew Dalke.  All rights reserved.
# Copyright 2002-2004 Brad Chapman.  All rights reserved.
# Copyright 2006-2020 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Represent a Sequence Record, a sequence with annotation."""


# NEEDS TO BE SYNCH WITH THE REST OF BIOPYTHON AND BIOPERL
# In particular, the SeqRecord and BioSQL.BioSeq.DBSeqRecord classes
# need to be in sync (this is the BioSQL "Database SeqRecord", see
# also BioSQL.BioSeq.DBSeq which is the "Database Seq" class)


from io import StringIO
from Bio import StreamModeError


_NO_SEQRECORD_COMPARISON = "SeqRecord comparison is deliberately not implemented. Explicitly compare the attributes of interest."


class _RestrictedDict(dict):

    def __init__(self, length):
        """Create an EMPTY restricted dictionary."""
        dict.__init__(self)
        self._length = int(length)

    def __setitem__(self, key, value):
        # The check hasattr(self, "_length") is to cope with pickle protocol 2
        # I couldn't seem to avoid this with __getstate__ and __setstate__
        if (
            not hasattr(value, "__len__")
            or not hasattr(value, "__getitem__")
            or (hasattr(self, "_length") and len(value) != self._length)
        ):
            raise TypeError(
                "We only allow python sequences (lists, tuples or strings) "
                f"of length {self._length}."
            )
        dict.__setitem__(self, key, value)

    def update(self, new_dict):
        # Force this to go via our strict __setitem__ method
        for (key, value) in new_dict.items():
            self[key] = value


class SeqRecord:

    def __init__(
        self,
        seq,
        id="<unknown id>",
        name="<unknown name>",
        description="<unknown description>",
        dbxrefs=None,
        features=None,
        annotations=None,
        letter_annotations=None,
    ):

        if id is not None and not isinstance(id, str):
            # Lots of existing code uses id=None... this may be a bad idea.
            raise TypeError("id argument should be a string")
        if not isinstance(name, str):
            raise TypeError("name argument should be a string")
        if not isinstance(description, str):
            raise TypeError("description argument should be a string")
        self._seq = seq
        self.id = id
        self.name = name
        self.description = description

        # database cross references (for the whole sequence)
        if dbxrefs is None:
            dbxrefs = []
        elif not isinstance(dbxrefs, list):
            raise TypeError("dbxrefs argument should be a list (of strings)")
        self.dbxrefs = dbxrefs

        # annotations about the whole sequence
        if annotations is None:
            annotations = {}
        elif not isinstance(annotations, dict):
            raise TypeError("annotations argument should be a dict")
        self.annotations = annotations

        if letter_annotations is None:
            # annotations about each letter in the sequence
            if seq is None:
                # Should we allow this and use a normal unrestricted dict?
                self._per_letter_annotations = _RestrictedDict(length=0)
            else:
                try:
                    self._per_letter_annotations = _RestrictedDict(length=len(seq))
                except TypeError:
                    raise TypeError(
                        "seq argument should be a Seq object or similar"
                    ) from None
        else:
            # This will be handled via the property set function, which will
            # turn this into a _RestrictedDict and thus ensure all the values
            # in the dict are the right length
            self.letter_annotations = letter_annotations

        # annotations about parts of the sequence
        if features is None:
            features = []
        elif not isinstance(features, list):
            raise TypeError(
                "features argument should be a list (of SeqFeature objects)"
            )
        self.features = features

    # TODO - Just make this a read only property?
    def _set_per_letter_annotations(self, value):
        if not isinstance(value, dict):
            raise TypeError(
                "The per-letter-annotations should be a (restricted) dictionary."
            )
        # Turn this into a restricted-dictionary (and check the entries)
        try:
            self._per_letter_annotations = _RestrictedDict(length=len(self.seq))
        except AttributeError:
            # e.g. seq is None
            self._per_letter_annotations = _RestrictedDict(length=0)
        self._per_letter_annotations.update(value)

    letter_annotations = property(
        fget=lambda self: self._per_letter_annotations,
        fset=_set_per_letter_annotations,
        doc="""Dictionary of per-letter-annotation for the sequence.

        For example, this can hold quality scores used in FASTQ or QUAL files.
        Consider this example using Bio.SeqIO to read in an example Solexa
        variant FASTQ file as a SeqRecord:

        >>> from Bio import SeqIO
        >>> record = SeqIO.read("Quality/solexa_faked.fastq", "fastq-solexa")
        >>> print("%s %s" % (record.id, record.seq))
        slxa_0001_1_0001_01 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN
        >>> print(list(record.letter_annotations))
        ['solexa_quality']
        >>> print(record.letter_annotations["solexa_quality"])
        [40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5]

        The letter_annotations get sliced automatically if you slice the
        parent SeqRecord, for example taking the last ten bases:

        >>> sub_record = record[-10:]
        >>> print("%s %s" % (sub_record.id, sub_record.seq))
        slxa_0001_1_0001_01 ACGTNNNNNN
        >>> print(sub_record.letter_annotations["solexa_quality"])
        [4, 3, 2, 1, 0, -1, -2, -3, -4, -5]

        Any python sequence (i.e. list, tuple or string) can be recorded in
        the SeqRecord's letter_annotations dictionary as long as the length
        matches that of the SeqRecord's sequence.  e.g.

        >>> len(sub_record.letter_annotations)
        1
        >>> sub_record.letter_annotations["dummy"] = "abcdefghij"
        >>> len(sub_record.letter_annotations)
        2

        You can delete entries from the letter_annotations dictionary as usual:

        >>> del sub_record.letter_annotations["solexa_quality"]
        >>> sub_record.letter_annotations
        {'dummy': 'abcdefghij'}

        You can completely clear the dictionary easily as follows:

        >>> sub_record.letter_annotations = {}
        >>> sub_record.letter_annotations
        {}

        Note that if replacing the record's sequence with a sequence of a
        different length you must first clear the letter_annotations dict.
        """,
    )

    def _set_seq(self, value):
        # TODO - Add a deprecation warning that the seq should be write only?
        if self._per_letter_annotations:
            if len(self) != len(value):
                # TODO - Make this a warning? Silently empty the dictionary?
                raise ValueError("You must empty the letter annotations first!")
            else:
                # Leave the existing per letter annotations unchanged:
                self._seq = value
        else:
            self._seq = value
            # Reset the (empty) letter annotations dict with new length:
            try:
                self._per_letter_annotations = _RestrictedDict(length=len(self.seq))
            except AttributeError:
                # e.g. seq is None
                self._per_letter_annotations = _RestrictedDict(length=0)

    seq = property(
        fget=lambda self: self._seq,
        fset=_set_seq,
        doc="The sequence itself, as a Seq or MutableSeq object.",
    )

    def __getitem__(self, index):
        if isinstance(index, int):
            # NOTE - The sequence level annotation like the id, name, etc
            # do not really apply to a single character.  However, should
            # we try and expose any per-letter-annotation here?  If so how?
            return self.seq[index]
        elif isinstance(index, slice):
            if self.seq is None:
                raise ValueError("If the sequence is None, we cannot slice it.")
            parent_length = len(self)
            try:
                from BioSQL.BioSeq import DBSeqRecord

                biosql_available = True
            except ImportError:
                biosql_available = False

            if biosql_available and isinstance(self, DBSeqRecord):
                answer = SeqRecord(
                    self.seq[index],
                    id=self.id,
                    name=self.name,
                    description=self.description,
                )
            else:
                answer = self.__class__(
                    self.seq[index],
                    id=self.id,
                    name=self.name,
                    description=self.description,
                )
            # TODO - The description may no longer apply.
            # It would be safer to change it to something
            # generic like "edited" or the default value.

            # Don't copy the annotation dict and dbxefs list,
            # they may not apply to a subsequence.
            # answer.annotations = dict(self.annotations.items())
            # answer.dbxrefs = self.dbxrefs[:]
            # TODO - Review this in light of adding SeqRecord objects?

            # TODO - Cope with strides by generating ambiguous locations?
            start, stop, step = index.indices(parent_length)
            if step == 1:
                # Select relevant features, add them with shifted locations
                # assert str(self.seq)[index] == str(self.seq)[start:stop]
                for f in self.features:
                    if f.ref or f.ref_db:
                        # TODO - Implement this (with lots of tests)?
                        import warnings

                        warnings.warn(
                            "When slicing SeqRecord objects, any "
                            "SeqFeature referencing other sequences (e.g. "
                            "from segmented GenBank records) are ignored."
                        )
                        continue
                    if (
                        start <= f.location.nofuzzy_start
                        and f.location.nofuzzy_end <= stop
                    ):
                        answer.features.append(f._shift(-start))

            # Slice all the values to match the sliced sequence
            # (this should also work with strides, even negative strides):
            for key, value in self.letter_annotations.items():
                answer._per_letter_annotations[key] = value[index]

            return answer
        raise ValueError("Invalid index")

    def __iter__(self):
        return iter(self.seq)

    def __contains__(self, char):
        return char in self.seq

    def __str__(self):
        lines = []
        if self.id:
            lines.append(f"ID: {self.id}")
        if self.name:
            lines.append(f"Name: {self.name}")
        if self.description:
            lines.append(f"Description: {self.description}")
        if self.dbxrefs:
            lines.append("Database cross-references: " + ", ".join(self.dbxrefs))
        lines.append(f"Number of features: {len(self.features)}")
        for a in self.annotations:
            lines.append(f"/{a}={str(self.annotations[a])}")
        if self.letter_annotations:
            lines.append(
                "Per letter annotation for: " + ", ".join(self.letter_annotations)
            )
        # Don't want to include the entire sequence
        lines.append(repr(self.seq))
        return "\n".join(lines)

    def __repr__(self):
        return (
            f"{self.__class__.__name__}(seq={self.seq!r}, id={self.id!r},"
            f" name={self.name!r}, description={self.description!r},"
            f" dbxrefs={self.dbxrefs!r})"
        )

    def format(self, format):
        # See also the __format__ method
        # See also the Bio.Align.Generic.Alignment class and its format()
        return self.__format__(format)

    def __format__(self, format_spec):
        if not format_spec:
            # Follow python convention and default to using __str__
            return str(self)
        from Bio import SeqIO

        # Easy case, can call string-building function directly
        if format_spec in SeqIO._FormatToString:
            return SeqIO._FormatToString[format_spec](self)

        # Harder case, make a temp handle instead
        handle = StringIO()
        try:
            SeqIO.write(self, handle, format_spec)
        except StreamModeError:
            raise ValueError(
                "Binary format %s cannot be used with SeqRecord format method"
                % format_spec
            ) from None
        return handle.getvalue()

    def __len__(self):
        return len(self.seq)

    def __lt__(self, other):
        raise NotImplementedError(_NO_SEQRECORD_COMPARISON)

    def __le___(self, other):
        raise NotImplementedError(_NO_SEQRECORD_COMPARISON)

    def __eq__(self, other):
        raise NotImplementedError(_NO_SEQRECORD_COMPARISON)

    def __ne__(self, other):
        raise NotImplementedError(_NO_SEQRECORD_COMPARISON)

    def __gt__(self, other):
        raise NotImplementedError(_NO_SEQRECORD_COMPARISON)

    def __ge__(self, other):
        raise NotImplementedError(_NO_SEQRECORD_COMPARISON)

    def __bool__(self):
        return True

    def __add__(self, other):
        if not isinstance(other, SeqRecord):
            # Assume it is a string or a Seq.
            # Note can't transfer any per-letter-annotations
            return SeqRecord(
                self.seq + other,
                id=self.id,
                name=self.name,
                description=self.description,
                features=self.features[:],
                annotations=self.annotations.copy(),
                dbxrefs=self.dbxrefs[:],
            )
        # Adding two SeqRecord objects... must merge annotation.
        answer = SeqRecord(
            self.seq + other.seq, features=self.features[:], dbxrefs=self.dbxrefs[:]
        )
        # Will take all the features and all the db cross refs,
        length = len(self)
        for f in other.features:
            answer.features.append(f._shift(length))
        del length
        for ref in other.dbxrefs:
            if ref not in answer.dbxrefs:
                answer.dbxrefs.append(ref)
        # Take common id/name/description/annotation
        if self.id == other.id:
            answer.id = self.id
        if self.name == other.name:
            answer.name = self.name
        if self.description == other.description:
            answer.description = self.description
        for k, v in self.annotations.items():
            if k in other.annotations and other.annotations[k] == v:
                answer.annotations[k] = v
        # Can append matching per-letter-annotation
        for k, v in self.letter_annotations.items():
            if k in other.letter_annotations:
                answer.letter_annotations[k] = v + other.letter_annotations[k]
        return answer

    def __radd__(self, other):
        if isinstance(other, SeqRecord):
            raise RuntimeError(
                "This should have happened via the __add__ of "
                "the other SeqRecord being added!"
            )
        # Assume it is a string or a Seq.
        # Note can't transfer any per-letter-annotations
        offset = len(other)
        return SeqRecord(
            other + self.seq,
            id=self.id,
            name=self.name,
            description=self.description,
            features=[f._shift(offset) for f in self.features],
            annotations=self.annotations.copy(),
            dbxrefs=self.dbxrefs[:],
        )

    def upper(self):
        return SeqRecord(
            self.seq.upper(),
            id=self.id,
            name=self.name,
            description=self.description,
            dbxrefs=self.dbxrefs[:],
            features=self.features[:],
            annotations=self.annotations.copy(),
            letter_annotations=self.letter_annotations.copy(),
        )

    def lower(self):
        return SeqRecord(
            self.seq.lower(),
            id=self.id,
            name=self.name,
            description=self.description,
            dbxrefs=self.dbxrefs[:],
            features=self.features[:],
            annotations=self.annotations.copy(),
            letter_annotations=self.letter_annotations.copy(),
        )

    def reverse_complement(
        self,
        id=False,
        name=False,
        description=False,
        features=True,
        annotations=False,
        letter_annotations=True,
        dbxrefs=False,
    ):
        from Bio.Seq import MutableSeq  # Lazy to avoid circular imports

        if "protein" in self.annotations.get("molecule_type", ""):
            raise ValueError("Proteins do not have complements!")
        if "RNA" in self.annotations.get("molecule_type", ""):
            if isinstance(self.seq, MutableSeq):
                # Does not currently have reverse_complement_rna method:
                answer = SeqRecord(self.seq.toseq().reverse_complement_rna())
            else:
                answer = SeqRecord(self.seq.reverse_complement_rna())
        else:
            # Default to DNA
            if isinstance(self.seq, MutableSeq):
                # Currently the MutableSeq reverse complement is in situ
                answer = SeqRecord(self.seq.toseq().reverse_complement())
            else:
                answer = SeqRecord(self.seq.reverse_complement())
        if isinstance(id, str):
            answer.id = id
        elif id:
            answer.id = self.id
        if isinstance(name, str):
            answer.name = name
        elif name:
            answer.name = self.name
        if isinstance(description, str):
            answer.description = description
        elif description:
            answer.description = self.description
        if isinstance(dbxrefs, list):
            answer.dbxrefs = dbxrefs
        elif dbxrefs:
            # Copy the old dbxrefs
            answer.dbxrefs = self.dbxrefs[:]
        if isinstance(features, list):
            answer.features = features
        elif features:
            # Copy the old features, adjusting location and string
            length = len(answer)
            answer.features = [f._flip(length) for f in self.features]
            # The old list should have been sorted by start location,
            # reversing it will leave it sorted by what is now the end position,
            # so we need to resort in case of overlapping features.
            # NOTE - In the common case of gene before CDS (and similar) with
            # the exact same locations, this will still maintain gene before CDS
            answer.features.sort(key=lambda x: x.location.start.position)
        if isinstance(annotations, dict):
            answer.annotations = annotations
        elif annotations:
            # Copy the old annotations,
            answer.annotations = self.annotations.copy()
        if isinstance(letter_annotations, dict):
            answer.letter_annotations = letter_annotations
        elif letter_annotations:
            # Copy the old per letter annotations, reversing them
            for key, value in self.letter_annotations.items():
                answer._per_letter_annotations[key] = value[::-1]
        return answer

    def translate(
        self,
        # Seq translation arguments:
        table="Standard",
        stop_symbol="*",
        to_stop=False,
        cds=False,
        gap=None,
        # SeqRecord annotation arguments:
        id=False,
        name=False,
        description=False,
        features=False,
        annotations=False,
        letter_annotations=False,
        dbxrefs=False,
    ):
        if "protein" == self.annotations.get("molecule_type", ""):
            raise ValueError("Proteins cannot be translated!")
        answer = SeqRecord(
            self.seq.translate(
                table=table, stop_symbol=stop_symbol, to_stop=to_stop, cds=cds, gap=gap
            )
        )
        if isinstance(id, str):
            answer.id = id
        elif id:
            answer.id = self.id
        if isinstance(name, str):
            answer.name = name
        elif name:
            answer.name = self.name
        if isinstance(description, str):
            answer.description = description
        elif description:
            answer.description = self.description
        if isinstance(dbxrefs, list):
            answer.dbxrefs = dbxrefs
        elif dbxrefs:
            # Copy the old dbxrefs
            answer.dbxrefs = self.dbxrefs[:]
        if isinstance(features, list):
            answer.features = features
        elif features:
            # Does not make sense to copy old features as locations wrong
            raise TypeError("Unexpected features argument %r" % features)
        if isinstance(annotations, dict):
            answer.annotations = annotations
        elif annotations:
            # Copy the old annotations
            answer.annotations = self.annotations.copy()
        # Set/update to protein:
        answer.annotations["molecule_type"] = "protein"
        if isinstance(letter_annotations, dict):
            answer.letter_annotations = letter_annotations
        elif letter_annotations:
            # Does not make sense to copy these as length now wrong
            raise TypeError(
                "Unexpected letter_annotations argument %r" % letter_annotations
            )
        return answer
