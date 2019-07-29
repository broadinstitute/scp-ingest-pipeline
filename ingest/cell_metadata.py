#! /usr/bin/python

from collections import defaultdict
from itertools import islice


class CellMetadata:
    """A class used to represent cell metadata.

    Attributes
    ----------
    headers : list of strings
        list of metadata names in metadata file
    metadata_types : list of strings
        list of metadata type designations from metadata file
    annotation_type : list of strings
        list of valid metadata type designations
    errors : dict
        dict of collected validation errors and warnings

    """
    def __init__(self, file_path):
        """
        Parameters
        ----------
        filepath : str
            The name of the metadata file

        """
        self.file = open(file_path, 'r')
        self.headers = self.file.readline().rstrip('\n').split('\t')
        self.metadata_types = self.file.readline().rstrip('\n').split('\t')
        self.annotation_type = ['group', 'numeric']
        self.errors = defaultdict(list)
        self.ontology = defaultdict(lambda: defaultdict(set))
        self.type = defaultdict(list)

    def validate_header_keyword(self):
        """Check metadata header row starts with NAME (case-insensitive).

        :return: boolean   True if valid, False otherwise
        """
        valid = False
        if self.headers[0].casefold() == 'NAME'.casefold():
            valid = True
            if self.headers[0] != 'NAME':
                # ToDO - capture warning below in error report
                print(
                    'Warning: metadata file keyword "NAME" provided as {x}'.
                    format(x=self.headers[0])
                )
        else:
            ### line below and similar in next method have autoformat oddities
            self.errors['format'].append(
                'Error: Metadata file header row malformed, missing NAME'
            )
        return valid

    def validate_type_keyword(self):
        """Check metadata second row starts with TYPE (case-insensitive).

        :return: boolean   True if valid, False otherwise
        """
        valid = False
        if self.metadata_types[0].casefold() == 'TYPE'.casefold():
            valid = True
            if self.metadata_types[0] != 'TYPE':
                # ToDO - capture warning below in error report
                ### investigate f-string formatting here
                print(
                    'Warning: metadata file keyword "TYPE" provided as {x}'.
                    format(x=self.metadata_types[0])
                )
        else:
            ### check black autoformatting on this long line
            self.errors['format'].append(
                'Error:  Metadata file TYPE row malformed, missing TYPE'
            )
        return valid

    def validate_type_annotations(self):
        """Check metadata second row contains only 'group' or 'numeric'.

        :return: boolean   True if valid, False otherwise
        """
        valid = False
        annot_err = False
        annots = []
        for t in self.metadata_types[1:]:
            if t not in self.annotation_type:
                annots.append(t)
                annot_err = True
        if annot_err:
            self.errors['format'].append(
                (
                    'Error: TYPE declarations should be "group" or "numeric"; '
                    'Please correct: {annots}'.format(
                        annots=', '.join(map(str, annots))
                    )
                )
            )
        else:
            valid = True
        return valid

    def validate_against_header_count(self, list):
        """Metadata header and type counts should match.

        :return: boolean   True if valid, False otherwise
        """
        valid = False
        if not len(self.headers) == len(list):
            self.errors['format'].append(
                str(
                    'Error: {x} TYPE declarations for {y} column headers'.
                    format(x=len(self.headers), y=len(list))
                )
            )
        else:
            valid = True
        return valid

    def validate_format(self):
        """Check all metadata file format criteria for file validity
        """
        self.validate_header_keyword()
        self.validate_type_keyword()
        self.validate_type_annotations()
        self.validate_against_header_count(self.metadata_types)
        if self.errors['format']:
            valid = False
        else:
            valid = True
        return valid

    def extract_txt(self, size: int = 500):
        """interate through file in 500 line chunks
        """
        while True:
            next_lines = list(islice(self.file, size))
            if not next_lines:
                break
            return (next_lines)
