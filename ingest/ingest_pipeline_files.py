import csv
import mimetypes

from clusters import Cluster
from mtx import Mtx


class IngestFiles:
    def __init__(self, file_path, allowed_file_types):
        self.file_type, self.file = sel.open_file(file_path)
        self.allowed_file_types = allowed_file_types

    # Inherited function
    def open_file(self, file_path):
        # Check file type
        file_type = get_file_type(file_path)
        # See if file type is allowed
        if file_type in self.self.allowed_file_types:
            # open file
            file_connections = {
                'text/csv': open,
                'text/plain': self.open_csv(open(file_path)),
                'text/tab-separated-values': open,
            }

        # Return file object
        return file_type, file_connections.get(file_type)(file_path)

    def extract(self):

    def get_file_type(self, file_path):
        return mime.guess_type(file_path)

    def open_csv(self, opened_file_object):
        csv.register_dialect('myDialect',
                             delimiter=',',
                             quoting=csv.QUOTE_ALL,
                             skipinitialspace=True)
        return csv.reader(opened_file_object, dialect='myDialect')

    def extract_csv(self):
        for row in self.file:
            yield next(row)
