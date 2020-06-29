import abc
import sys

from command import Command
from expression_files.dense import Dense

sys.path.append("../ingest")


class DenseCommand(Command):
    """
    A concrete / specific Command class that calls a
    specific or an appropriate action of a method from
    Dense, a receiver class.
    Args:
        dense (Dense): Receiver class to be attached to the command
    """

    def __init__(self, dense: Dense):
        self._dense = dense

    def execute_ingest(self):
        """
        Delegate ingest methods of the Dense receiver.
        """
        for gene_docs, data_array_documents in self._dense.transform():
            load_status = self._dense.load_expression_file(
                gene_docs, data_array_documents)
        self._dense.close()
        return 0
