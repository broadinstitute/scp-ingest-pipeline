import abc


class Ingestor(metaclass=abc.ABCMeta):
    """
    The command interface that declares a method (execute) for a particular
    action.
    """
    @abc.abstractmethod
    def execute_ingest(self, ):
        pass

    """
    Returns true/false for whether this ingestor can handle the given file type
    type should correspond to the
    """
    @abc.abstractmethod
    def matches_file_type(self, type):
        pass
