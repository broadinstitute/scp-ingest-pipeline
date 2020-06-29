import abc


class Command(metaclass=abc.ABCMeta):
    """
    The command interface that declares a method (execute) for a particular
    action.
    """
    @abc.abstractmethod
    def execute(self):
        pass
