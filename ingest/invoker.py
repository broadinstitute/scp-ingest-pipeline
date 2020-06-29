from command import Command

class IngestInvoker:
    """
    The invoker has a reference to the Command
    and the command executes the appropriate action of the receiver.
    """

    def __init__(self, command: Command):
        self._command = command
        self._command_list = []  # type: List[Command]

    def add_command_to_list(self, command: Command):
        self._command_list.append(command)

    def execute_commands(self):
        """
        Execute all the saved commands, then empty the list.
        """
        for cmd in self._command_list:
            cmd.execute()

        self._command_list.clear()

    def invoke_ingest(self):
        """
        Invokes a command and executes
        the appropriate action (command.execute())
        """
        self._command.execute_ingest()
