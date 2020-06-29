from command import Command

class IngestInvoker:
    """
    Has a reference to the Command, and can execute the method.
    Notice how the command.execute() is never directly called,
    but always through the invoker.
    The action invoked is decoupled from the action performed
    by the Receiver.
    The Invoker (self) invokes a Command (LunchCommand),
    and the Command executes the appropriate action (command.execute())
    """

    def __init__(self, command: Command):
        self._command = command
        self._command_list = []  # type: List[Command]

    def set_command(self, command: Command):
        self.command = command

    def get_command(self):
        print(self.command.__class__.__name__)

    def add_command_to_list(self, command: Command):
        self._command_list.append(command)

    def execute_commands(self):
        """
        Execute all the saved commands, then empty the list.
        """
        for cmd in self._command_list:
            cmd.execute()

        self._command_list.clear()

    def execute(self):
        self._command.execute()
