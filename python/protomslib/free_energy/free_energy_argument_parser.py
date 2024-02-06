import argparse


class FEArgumentParser(argparse.ArgumentParser):
    """Thin wrapper around argparse.ArgumentParser designed to allow
    specification of individual arguments that cannot be used at the
    same time as one another."""

    def __init__(self, *args, **kwargs):
        # use a dict to store clashes for arguments
        self.clashes = {}
        # inherit any clashes belonging to parent argument parsers
        for parent in kwargs.get("parents", {}):
            self.clashes.update(parent.clashes)

        argparse.ArgumentParser.__init__(self, *args, **kwargs)

    def parse_args(self, *args, **kwargs):
        """Thin wrapper around ArgumentParser.parse_args that checks
        for clashing arguments."""
        parsed = argparse.ArgumentParser.parse_args(self, *args, **kwargs)
        self.check_clashes(parsed)
        return parsed

    def add_argument(self, *args, **kwargs):
        """Thin wrapper around ArgumentParser.add_argument that
        uses additional keyword argument 'clashes'. Clashes should
        be a list that contains the names of other arguments which
        cannot be provided at the same time as this one."""
        # figure out what the argument will be called in the parser namespace
        try:
            name = kwargs["dest"]
        except KeyError:
            name = args[-1].strip("-")

        # store clashes if these are provided, otherwise empty list
        try:
            self.clashes[name] = kwargs.pop("clashes")
        except KeyError:
            self.clashes[name] = []

        argparse.ArgumentParser.add_argument(self, *args, **kwargs)

    def check_clashes(self, parsed):
        """Check that arguments in Namespace parsed are compatible,
        according to provided argument clashes."""
        for arg in self.clashes:
            # if parsed.arg has its default value it WAS NOT used so ignore
            if getattr(parsed, arg, None) == self.get_default(arg):
                continue
            for clash in self.clashes[arg]:
                # check clashes for this argument
                # if clashing argument has non-default value it WAS used
                # so complain and quit
                if getattr(parsed, clash) != self.get_default(clash):
                    self.error(
                        "Cannot provide both --%s and --%s arguments"
                        % (arg, clash)
                    )
