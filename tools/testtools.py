
import filecmp
import os
import itertools

import sys


class CompareTools:
    """
    Compare files to reference data
    """

    def __init__(self, refdir, verbose=False):
        """
        Create file comparator

        Parameters
        ----------
        refdir : directory containing reference files
        verbose : whether to print details to stdout
        """
        self.refdir = refdir
        self.verbose = verbose

        self.comparetools = {
            "accept"      : self.diff_filecmp,
            "warning"     : self.diff_filecmp,
            "cmd"         : self.diff_filecmp,  # Use ign_starts_with if paths vary on different systems
            "info"        : self.diff_text_ign_starts_with,
            ".pdb"        : self.diff_text_try_number,
            "restart"     : self.diff_text_try_number,
            "restart.prev": self.diff_text_try_number,
            "results"     : self.diff_text_try_number,
            "results_inst": self.diff_text_try_number,
            "UNKNOWN"     : self.diff_filecmp
        }

        self.ign_starts_with = {
            "cmd" : {"parfile"},
            "info": {"#",
                     "protoms3 started at",
                     "Reading parameter file",
                     "Starting simulation at",
                     "protoms3 completed at",
                     "These moves took"}
        }

    def compare(self, filetuple):
        """
        Check whether file contents match using appropriate method

        Parameters
        ----------
        filetuple : tuple containing relative filepath and file type from self.comparetools dictionary

        Returns
        -------
        boolean
            Whether files match
        """
        file, type = filetuple
        reffile = os.path.join(self.refdir, file)
        if type not in self.comparetools:
            type = os.path.splitext(type)[1]
            if type not in self.comparetools:
                print("Unrecognised type, using filecmp: {0}".format(type))
                type = "UNKNOWN"

        if self.comparetools[type](file, reffile, type):
            if self.verbose:
                print("File matched reference: {0}".format(file))
            return True
        else:
            print("File did not match reference: {0}".format(file))
            return False

    def diff_filecmp(self, file1, file2, type):
        """
        Check whether file contents match using filecmp.cmp()

        Parameters
        ----------
        file1 : path to first file to compare
        file2 : path to second file to compare
        type : file type being compared

        Returns
        -------
        boolean
            Whether files match
        """
        return filecmp.cmp(file1, file2)

    def diff_text(self, file1, file2, type, ign_white=False, ign_time_path=False, ign_starts_with=None):
        """
        Check whether text file contents match line by line

        Parameters
        ----------
        file1 : path to first file to compare
        file2 : path to second file to compare
        type : file type being compared

        Returns
        -------
        boolean
            Whether files match
        """
        with open(file1) as f1, open(file2) as f2:
            for l1, l2 in itertools.izip(f1, f2):
                skip = False

                if ign_starts_with is not None:
                    for start in ign_starts_with:
                        if l1.startswith(start):
                            skip = True
                            break
                    if skip:
                        continue

                if ign_white:
                    l1 = l1.translate(None, string.whitespace)
                    l2 = l2.translate(None, string.whitespace)

                if not l1 == l2:
                    if self.verbose:
                        print(l1)
                        print(l2)
                    return False

        return True

    def diff_text_ign_white(self, file1, file2, type):
        """
        Check whether text file contents match line by line after stripping all whitespace

        Parameters
        ----------
        file1 : path to first file to compare
        file2 : path to second file to compare
        type : file type being compared

        Returns
        -------
        boolean
            Whether files match
        """
        return self.diff_text(file1, file2, type, ign_white=True)

    def diff_text_ign_starts_with(self, file1, file2, type):
        """
        Check whether text file contents match ignoring lines starting with a set of strings

        Parameters
        ----------
        file1 : path to first file to compare
        file2 : path to second file to compare
        type : file type being compared

        Returns
        -------
        boolean
            Whether files match
        """
        return self.diff_text(file1, file2, type, ign_starts_with=self.ign_starts_with[type])

    def diff_text_try_number(self, file1, file2, type):
        """
        Check whether text file contents match after trying to convert tokens to float

        Parameters
        ----------
        file1 : path to first file to compare
        file2 : path to second file to compare
        type : file type being compared

        Returns
        -------
        boolean
            Whether files match
        """
        def convert_line(line):
            for tok in line.split():
                try:
                    yield float(tok)
                except ValueError:
                    yield tok

        with open(file1) as f1, open(file2) as f2:
            for l1, l2 in itertools.izip(f1, f2):
                tok1 = convert_line(l1)
                tok2 = convert_line(l2)
                sentinel = object()
                if not all(a == b for a, b in itertools.izip_longest(tok1, tok2, fillvalue=sentinel)):
                    if self.verbose:
                        print(l1)
                        print(l2)
                    return False
        return True

    def diff_null(self, file1, file2, type):
        """
        Null file checker, always return True

        Parameters
        ----------
        file1 : path to first file to compare
        file2 : path to second file to compare
        type : file type being compared

        Returns
        -------
        boolean
            Always return True
        """
        return True


if __name__ == "__main__":
    tool = CompareTools(sys.argv[1], verbose=True)
    filetype = os.path.basename(sys.argv[2])
    filetuple = (sys.argv[2], filetype)
    tool.compare(filetuple)
