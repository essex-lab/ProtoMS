from __future__ import print_function

import unittest
import subprocess
import os
import shutil
import itertools
import string
import filecmp
# import tempfile


class BaseTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        super(BaseTest, cls).setUpClass()

        try:
            protoms_env = os.environ["PROTOMSHOME"]
        except KeyError:
            raise Exception("PROTOMSHOME environment variable is not set.")

        cls.setup_executable = os.path.join(protoms_env, "protoms.py")
        cls.simulation_executable = os.path.join(protoms_env, "build", "protoms3")

        # Temporary directory - files will be deleted after test
        # cls.test_dir = tempfile.mkdtemp(prefix=".")
        # os.chdir(cls.test_dir)

        cls.all_files = []

        try:
            cls.ref_dir
            cls.copy_files

            cls.all_files.extend(cls.copy_files)
        except AttributeError:
            raise AttributeError("You must provide all required data when defining a test case. Check framework.py")

        try:
            cls.setup_args
            cls.setup_output_files
            cls.do_test_setup = True

            cls.all_files.extend(cls.setup_output_files)
        except AttributeError:
            cls.do_test_setup = False

        try:
            cls.simulation_args
            cls.simulation_output_directory
            cls.simulation_output_files
            cls.do_test_simulation = True

            output_files = [os.path.join(cls.simulation_output_directory, f) for f in cls.simulation_output_files]
            cls.all_files.extend(output_files)
        except AttributeError:
            cls.do_test_simulation = False

        if not cls.do_test_setup and not cls.do_test_simulation:
            raise Exception("Test case runs neither setup nor simulation, it is defined incorrectly.")

        cls.helper_clean_files()

        cls.full_ref_dir = os.path.join(protoms_env, cls.ref_dir)
        print(cls.full_ref_dir)

    @classmethod
    def tearDownClass(cls):
        super(BaseTest, cls).tearDownClass()
        cls.helper_clean_files()

    @classmethod
    def helper_clean_files(cls):
        for filename in cls.all_files:
            try:
                os.remove(filename)
            except OSError:
                pass
        try:
            if os.path.isdir(cls.simulation_output_directory) and not os.listdir(cls.simulation_output_directory):
                os.rmdir(cls.simulation_output_directory)
        except AttributeError:
            pass

    def helper_subprocess_call(self, args):
        return_code = subprocess.call(args)
        self.assertEqual(0, return_code)
        if return_code == 0:
            print("ProtoMS call successful")

    def helper_check_output(self, stage, output_files):
        compare_tools = CompareTools(self.full_ref_dir, verbose=True)

        for filename in output_files:
            self.assertTrue(os.path.exists(filename),
                            "Expected {0} output file {1} is missing".format(stage, filename))

        for filename in output_files:
            file_match = compare_tools.compare(filename)
            if not file_match:
                print(output_files)
                self.all_files.remove(filename)
            self.assertTrue(file_match, "Content mismatch between {0} output and reference for file {1}".format(stage, filename))

    def test(self):
        # os.chdir(self.test_dir)

        for filename in self.copy_files:
            shutil.copy(os.path.join(self.full_ref_dir, filename), ".")

        if self.do_test_setup:
            print("\nTEST_SETUP_RUN\n")
            args = [self.setup_executable] + self.setup_args
            self.helper_subprocess_call(args)

            print("\nTEST_SETUP_OUTPUT\n")
            self.helper_check_output("setup", self.setup_output_files)

        if self.do_test_simulation:
            print("\nTEST_SIMULATION_RUN\n")
            args = [self.simulation_executable] + self.simulation_args
            self.helper_subprocess_call(args)

            print("\nTEST_SIMULATION_OUTPUT\n")
            sim_out_files = [os.path.join(self.simulation_output_directory, filename) for filename in self.simulation_output_files]
            self.helper_check_output("simulation", sim_out_files)


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

        # Which comparison tool to use for each file format
        self.comparetools = {
            "accept"      : self.diff_text_try_number,
            "warning"     : self.diff_text_try_number,
            ".cmd"        : self.diff_text,  # Use ign_starts_with if paths vary on different systems
            "info"        : self.diff_text_ign_starts_with,
            ".pdb"        : self.diff_text_try_number,
            "restart"     : self.diff_text_try_number,
            "restart.prev": self.diff_text_try_number,
            "results"     : self.diff_text_try_number,
            "results_inst": self.diff_text_try_number,
            None          : self.diff_filecmp
        }

        # Lines beginning with these strings will be ignored when using diff_text_ign_starts_with
        self.ign_starts_with = {
            "cmd" : {"parfile"},
            "info": {"#",
                     "protoms3 started at",
                     "Reading parameter file",
                     "Starting simulation at",
                     "protoms3 completed at",
                     "These moves took"}
        }

    def determine_file_type(self, filename):
        """
        Determine the type of the tested file from the dictionary self.comparetools

        Parameters
        ----------
        filename : name of file for which to determine type

        Returns
        -------
        str
            Type of file, from dictionary self.comparetools
        """
        filetype = os.path.basename(filename)

        if filetype not in self.comparetools:
            filetype = os.path.splitext(filename)[1]
            if not filetype:
                filetype = os.path.basename(filename)

        if filetype not in self.comparetools:
            if self.verbose:
                print("Unrecognised type for file {0}, using default filecmp.".format(filename))
            filetype = None

        return filetype

    def compare(self, filename):
        """
        Check whether file contents match using appropriate method

        Parameters
        ----------
        filename : name of file to compare against reference

        Returns
        -------
        boolean
            Whether files match
        """
        reffile = os.path.join(self.refdir, filename)

        filetype = self.determine_file_type(filename)

        if self.comparetools[filetype](filename, reffile, filetype):
            if self.verbose:
                print("File matched reference: {0}".format(filename))
            return True
        else:
            print("File did not match reference: {0}".format(filename))
            return False

    def diff_filecmp(self, file1, file2, filetype):
        """
        Check whether file contents match using filecmp.cmp()

        Parameters
        ----------
        file1 : path to first file to compare
        file2 : path to second file to compare
        filetype : file type being compared

        Returns
        -------
        boolean
            Whether files match
        """
        return filecmp.cmp(file1, file2)

    def diff_text(self, file1, file2, filetype, ign_white=False, ign_time_path=False, ign_starts_with=None):
        """
        Check whether text file contents match line by line

        Parameters
        ----------
        file1 : path to first file to compare
        file2 : path to second file to compare
        filetype : file type being compared

        Returns
        -------
        boolean
            Whether files match
        """
        if ign_time_path:
            raise NotImplementedError("Text diff with time and path exceptions is not yet implemented.")

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

    def diff_text_ign_white(self, file1, file2, filetype):
        """
        Check whether text file contents match line by line after stripping all whitespace

        Parameters
        ----------
        file1 : path to first file to compare
        file2 : path to second file to compare
        filetype : file type being compared

        Returns
        -------
        boolean
            Whether files match
        """
        return self.diff_text(file1, file2, filetype, ign_white=True)

    def diff_text_ign_starts_with(self, file1, file2, filetype):
        """
        Check whether text file contents match ignoring lines starting with a set of strings

        Parameters
        ----------
        file1 : path to first file to compare
        file2 : path to second file to compare
        filetype : file type being compared

        Returns
        -------
        boolean
            Whether files match
        """
        return self.diff_text(file1, file2, filetype, ign_starts_with=self.ign_starts_with[filetype])

    def diff_text_try_number(self, file1, file2, filetype):
        """
        Check whether text file contents match after trying to convert tokens to float

        Parameters
        ----------
        file1 : path to first file to compare
        file2 : path to second file to compare
        filetype : file type being compared

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

    def diff_null(self, file1, file2, filetype):
        """
        Null file checker, always return True

        Parameters
        ----------
        file1 : path to first file to compare
        file2 : path to second file to compare
        filetype : file type being compared

        Returns
        -------
        boolean
            Always return True
        """
        return True
