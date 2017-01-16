from __future__ import print_function

import unittest
import subprocess
import os
import shutil
import itertools
import string
import filecmp
import errno


class TestDefinitionError(BaseException):
    pass


class BaseTest(unittest.TestCase):
    """
    Base class for ProtoMS functional(1) tests.

    The test workflow consists of 3 stages:
        1. Input files are copied from a reference directory
        2. A ProtoMS executable is called using a provided set of arguments
        3. Output files are compared against the output files provided in the reference directory

    A new ProtoMS functional test should inherit from this class and set
    the following attributes:
        executable : string
            ProtoMS executable to be run during the test - relative to PROTOMSHOME
        ref_dir : string
            The directory containing the test input and reference files - relative to PROTOMSHOME/tests
        input_files : list of string
            List of files to be copied from 'ref_dir' before running 'executable'
        args : list of string
            The arguments to be provided to 'executable'
        output_files : list of string
            Files to compare against reference at end of test - differences will fail test

    Additionally, the test may set the following attributes:
        output_directories : list of string
            List of directories in which to check output files - relative to working directory
            Used f.e. when running the Fortran executable which stores output in an output directory
            Multiple directories are listed f.e. in the case of replica exchange simulations
            Default value is ["."], i.e. only the test working directory
        mpi_processes : int
            Number of MPI processes to use when calling 'executable'
            Providing a value <= 0 disables MPI, default value is 0 (no MPI)

    1) https://en.wikipedia.org/wiki/Functional_testing
    """

    # Required class attributes used to define a new test
    executable = None
    """ProtoMS executable to be run during simulation phase - location relative to PROTOMSHOME"""
    ref_dir = None
    """Directory containing the test reference files"""
    input_files = None
    """Files to copy to the test's working directory at the beginning of the test"""
    args = None
    """List of arguments provided to the executable during the test"""
    output_files = None
    """Output files to check against reference at the end of the test"""

    # Optional class attributes used to define a new test
    output_directories = ["."]
    """List of names of the output directories - default is the test's working directory"""
    mpi_processes = 0
    """Number of MPI processes to use - default is no MPI"""

    @classmethod
    def setUpClass(cls):
        super(BaseTest, cls).setUpClass()

        # Create a directory for each test to hold test files
        # Allows you to check them if the test fails
        cls._test_dir = os.path.join("test_files", cls.__name__)
        try:
            os.makedirs(cls._test_dir)
        except OSError as e:
            if errno.EEXIST != e.errno:
                raise e
        cls._start_dir = os.getcwd()
        os.chdir(cls._test_dir)

        try:
            protoms_env = os.environ["PROTOMSHOME"]
        except KeyError:
            raise EnvironmentError("PROTOMSHOME environment variable is not set.")

        # These are required for every test
        if None in (cls.ref_dir, cls.input_files, cls.executable, cls.args, cls.output_files):
            raise TestDefinitionError("Missing required test attributes in definition of test '{0}'. Check framework.py".format(cls.__name__))

        cls.executable = os.path.join(protoms_env, cls.executable)

        cls._all_files = cls.input_files[:]
        cls._full_ref_dir = os.path.join(protoms_env, cls.ref_dir)

        output_files = [os.path.join(d, f) for f in cls.output_files for d in cls.output_directories]
        cls._all_files.extend(output_files)

        # Delete files left over from previous test runs
        cls._helper_clean_files()

    @classmethod
    def tearDownClass(cls):
        super(BaseTest, cls).tearDownClass()
        # cls.helper_clean_files()
        os.chdir(cls._start_dir)

    @classmethod
    def _helper_clean_files(cls):
        for filename in cls._all_files:
            try:
                os.remove(filename)
            except OSError:
                pass

        try:
            for dir in cls.output_directories:
                if dir != "." and os.path.isdir(dir) and not os.listdir(dir):
                    os.rmdir(dir)
        except TypeError:
            pass

    def _helper_subprocess_call(self, args):
        return_code = subprocess.call(args)
        self.assertEqual(0, return_code)
        if return_code == 0:
            print("ProtoMS call successful")

    def _helper_check_output(self, output_files):
        compare_tools = CompareTools(self._full_ref_dir, verbose=True)

        for filename in output_files:
            self.assertTrue(os.path.exists(filename),
                            "Expected output file {0} is missing".format(filename))
            self.assertTrue(os.path.exists(os.path.join(self._full_ref_dir, filename)),
                            "Reference output file {0} is missing".format(filename))

        for filename in output_files:
            file_match = compare_tools.compare(filename)
            if not file_match:
                self._all_files.remove(filename)
            self.assertTrue(file_match, "Content mismatch between output and reference for file {0}".format(filename))

    def test(self):
        for filename in self.input_files:
            try:
                shutil.copy(os.path.join(self._full_ref_dir, filename), ".")
            except IOError:
                raise IOError("The required reference input file {0} could not be copied".format(filename))

        print("\nTEST_RUN\n")
        args = [self.executable] + self.args
        if self.mpi_processes > 0:
            args = ["mpirun", "-np", str(self.mpi_processes)] + args
        self._helper_subprocess_call(args)

        print("\nTEST_OUTPUT\n")
        out_files = [os.path.join(d, f) for f in self.output_files for d in self.output_directories]
        self._helper_check_output(out_files)


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

                if l1 == l2:
                    continue

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
                if l1 == l2:
                    continue

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
