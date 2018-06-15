import logging
import numpy as np
import os
import subprocess

logger = logging.getLogger('protoms')


def _is_float(num):
    """
    Check whether a string is convertible to a float

    Parameters
    ----------
    num : string
      the string which might be convertible to float

    Returns
    -------
    boolean
      whether the string is convertible to float
    """
    try:
        float(num)
    except (ValueError, TypeError):
        return False
    return True


def _get_prefix(filename):
    """
    Remove extension (including period from a filename)

    Parameters
    ----------
    filename : string
      the filename to modify

    Returns
    -------
    string
      the filename without extension
    """
    h, t = os.path.splitext(filename)
    return h


def _locate_file(filename, folders):
    """
    Find a file

    Tries to find the file as it is given or in
    any of the folders given

    Parameters
    ----------
    filename : string
      the name of the file to find
    folders : list of strings
      folders to search for the file

    Returns
    -------
    string or None
      the full filename or None if the file could not be found
    """
    # Try to see if the filename as given exists
    if os.path.isfile(filename):
        return filename
    else:
        # Remove everything but the actual filename
        h, t = os.path.split(filename)
        # Loop over all folder and try to find it there
        for f in folders:
            test = os.path.join(f, t)
            if os.path.isfile(test):
                return test
    # If we haven't found it up to now, give up and return None
    return None


def rotmat_x(alpha):
    """
    3D rotation matrix for rotations about x-axis.

    Parameters
    ----------
    alpha : float or integer
      the angle (in radians) by which rotation is performed about the x-axis.

    Returns
    -------
    numpy array
      the rotation matrix for the desired angle.
    """
    return (np.mat([[1.0, 0.0, 0.0], [0.0, np.cos(alpha), -np.sin(alpha)],
                    [0.0, np.sin(alpha), np.cos(alpha)]]))


def rotmat_y(beta):
    """
    3D rotation matrix for rotations about x-axis.

    Parameters
    ----------
    beta : float or integer
      the angle (in radians) by which rotation is performed about the y-axis.

    Returns
    -------
    numpy array
      the rotation matrix for the desired angle.
    """
    return (np.mat([[np.cos(beta), 0.0, np.sin(beta)], [0.0, 1.0, 0.0],
                    [-np.sin(beta), 0.0, np.cos(beta)]]))


def rotmat_z(gamma):
    """
    3D rotation matrix for rotations about x-axis.

    Parameters
    ----------
    gamma : float or integer
      the angle (in radians) by which rotation is performed about the z-axis.

    Returns
    -------
    numpy array
      the rotation matrix for the desired angle.
    """
    return (np.mat([[np.cos(gamma), -np.sin(gamma), 0.0],
                    [np.sin(gamma), np.cos(gamma), 0.0], [0.0, 0.0, 1.0]]))


def _cleanup(tarlist):
    """
    Clean up extra files

    Parameters
    ----------
    tarlist : list of string
      the files to be cleaned up
    """
    tarlist2 = []
    for filename in tarlist:
        if filename in tarlist2:
            continue
        if filename.find(os.environ["PROTOMSHOME"]) == 0:
            continue
        tarlist2.append(filename)

    logger.info("")
    logger.info("Cleaning up and saving extra files to prep_files.tar")
    logger.debug("The files are: %s" % " ".join(tarlist2))
    subprocess.call(
        "tar -cf prep_files.tar %s" % " ".join(tarlist2), shell=True)
    subprocess.call("rm -f %s" % " ".join(tarlist2), shell=True)
