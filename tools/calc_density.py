# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross
"""
Routine to calculate a density from a set of PDB files

This module defines a these public functions:
calc_density
writeDX

Can be executed from the command line as a stand-alone program
"""

import logging

import numpy as np
import six
from scipy.stats import norm

from protomslib import simulationobjects

logger = logging.getLogger('protoms')


def _fill_gauss(coord, grid, edges, spacing, std):
    """
    Fill a grid using Gaussian smoothing

    Parameters
    ----------
    coord : Numpy array
      the Cartesian coordinates to put on the grid
    grid  : Numpy array
      the 3D grid. Will be modified
    edges : list of Numpy array
      the edges of the grid
    spacing : float
      the grid spacing
    std  : float
      the sigma of the Gaussian distribution
    """
    # Maximum coordinate
    maxxyz = np.minimum(coord + 3 * std,
                        np.array([edges[0][-1], edges[1][-1], edges[2][-1]]))

    # Iterater over 3 standard deviations
    x = max(coord[0] - 3 * std, edges[0][0])
    gx = norm.pdf(x, coord[0], std)
    while x <= maxxyz[0]:
        y = max(coord[1] - 3 * std, edges[1][0])
        gy = norm.pdf(y, coord[1], std)
        while y <= maxxyz[1]:
            z = max(coord[2] - 3 * std, edges[2][0])
            gz = norm.pdf(z, coord[2], std)
            while z <= maxxyz[2]:
                # Increase the grid with a Gaussian probability density
                v = _voxel(np.array([x, y, z]), edges)
                grid[v[0], v[1], v[2]] = grid[v[0], v[1], v[2]] + gx * gy * gz

                z = z + spacing
                gz = norm.pdf(z, coord[2], std)
            y = y + spacing
            gy = norm.pdf(y, coord[1], std)
        x = x + spacing
        gx = norm.pdf(x, coord[0], std)


def _fill_sphere(coord, grid, edges, spacing, radius):
    """
    Fill a grid using spherical smoothing

    Parameters
    ----------
    coord : Numpy array
      the Cartesian coordinates to put on the grid
    grid  : Numpy array
      the 3D grid. Will be modified
    edges : list of Numpy array
      the edges of the grid
    spacing : float
      the grid spacing
    radius  : float
      the radius of the smoothing
  """
    # Maximum coordinate
    maxxyz = np.minimum(coord + radius,
                        np.array([edges[0][-1], edges[1][-1], edges[2][-1]]))

    # Iterate over the sphere
    rad2 = radius**2
    x = max(coord[0] - radius, edges[0][0])
    while x <= maxxyz[0]:
        y = max(coord[1] - radius, edges[1][0])
        while y <= maxxyz[1]:
            z = max(coord[2] - radius, edges[2][0])
            while z <= maxxyz[2]:
                # Check if we are on the sphere
                r2 = (x - coord[0])**2 + (y - coord[1])**2 + (z - coord[2])**2
                if r2 <= rad2:
                    # Increase grid with one
                    v = _voxel(np.array([x, y, z]), edges)
                    grid[v[0], v[1], v[2]] = grid[v[0], v[1], v[2]] + 1
                z = z + spacing
            y = y + spacing
        x = x + spacing


def _init_grid(xyz, spacing, padding):
    """
    Initialize a grid based on a list of x,y,z coordinates

    Parameters
    ----------
    xyz  : Numpy array
      Cartesian coordinates that should be covered by the grid
    spacing : float
      the grid spacing
    padding : float
      the space to add to minimum extent of the coordinates

    Returns
    -------
    Numpy array
      the grid
    list of Numpy arrays
      the edges of the grid
  """

    origin = np.floor(xyz.min(axis=0)) - padding
    tr = np.ceil(xyz.max(axis=0)) + padding
    length = tr - origin
    shape = np.array([int(l / spacing + 0.5) + 1 for l in length], dtype=int)
    grid = np.zeros(shape)
    edges = [np.linspace(origin[i], tr[i], shape[i]) for i in range(3)]
    return grid, edges


def _voxel(coord, edges):
    """
    Wrapper for the numpy digitize function to return the grid coordinates
    """
    return np.array(
        [np.digitize(coord, edges[i])[i] for i in range(3)], dtype=int)


def writeDX(grid, origin, spacing, filename):
    """
    Write the grid to file in DX-format

    Parameters
    ----------
    grid : Numpy array
      the 3D grid
    origin : NumpyArray
      the bottom-left coordinate of the grid
    spacing  : float
      the grid spacing
    filename : string
      the name of the DX file
  """
    f = open(filename, 'w')
    f.write("object 1 class gridpositions counts %5d%5d%5d\n" %
            (grid.shape[0], grid.shape[1], grid.shape[2]))
    f.write("origin %9.4f%9.4f%9.4f\n" % (origin[0], origin[1], origin[2]))
    f.write("delta %10.7f 0.0 0.0\n" % spacing)
    f.write("delta 0.0 %10.7f 0.0\n" % spacing)
    f.write("delta 0.0 0.0 %10.7f\n" % spacing)
    f.write("object 2 class gridconnections counts %5d%5d%5d\n" %
            (grid.shape[0], grid.shape[1], grid.shape[2]))
    f.write(
        "object 3 class array type double rank 0 items  %10d data follows\n" %
        (grid.shape[0] * grid.shape[1] * grid.shape[2]))
    cnt = 0
    for x in range(grid.shape[0]):
        for y in range(grid.shape[1]):
            for z in range(grid.shape[2]):
                f.write("%19.10E" % grid[x, y, z])
                cnt = cnt + 1
                if cnt >= 3:
                    cnt = 0
                    f.write("\n")
    if cnt > 0:
        f.write("\n")
    f.write('attribute "dep" string "positions"\n')
    f.write('object "regular positions regular connections" class field\n')
    f.write('component "positions" value 1\n')
    f.write('component "connections" value 2\n')
    f.write('component "data" value 3\n')
    f.close()


def calc_density(pdbfiles,
                 molname,
                 atomname,
                 padding=2.0,
                 extent=1.0,
                 spacing=0.5,
                 norm=None,
                 smoothing="sphere"):
    """
    Calculate the density from a set of PDB files

    Parameters
    ----------
    pdbfiles : PDBSet
      the PDB files
    residue : string
      the residue to extract
    atom : string
      the name of the atom to make a density from
    padding : float, optional
      extra space to add to the minimum extent
    extent : float, optional
      the extent of the smoothing
    spacing : float, optional
      the spacing of the grid
    norm : int, optional
      number to normalize the density by
    smoothing : string, optional
      kind of smoothing, can be either sphere or gauss

    Returns
    -------
    NumpyArray of integers
      the number of atom extracted from each PDB structure
    NumpyArray
      the 3D density
    dictionary
      grid properties, keys = spacing, min, and max
  """
    residue = molname.lower()
    atom = atomname.lower()

    # Extract coordinates from PDB-files
    xyz = []
    nextract = 0.0
    nfound = []
    for pdb in pdbfiles.pdbs:
        found = 0
        for i, res in six.iteritems(pdb.residues):
            if res.name.lower() != molname:
                continue
            for atom in res.atoms:
                if atom.name.strip().lower() == atomname:
                    xyz.append(atom.coords)
                    found = found + 1
        for i, sol in six.iteritems(pdb.solvents):
            if sol.name.lower() != molname:
                continue
            for iatom in sol.atoms:
                if iatom.name.strip().lower() == atomname:
                    xyz.append(iatom.coords)
                    found = found + 1
        nfound.append(found)
    xyz = np.array(xyz)
    nfound = np.array(nfound)
    if nfound.sum() == 0:
        return nfound, None, None

    # Create a grid
    grid, edges = _init_grid(xyz, spacing, padding)

    # Put all the xyz coordinates on the grid using some smoothing
    for coord in xyz:
        if smoothing == "sphere":
            _fill_sphere(coord, grid, edges, spacing, extent)
        else:
            _fill_gauss(coord, grid, edges, spacing, extent)

    # Normalize grid
    if norm is None or norm <= 0:
        grid = grid / len(pdbfiles.pdbs)
    else:
        grid = grid / norm

    prop = {}
    prop["spacing"] = spacing
    prop["min"] = [e[0] for e in edges]
    prop["max"] = [e[-1] for e in edges]
    return np.array(nfound), grid, prop


def get_arg_parser():
    import argparse

    # Setup a parser of the command-line arguments
    parser = argparse.ArgumentParser(
        description="Program to discretize atoms on a 3D grid", )
    parser.add_argument('-f', '--files', nargs="+", help="the input PDB-files")
    parser.add_argument(
        '-o',
        '--out',
        help="the name of the output grid-file in DX-format, "
             "default='grid.dx'",
        default="grid.dx")
    parser.add_argument(
        '-r',
        '--residue',
        help="the name of the residue to extract, default='wa1'",
        default="wa1")
    parser.add_argument(
        '-a',
        '--atom',
        help="the name of the atom to extract, default='o00'",
        default="o00")
    parser.add_argument(
        '-p',
        '--padding',
        type=float,
        help="the amount to increase the minimum box in each direction, "
             "default=2 A",
        default=2.0)
    parser.add_argument(
        '-s',
        '--spacing',
        type=float,
        help="the grid resolution, default=0.5 A",
        default=0.5)
    parser.add_argument(
        '-e',
        '--extent',
        type=float,
        help="the size of the smoothing, i.e. the extent of an atom, "
             "default=1A",
        default=1.0)
    parser.add_argument(
        '-n',
        '--norm',
        type=float,
        help="number used to normalize the grid, if not specified the "
             "number of input files is used"
    )
    parser.add_argument(
        '-t',
        '--type',
        choices=["sphere", "gaussian"],
        help="the type  of coordinate smoothing, should be either"
             " 'sphere', 'gaussian'",
        default="sphere")
    parser.add_argument(
        '--skip',
        type=int,
        help="the number of blocks to skip to calculate the density. "
             "default is 0. Skip must be greater or equal to 0",
        default=0)
    parser.add_argument(
        '--max',
        type=int,
        help="the upper block to use. default is 99999 which should "
             "make sure you will use all the available blocks. "
             "max must be greater or equal to 0",
        default=99999)
    return parser


#
# If this is run from the command-line
#
if __name__ == "__main__":

    args = get_arg_parser().parse_args()

    if not args.files:
        print("No input files! Nothing to do, so exit.")
        quit()

    # Fix negative values of skip and max
    if args.max < 0:
        args.max = 99999
    if args.skip <= 0:
        args.skip = -1

    # Read in PDB files
    if len(args.files) == 1:
        pdbfiles = simulationobjects.PDBSet()
        pdbfiles.read(args.files[0], skip=args.skip, readmax=args.max)
    else:
        pdbfiles = simulationobjects.PDBSet()
        for filename in args.files[args.skip:args.max + 1]:
            pdb = simulationobjects.PDBFile(filename=filename)
            pdbfiles.pdbs.append(pdb)

    nextracted, grid, prop = calc_density(
        pdbfiles,
        args.residue.lower(),
        args.atom.lower(),
        padding=args.padding,
        extent=args.extent,
        spacing=args.spacing,
        norm=args.norm,
        smoothing=args.type)

    print("Extracted atoms in each PDB: %s" % ", ".join("%d" % i
                                                        for i in nextracted))
    print("Extracted %.3f on average" % nextracted.mean())

    if grid is not None:
        writeDX(grid, prop["min"], args.spacing, args.out)
