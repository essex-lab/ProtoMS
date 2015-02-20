
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string>
#include "SimTKOpenMMCommon.h"
#include "RealVec.h"
#include "PME.h"
using namespace OpenMM;

extern "C" {

static RealOpenMM** allocateRealArray(int length, int width) {
  RealOpenMM** array = new RealOpenMM*[length];
  for (int i = 0; i < length; ++i)
    array[i] = new RealOpenMM[width];
  return array;
}

static void disposeRealArray(RealOpenMM** array, int size) {
  if (array) {
    for (int i = 0; i < size; ++i)
      delete[] array[i];
    delete[] array;
  }
}

/*
  Routine to calculate the energy of the mesh in a SPME calculation
  fnatoms - number of atoms
  fbox    - box sizes in nm
  fx      - X-coordinates in nm
  fy      - Y-coordinates in nm
  fz      - Z-coordinates in nm
  fcharge - charges
  fcutoff - non-bonded cut-off in nm
  ftol    - tolerance
  verbose - if larger than 0 the routine will print detailed information to standard output
  
  Returns the energy in kJ/mol

*/
double mesh_energy_(int *fnatoms,double *fbox,double *fx, double *fy, double *fz, double *fcharge, double *fcutoff,double *ftol,double *falpha, int *verbose)
{
  RealOpenMM cutoff = *fcutoff;
  RealOpenMM tolerance = *ftol;
  RealOpenMM alpha  = *falpha;
  int natoms = *fnatoms;

  RealOpenMM box[3];
  for (int i=0;i<3;i++)
    box[i] = fbox[i];

  // Read coordinates and charges
  RealOpenMM** charges = allocateRealArray(natoms,3);
  std::vector<RealVec> coordinates;
  std::vector<RealVec> forces;
  for (int i=0;i<natoms;i++) {
    charges[i][2] = fcharge[i];
    coordinates.push_back(RealVec(fx[i],fy[i],fz[i]));
    forces.push_back(RealVec(0.0,0.0,0.0));
  }  

  if (*verbose > 0) {
    std::cout << "Number of atoms: " << natoms << "\n";
    std::cout << "Box size (nm): " << box[0] << ", " << box[1] << ", " << box[2] << "\n";
  }

  // Calculate Ewald parameter alpha/beta
  //alpha = (1.0/cutoff)*std::sqrt(-log(2.0*tolerance));
  if (*verbose > 0) std::cout << "Ewald width: " << 1.0/alpha << "\n";
  if (*verbose > 0) std::cout << "Ewald tolerance: " << tolerance << "\n";

  // Calculate mesh size
  int mesh[3];
  mesh[0] = (int) ceil(2*alpha*box[0]/(3*pow(tolerance, 0.2)));
  mesh[1] = (int) ceil(2*alpha*box[1]/(3*pow(tolerance, 0.2)));
  mesh[2] = (int) ceil(2*alpha*box[2]/(3*pow(tolerance, 0.2)));
  if (*verbose > 0)  std::cout << "Mesh: " << mesh[0] << ", " << mesh[1] << ", " << mesh[2] << "\n";
 
  // Initialize and perform the reciprocal part of PME
  pme_t pmedata;
  RealOpenMM energy = 0.0;
  RealOpenMM virial[3][3];
  pme_init(&pmedata,alpha,natoms,mesh,4,1);
  pme_exec(pmedata,coordinates,forces,charges,box,&energy,virial);
  pme_destroy(pmedata);

  if (*verbose > 0)  std::cout << "Reciprocal energy: " << energy << "\n";

  disposeRealArray(charges,natoms);
  return energy;
}

}
