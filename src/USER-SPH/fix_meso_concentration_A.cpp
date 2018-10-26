/* ----------------------------------------------------------------------
 LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
 http://lammps.sandia.gov, Sandia National Laboratories
 Steve Plimpton, sjplimp@sandia.gov

 Copyright (2003) Sandia Corporation.  Under the terms of Contract
 DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 certain rights in this software.  This software is distributed under
 the GNU General Public License.

 See the README file in the top-level LAMMPS directory.
 ------------------------------------------------------------------------- */

#include <stdio.h>
#include <string.h>
#include "fix_meso_concentration_A.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "integrate.h"
#include "respa.h"
#include "memory.h"
#include "error.h"
#include "pair.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMesoConcentrationA::FixMesoConcentrationA(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg) {  
  if ((atom->e_flag != 1) || (atom->rho_flag != 1))
    error->all(FLERR,
        "fix meso/concentration command requires atom_style with both energy and density, e.g. meso");

  if (narg != 3)
    error->all(FLERR,"Illegal number of arguments for fix meso/concentrationA command");

  time_integrate = 0;

  // find the concentration property
  int fcA;
  int icA = atom->find_custom("cA", fcA);
  if (icA < 0)
    error->all(FLERR,
        "Can't find property cA for fix meso/concentrationA");
  cA = atom->dvector[icA];

  // find the concentration property
  int fdcA;
  int idcA = atom->find_custom("dcA", fdcA);
  if (icA < 0)
    error->all(FLERR,
        "Can't find property dcA for fix meso/concentrationA");
  dcA = atom->dvector[idcA];

  // find the solid-liquid interaction
  int fdmA;
  int idmA = atom->find_custom("dmA", fdmA);
  if (idmA < 0)
    error->all(FLERR,
        "Can't find property dmA for fix meso/concentrationA");
  dmA = atom->dvector[idmA];

  // find the mass of A property
  int fmA;
  int imA = atom->find_custom("mA", fmA);
  if (imA < 0)
    error->all(FLERR,
        "Can't find property mA for fix meso/concentrationA");
  mA = atom->dvector[imA];

  // Set the communication size
  comm_forward = 2;
}

/* ---------------------------------------------------------------------- */

FixMesoConcentrationA::~FixMesoConcentrationA() {
}

/* ---------------------------------------------------------------------- */

int FixMesoConcentrationA::setmask() {
  int mask = 0;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMesoConcentrationA::init() {
  dtcA = update->dt;
}

/* ---------------------------------------------------------------------- */

void FixMesoConcentrationA::final_integrate() {
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int *type = atom->type;

  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  int nall = nlocal + atom->nghost;

  for (int i = 0; i < nall; i++) {
    if (mask[i] & groupbit) {
      cA[i] += dtcA*dcA[i];
      // Only update mass for solid particles
      if (type[i] == 2)
	mA[i] += dtcA*dmA[i];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixMesoConcentrationA::reset_dt() {
  dtcA = update->dt;
}
