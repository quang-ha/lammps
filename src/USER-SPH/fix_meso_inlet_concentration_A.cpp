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
#include "fix_meso_inlet_concentration_A.h"
#include <math.h>
#include <stdlib.h>
#include "atom.h"
#include "atom_vec.h"
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
#include "region.h"
#include "domain.h"
#include "lattice.h"
#include "pair.h"
#include "modify.h"
#include "assert.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMesoInletConcentrationA::FixMesoInletConcentrationA(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg) {  
  if ((atom->e_flag != 1) || (atom->rho_flag != 1))
    error->all(FLERR,
        "fix meso/concentration command requires atom_style with both energy and density, e.g. meso");

  // This fix require 5 variables
  int nnarg = 4;

  if (narg < nnarg)
    error->all(FLERR,"Illegal number of arguments for fix meso/concentration command");

  time_integrate = 0;

  // find the concentration property
  int fcA;
  int icA = atom->find_custom("cA", fcA);
  if (icA < 0)
    error->all(FLERR,
	       "Can't find property cA for fix meso/concAdiffusion"); 
  cA = atom->dvector[icA];

  // find the local concentration property
  int fdcA;
  int idcA = atom->find_custom("dcA", fdcA);
  if (idcA < 0)
    error->all(FLERR,
	       "Can't find property dcA for fix meso/concAdiffusion");  
  dcA = atom->dvector[idcA];

  // Required args
  int m = 3;
  // Get the inlet concentration
  cAin = atof(arg[m++]);

  // Check the region
  iregion = -1;
  idregion = NULL;
  scaleflag = 1;
  
  // Read options from 
  options(narg-nnarg,&arg[nnarg]);
  // error checks on region and its extent being inside simulation box

  if (iregion == -1) error->all(FLERR,"Must specify a region in fix meso/inletconcentrationA");
  if (domain->regions[iregion]->bboxflag == 0)
    error->all(FLERR,"Fix phase_change_insert_random region does not support a bounding box");
  if (domain->regions[iregion]->dynamic_check())
    error->all(FLERR,"Fix phase_change_insert_random region cannot be dynamic");

  xlo = domain->regions[iregion]->extent_xlo;
  xhi = domain->regions[iregion]->extent_xhi;
  ylo = domain->regions[iregion]->extent_ylo;
  yhi = domain->regions[iregion]->extent_yhi;
  zlo = domain->regions[iregion]->extent_zlo;
  zhi = domain->regions[iregion]->extent_zhi;
  
  if (domain->triclinic == 0) {
    if (xlo < domain->boxlo[0] || xhi > domain->boxhi[0] ||
	ylo < domain->boxlo[1] || yhi > domain->boxhi[1] ||
	zlo < domain->boxlo[2] || zhi > domain->boxhi[2])
      error->all(FLERR,"Phase change region extends outside simulation box");
  } else {
    if (xlo < domain->boxlo_bound[0] || xhi > domain->boxhi_bound[0] ||
	ylo < domain->boxlo_bound[1] || yhi > domain->boxhi_bound[1] ||
	zlo < domain->boxlo_bound[2] || zhi > domain->boxhi_bound[2])
      error->all(FLERR,"Phase change region extends outside simulation box");
  }
  
  // // set comm size needed by this fix
  // comm_forward = 1;
}

/* ---------------------------------------------------------------------- */

int FixMesoInletConcentrationA::setmask() {
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMesoInletConcentrationA::init() {
  // set index and check validity of region
  iregion = domain->find_region(idregion);
  if (iregion == -1)
    error->all(FLERR,"Region ID for fix meso/inletconcentrationA does not exist");
}

/* ----------------------------------------------------------------------
 allow for both per-type and per-atom mass
 ------------------------------------------------------------------------- */

void FixMesoInletConcentrationA::pre_exchange() {
  double **x = atom->x;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int i;

  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if ((x[i][0] >= xlo) && (x[i][0] <= xhi) &&
	  (x[i][1] >= ylo) && (x[i][1] <= yhi) &&
	  (x[i][2] >= zlo) && (x[i][2] <= zhi))
	{
	  cA[i] = cAin;
	}
    }
  }

  
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of input line                                      
------------------------------------------------------------------- */

void FixMesoInletConcentrationA::options(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR,"Illegal fix indent command");

  int iarg = 0;
  // default number of attempts
  maxattempt = 10;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix meso/inletconcentrationA command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
	error->all(FLERR,"Region ID for fix meso/inletconcentrationA does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"attempt") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix meso/inletconcentrationA command");
      maxattempt = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix meso/inletconcentrationA command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0)
	error->all(FLERR,"Illegal fix meso/inletconcentrationA command: 'units lattice' "
		   "is not implemented");
      else error->all(FLERR,"Illegal fix meso/inletconcentrationA command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix meso/inletconcentrationA command");
  }
}
