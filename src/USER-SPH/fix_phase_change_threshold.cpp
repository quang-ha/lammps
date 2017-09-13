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
#include "assert.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_phase_change_threshold.h"
#include "sph_kernel_quintic.h"
#include "sph_energy_equation.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "pair.h"
#include "random_park.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "comm.h"
#include "domain.h"
#include "lattice.h"
#include "region.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */
FixPhaseChangeThreshold::FixPhaseChangeThreshold(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  // communicate energy change due to phase change
  int nnarg = 5;
  if (narg < nnarg) error->all(FLERR, "Illegal fix phase_change_threshold command");

  restart_global = 1; // not sure what this do
  time_depend = 1;    // not sure what this do

  // required args
  int m = 3; // arguments start from 3
  Tc = atof(arg[m++]);
  from_type = atoi(arg[m++]);
  to_type = atoi(arg[m++]);
  cutoff = atof(arg[m++]);

  // check region id
  iregion = -1;
  idregion = NULL;
  scaleflag = 1;

  // read options from end of input line
  options(narg-nnarg, &arg[nnarg]);

  // error checks on region and its extent being inside simulation box
  if (iregion == -1) error->all(FLERR,"Must specify a region in fix phase_change_insert_random");
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
  
}
/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */
