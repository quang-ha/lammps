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

  // set up reneighbouring
  force_reneighbor = 1;
  next_reneighbor = update->ntimestep+1;
  nfirst = next_reneighbor;
}
/* ---------------------------------------------------------------------- */
FixPhaseChangeThreshold::~FixPhaseChangeThreshold()
{
  delete [] idregion;
}
/* ---------------------------------------------------------------------- */
int FixPhaseChangeThreshold::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}
/* ---------------------------------------------------------------------- */
void FixPhaseChangeThreshold::init()
{
  // set index and check validity of region
  iregion = domain->find_region(id_region);
  if (iregion == -1)
    error->all(FLERR,"Region ID for fix phase_change_threshold does not exist.");

  // Need a full neighbour list, built whenever re-neighbouring occurs
  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}
/* ---------------------------------------------------------------------- */
void FixPhaseChangeThreshold::init_list(int, NeighList *ptr)
{
  list = ptr;
}
/* ---------------------------------------------------------------------- */
/* perform phase change */
/* ---------------------------------------------------------------------- */
void FixPhaseChangeThreshold::pre_exchange()
{
  // just return if should not be called on this time step
  if (nexT_reneighbor != update->ntimestep) return;
 
  // TODO: Check if atom is in sub box below or above it

  int nins = 0;
  int nlocal = atom->nlocal;
  int* numneigh = list->numneigh;
  double **x = atom->x;
  double **v = atom->v;
  double **vest = atom->vest;
  double *rmas = atom->rmas;
  double *rho = atom->rho;
  double *cv = atom->cv;
  double *e = atom->e;
  dmass = atom->drho;
  int *type = atom->type;

  int nall;
  if (force->newtom) nall = atom->nlocal + atom->nghost;
  else nall = atom->nlocal;

  for (int i=0; i<nall; i++)
    {
      dmass[i] = 0.0;
    }

  for (int i=0; i<nlocal; i++)
    {
      double Ti = sph_energy2t(e[i], cv[i]);
      if ( (Ti < Tc) && type[i] != to_type)
	{
	  atom->type[i] = to_type;
	}
    }
  
}
/* ---------------------------------------------------------------------- */
/* pack entire state of Fix into one write */
/* ---------------------------------------------------------------------- */
void FixPhaseChangeThreshold::write_restart(FILE *fp)
{
  int n = 0;
  double list[4];
  list[n++] = next_reneighbor;

  if (comm->me == 0)
    {
      int size = n*sizeof(double);
      fwrite(&size, sizeof(int), 1, fp);
      fwrite(list, sizeof(double), n, fp);
    }
}
/* ---------------------------------------------------------------------- */
/* use state info from restart file to restart the fix /*
/* ---------------------------------------------------------------------- */
void FixPhaseChangeThreshold::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;

  seed = static_case<int> (list[n++]);
  nfirst = static_cast<int> (list[n++]);
  next_reneighbor = static_cast<int> (list[n++]);

  random->reset(seed);
}
