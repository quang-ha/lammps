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
#include "fix_meso_precipitation_A.h"
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

/* TODO: Make this fix run in parallel */

/* ---------------------------------------------------------------------- */

FixMesoPrecipitationA::FixMesoPrecipitationA(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg) {

  if ((atom->e_flag != 1) || (atom->rho_flag != 1))
    error->all(FLERR,
        "fix meso/precipitation command requires atom_style with both energy and density, e.g. meso");

  if (narg != 5)
    error->all(FLERR,"Illegal number of arguments for fix meso/precipitation command");

  time_integrate = 0;

  // Obtain the neighbor cutoff distance
  int m = 3;
  mAthres = atof(arg[m++]);
  cAeq = atof(arg[m++]);
  if (mAthres <= 0) {
    error->all(FLERR,"Illegal value for mass threshold");
  }

  if (cAeq <= 0) {
    error->all(FLERR,"Illegal value for equilibrium concentration");
  }
  
  // find the concentration property
  int fcA;
  int icA = atom->find_custom("cA", fcA);
  if (icA < 0)
    error->all(FLERR,
        "Can't find property cA for fix meso/precipitationA"); 
  cA = atom->dvector[icA];

  // find the mass of A property
  int fmA;
  int imA = atom->find_custom("mA", fmA);
  if (imA < 0)
    error->all(FLERR,
        "Can't find property mA for fix meso/precipitationA");
  mA = atom->dvector[imA];

  // set comm size by this fix
  comm_reverse = 1;
}

/* ---------------------------------------------------------------------- */

FixMesoPrecipitationA::~FixMesoPrecipitationA() {
}

/* ---------------------------------------------------------------------- */

int FixMesoPrecipitationA::setmask() {
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMesoPrecipitationA::init() {
  // need a full neighbor list, built whenever re-neighboring occurs
  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ---------------------------------------------------------------------- */

void FixMesoPrecipitationA::init_list(int, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixMesoPrecipitationA::end_of_step()
{
  int i, j, ii, jj, itype, jtype, jshortest, inum, jnum;
  int *ilist, *jlist, *numneigh, **firstneigh;

  double delx, dely, delz;
  double shortest, rsq, r;
  
  double **x = atom->x;
  double **v = atom->v;
  int *type = atom->type;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  bool foundShortest;

  // Create memory for ischangecA
  memory->create(ischangecA, nall, "fix:ischangecA");
  
  // Initilise all phase change variables
  for (i=0; i < nall; i++) {
    ischangecA[i] = 0.0;
  }

  for (ii=0; ii<inum; ii++) {
    i = ilist[ii];
    itype = type[i];

    if (i<nlocal) {
      // Only deal with local solid particles
      if (itype == 2) {
        if (mA[i] >= mAthres) { // precipitation
          // Check neighbouring atoms
          jnum = numneigh[i];
          jlist = firstneigh[i];

          // Then need to find the closest fluid particles
          shortest = 100000.0;
          foundShortest = false;

          for (jj = 0; jj < jnum; jj++) {
            j = jlist[jj];
            j &= NEIGHMASK;
            jtype = type[j];

            delx = x[i][0] - x[j][0];
            dely = x[i][1] - x[j][1];
            delz = x[i][2] - x[j][2];
            r = sqrt(delx*delx + dely*dely + delz*delz);

            // To avoid the ghost region effect, limit the range for search
            if ((r<shortest) && (jtype==1)) {
              // check that we are not triggering precipitation in ghost zone
              // Also, to avoid one precipitation triggering many
              foundShortest = true;
              shortest = r;
              jshortest = j;
            }
          } // for loop to find closest fluid

          // if there is a closest liquid particle
          if (foundShortest) {
            // Trigger phase change for i and j
            ischangecA[i] += 1.0;
            ischangecA[jshortest] += 1.0;
          }
        }
        else if (mA[i] <= -mAthres) { // Dissolution
          ischangecA[i] += 1.0;
        }
      }
    }
  }

  // Send the phase change status in reverse first
  comm->reverse_comm_fix(this);

  // Then loop through the local atoms to see if there is a trigger
  for (i = 0; i < nlocal; i++) {
    if (ischangecA[i] > 0.01) {
      // Phase change happening
      if (type[i] == 1) {
        // Changing liquid into solid
        mA[i] = 0.0;
        type[i] = 2; // convert the liquid to the solid
        cA[i] = 0.0; // concentration is 0.0
        //   v[jshortest][0] = 0.0; // set velocity to 0.0
        //   v[jshortest][1] = 0.0;
        //   v[jshortest][2] = 0.0;
      }
      else if (type[i] == 2) {
        // changing solid into liquid
        if (mA[i] < -mAthres) {
          mA[i] = 0.0; // mass becomes 0.0 for dissolution
          type[i] = 1; // convert solid to liquid
          cA[i] = cAeq; // concentration reach back to equilibrium
        }
        else if (mA[i] >= mAthres) {
          // Only causing the particle to reduce the mass
          mA[i] = mA[i] - mAthres;
        }
      }
    }
  }

  // Destroy the memory created for ischangecA
  memory->destroy(ischangecA);
}

/* ---------------------------------------------------------------------- */

int FixMesoPrecipitationA::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;
  int *type = atom->type;
  
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = ischangecA[i];
  }
 return m;
}

/* ---------------------------------------------------------------------- */

void FixMesoPrecipitationA::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;
  int *type = atom->type;
  
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    ischangecA[j] += buf[m++];
  }
}
