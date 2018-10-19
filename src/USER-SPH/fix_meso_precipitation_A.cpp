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

  if (narg != 3)
    error->all(FLERR,"Illegal number of arguments for fix meso/precipitation command");

  time_integrate = 0;
  
  // find the concentration property
  int fcA;
  int icA = atom->find_custom("cA", fcA);
  if (icA < 0)
    error->all(FLERR,
        "Can't find property cA for fix meso/precipitationA"); 
  cA = atom->dvector[icA];

  // find the equlibrium concentration property
  int fcAeq;
  int icAeq = atom->find_custom("cAeq", fcAeq);
  if (icAeq < 0)
    error->all(FLERR,
        "Can't find property cAeq for fix meso/precipitationA");
  cAeq = atom->dvector[icAeq];
  
  // find the concentration property
  int fdcA;
  int idcA = atom->find_custom("dcA", fdcA);
  if (icA < 0)
    error->all(FLERR,
        "Can't find property dcA for fix meso/precipitationA"); 
  dcA = atom->dvector[idcA];

  // find the solid-liquid interaction
  int fRA;
  int iRA = atom->find_custom("RA", fRA);
  if (iRA < 0)
    error->all(FLERR,
        "Can't find property RA for fix meso/precipitationA");
  RA = atom->dvector[iRA];

  // find the solid-liquid interaction
  int fdmA;
  int idmA = atom->find_custom("dmA", fdmA);
  if (idmA < 0)
    error->all(FLERR,
        "Can't find property dmA for fix meso/precipitationA");
  dmA = atom->dvector[idmA];

  // find the mass of A property
  int fmA;
  int imA = atom->find_custom("mA", fmA);
  if (imA < 0)
    error->all(FLERR,
        "Can't find property mA for fix meso/precipitationA");
  mA = atom->dvector[imA];
  
  // find the mass threshold property
  int fmAthres;
  int imAthres = atom->find_custom("mAthres", fmAthres);
  if (imAthres < 0)
    error->all(FLERR,
        "Can't find property mAthres for fix meso/precipitationA");
  mAthres = atom->dvector[imAthres];

  // Allocate memory for the checking on the phase change trigger
  allocate();

  // set comm size by this fix
  comm_forward = 5;
  comm_reverse = 5;
}

/* ---------------------------------------------------------------------- */

FixMesoPrecipitationA::~FixMesoPrecipitationA() {
  destroy();
}

/* ---------------------------------------------------------------------- */

int FixMesoPrecipitationA::setmask() {
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMesoPrecipitationA::init() {
  dtv = update->dt;
  dtf = 0.5*update->dt*force->ftm2v;

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

/* ----------------------------------------------------------------------
 allow for both per-type and per-atom mass
 ------------------------------------------------------------------------- */

void FixMesoPrecipitationA::initial_integrate(int vflag) {  
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int *type = atom->type;
  
  int i;

  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      cA[i] += dtf * dcA[i]; // half-step update of particle precipitation
      if (type[i] == 2) // Only update mass for solid particles
	mA[i] += dtf*dmA[i];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixMesoPrecipitationA::final_integrate() {
  
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int *type = atom->type;
  
  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      cA[i] += dtf*dcA[i];
      if (type[i] == 2) // Only update mass for solid particles
	mA[i] += dtf*dmA[i];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixMesoPrecipitationA::end_of_step()
{
  int i, j, ii, jj, itype, jtype, jshortest;

  double delx, dely, delz;
  double shortest, xtmp, ytmp, ztmp, rsq, r;
  
  double **x = atom->x;
  double **v = atom->v;
  int *type = atom->type;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int* numneigh = list->numneigh;

  bool foundShortest;
  
  for (i = 0; i < nlocal; i++) {
    itype = type[i];

    // Only deal with solid particles
    if (itype == 2)
      {
	if (mA[i] >= mAthres[i]) // precipitation
	  {
	    // Keep the position of the solid particle
	    xtmp = x[i][0];
	    ytmp = x[i][1];
	    ztmp = x[i][2];
	    
	    // Check neighbouring atoms
	    int** firstneigh = list->firstneigh;
	    int jnum = numneigh[i];
	    int* jlist = firstneigh[i];
    	    
	    // Then need to find the closest fluid particles
	    shortest = 1000.0;
	    foundShortest = false;

	    for (jj = 0; jj < jnum; jj++) {
	      j = jlist[jj];
	      j &= NEIGHMASK;	    
	      jtype = type[j];
	      
	      if (jtype == 1) // if liquid particle then check for shortest distance
		{
		  delx = xtmp - x[j][0];
		  dely = ytmp - x[j][1];
		  delz = ztmp - x[j][2];
		  r = sqrt(delx*delx + dely*dely + delz*delz);
		  
		  if (r < shortest)
		    {
		      foundShortest = true;
		      shortest = r;
		      jshortest = j;
		    }
		} // ifliquid
	    } // for loop to find closest fluid
	    
	    // if there is a closest liquid particle
	    if (foundShortest)
	      {
		// printf("WARNING: Precipitation particles! \n");
	        mA[i] = mA[i] - mAthres[i];
		mA[jshortest] = 0.0;
		type[jshortest] = 2; // convert the liquid to the solid
		cA[jshortest] = 0.0; // concentration is 0.0
		dcA[jshortest] = 0.0; // change in concentration is also 0.0
		v[jshortest][0] = 0.0; // set velocity to 0.0
		v[jshortest][1] = 0.0;
		v[jshortest][2] = 0.0;
	      }
	  }
	else if (mA[i] <= -mAthres[i]) // Dissolution
	  {
	    // printf("WARNING: Dissolution particles! \n");
	    mA[i] = 0.0; // mass becomes 0.0
	    dmA[i] = 0.0; // change in mass also becomes 0.0
	    type[i] = 1; // convert solid to liquid
	    cA[i] = cAeq[i]; // concentration reach back to equilibrium
	    dcA[i] = 0.0; // change of concentration to 0.0
	  }
      }
  }
  // comm->forward_comm_fix(this);
  // comm->reverse_comm_fix(this);
}

/* ---------------------------------------------------------------------- */

void FixMesoPrecipitationA::reset_dt() {
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}

/* ---------------------------------------------------------------------- */

void FixMesoPrecipitationA::allocate() {
  int nall = atom->nlocal + atom->nghost;
  memory->create(ischangecA, nall, "fix:ischangecA");
}

/* ---------------------------------------------------------------------- */

void FixMesoPrecipitationA::destroy() {
  memory->destroy(ischangecA);
}

/* ---------------------------------------------------------------------- */

int FixMesoPrecipitationA::pack_forward_comm(int n, int *list, double *buf,
					     int pbc_flag, int *pbc)
{
  int i, j, m;
  int *type = atom->type;
  
  m = 0;
  for (i = 0; i < n; i++)
    {
      j = list[i];
      buf[m++] = mA[j];
      buf[m++] = dmA[j];
      buf[m++] = cA[j];
      buf[m++] = dcA[j];
      buf[m++] = type[j];
    }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixMesoPrecipitationA::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;
  int *type = atom->type;
  
  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    {
      mA[i] = buf[m++];
      dmA[i] = buf[m++];
      cA[i] = buf[m++];
      dcA[i] = buf[m++];
    }
}

/* ---------------------------------------------------------------------- */

int FixMesoPrecipitationA::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;
  int *type = atom->type;
  
  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    {
      buf[m++] = dmA[i];
      buf[m++] = dcA[i];
    }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixMesoPrecipitationA::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;
  int *type = atom->type;
  
  m = 0;
  for (i = 0; i < n; i++)
    {
      j = list[i];
      dmA[j] = buf[m++];
      dcA[j] = buf[m++];
    }
}
