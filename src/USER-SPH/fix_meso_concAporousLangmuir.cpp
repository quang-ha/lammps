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
#include "fix_meso_concAporousLangmuir.h"
#include "sph_kernel_quintic.h"
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
#include "domain.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* TODO: Make this fix run in parallel */

/* ---------------------------------------------------------------------- */

FixMesoConcAPorousLangmuir::FixMesoConcAPorousLangmuir(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg) {

  if ((atom->e_flag != 1) || (atom->rho_flag != 1))
    error->all(FLERR,
        "fix meso/concAporousLangmuir command requires atom_style with both energy and density, e.g. meso");

  if (narg != 4)
    error->all(FLERR,"Illegal number of arguments for fix meso/concAporousLangmuir command");

  // required args
  int m = 3;
  h = atof(arg[m++]);

  time_integrate = 0;

  // find the concentration property
  int fcA;
  int icA = atom->find_custom("cA", fcA);
  if (icA < 0)
    error->all(FLERR,
        "Can't find property cA for fix meso/concAporousLangmuir");
  cA = atom->dvector[icA];

  // find the local concentration property
  int fdcA;
  int idcA = atom->find_custom("dcA", fdcA);
  if (idcA < 0)
    error->all(FLERR,
        "Can't find property dcA for fix meso/concAporousLangmuir");
  dcA = atom->dvector[idcA];

  // find the mass fraction property
  int fyA;
  int iyA = atom->find_custom("yA", fyA);
  if (iyA < 0)
    error->all(FLERR,
        "Can't find property yA for fix meso/concAporousLangmuir");
  yA = atom->dvector[iyA];

  // find the change in the mass fraction
  int fdyA;
  int idyA = atom->find_custom("dyA", fdyA);
  if (idyA < 0)
    error->all(FLERR,
        "Can't find property dyA for fix meso/concAporousLangmuir");
  dyA = atom->dvector[idyA];

  // find the maximum mass fraction property
  int fyAmax;
  int iyAmax = atom->find_custom("yAmax", fyAmax);
  if (iyAmax < 0)
    error->all(FLERR,
        "Can't find property yAmax for fix meso/concAporousLangmuir");
  yAmax = atom->dvector[iyAmax];

  // find the local diffusivity constant
  int fDA;
  int iDA = atom->find_custom("DA", fDA);
  if (iDA < 0)
    error->all(FLERR,
        "Can't find property DA for fix meso/concAporousLangmuir");
  DA = atom->dvector[iDA];

  // find the adsorption rate
  int fkAa;
  int ikAa = atom->find_custom("kAa", fkAa);
  if (ikAa < 0)
    error->all(FLERR,
        "Can't find property kAa for fix meso/concAporousLangmuir");
  kAa = atom->dvector[ikAa];

  // find the desorption rate
  int fkAd;
  int ikAd = atom->find_custom("kAd", fkAd);
  if (ikAd < 0)
    error->all(FLERR,
        "Can't find property kAd for fix meso/concAporousLangmuir");
  kAd = atom->dvector[ikAd];

  // find the normalised concentration of adsorbed species
  int fthetaA;
  int ithetaA = atom->find_custom("thetaA", fthetaA);
  if (ithetaA < 0)
    error->all(FLERR,
        "Can't find property thetaA for fix meso/concAporousLangmuir");
  thetaA = atom->dvector[ithetaA];

  // Find the maximum surface concentration
  int fsA;
  int isA = atom->find_custom("sA", fsA);
  if (isA < 0)
    error->all(FLERR,
        "Can't find property sA for fix meso/concAporousLangmuir");
  sA = atom->dvector[isA];

  // Find the surface area of solid's pore
  int fAs;
  int iAs = atom->find_custom("As", fAs);
  if (iAs < 0)
    error->all(FLERR,
        "Can't find property As for fix meso/concAporousLangmuir");
  As = atom->dvector[iAs];

  // Find the pore volume of solid's pore
  int fVp;
  int iVp = atom->find_custom("Vp", fVp);
  if (iVp < 0)
    error->all(FLERR,
        "Can't find property Vp for fix meso/concAporousLangmuir");
  Vp = atom->dvector[iVp];
}

/* ---------------------------------------------------------------------- */

int FixMesoConcAPorousLangmuir::setmask() {
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMesoConcAPorousLangmuir::init() {
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

void FixMesoConcAPorousLangmuir::init_list(int, NeighList *ptr)
{
  list = ptr;
}

/* ----------------------------------------------------------------------
 allow for both per-type and per-atom mass
 ------------------------------------------------------------------------- */

void FixMesoConcAPorousLangmuir::initial_integrate(int vflag) {
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int *type = atom->type;

  int i;

  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      cA[i] += dtf * dcA[i]; // half-step update of particle concentraion
      if (type[i] == 2) // Only update mass fraction for solid particles
        yA[i] += dtf*dyA[i];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixMesoConcAPorousLangmuir::final_integrate() {
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int *type = atom->type;

  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      cA[i] += dtf*dcA[i];
      if (type[i] == 2) // Only update mass fraction for solid particels
        yA[i] += dtf*dyA[i];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixMesoConcAPorousLangmuir::end_of_step()
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double ni, nj;

  double delx, dely, delz;
  double xtmp, ytmp, ztmp;
  double xNij, yNij, zNij, Nij;

  int *ilist, *jlist, *numneigh, **firstneigh;
  double imass, jmass, ih, ihsq;
  double rsq, wf, wfd;

  double **x = atom->x;
  double **v = atom->v;

  double *rmass = atom->rmass;
  double *rho = atom->rho;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double **cg = atom->colorgradient;

  int *mask = atom->mask;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];

    // Only update maximum mass fraction for solid particles
    if (itype == 2) {
      imass = rmass[i];
      // Keep the position of the solid particle
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];

      // Check neighbouring atoms
      int** firstneigh = list->firstneigh;
      int jnum = numneigh[i];
      int* jlist = firstneigh[i];

      // variable to check if there are fluid particles in the support
      bool isfluidin = false;
      // Then need to find the closest fluid particles
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;
        jtype = type[j];

        // Calculate the distance between the particles
        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx * delx + dely * dely + delz * delz;

        // Check if j is within the support kernel
        if (rsq < h) {
          ih = 1.0/h;

          // kernel function
          if (domain->dimension == 3) {
            wfd = sph_dw_quintic3d(sqrt(rsq)*ih);
            wfd = wfd*ih*ih*ih*ih;
            wf = sph_kernel_quintic3d(sqrt(rsq)*ih)*ih*ih*ih;
            } else {
            wfd = sph_dw_quintic2d(sqrt(rsq)*ih);
            wfd = wfd*ih*ih*ih;
            wf = sph_kernel_quintic2d(sqrt(rsq)*ih)*ih*ih;
          }

          // Perform interaction calculation
          if (jtype == 1) { // only for fluid particles
            isfluidin = true;
            jmass = rmass[j];
            // Calculate the normal vector of colour gradient
            // Since the colour gradient is pointing AWAY from the surface, this needs to be modified
            // to match the definition from Ryan's paper
            xNij = -cg[i][0] + cg[j][0];
            yNij = -cg[i][1] + cg[j][1];
            zNij = -cg[i][2] + cg[j][2];
            // Check if Nijsq is zero
            double Nijsq = sqrt(xNij*xNij + yNij*yNij + zNij*zNij);
            Nij = (Nijsq == 0.0) ? 0.0 : (xNij*delx + yNij*dely + zNij*delz);
            // Calculate the exchange in concentration
            ni = rho[i] / imass;
            nj = rho[j] / jmass;
            yAmax[i] = yAmax[i] + (Nij*wfd*sA[i])/(ni*nj*imass);
          } // jtype fluid
        } // loop inside support kernel
      } // for loop jj
      // Calculate the absorbed concentration
      if (isfluidin)
        thetaA[i] = (yAmax[i] == 0.0) ? 0.0 : (yA[i]/yAmax[i]);
      else {
        // within solid domain, use Darcy-scale model
        yAmax[i] = sA[i]/rho[i];
        thetaA[i] = yA[i]/yAmax[i];
      }
    } // if i type is solid
  } // loop through i
  // comm->forward_comm_fix(this);
  // comm->reverse_comm_fix(this);
}

/* ---------------------------------------------------------------------- */

void FixMesoConcAPorousLangmuir::reset_dt() {
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}

/* ---------------------------------------------------------------------- */
