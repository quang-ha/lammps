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

/* ---------------------------------------------------------------------- */

FixMesoConcAPorousLangmuir::FixMesoConcAPorousLangmuir(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg) {

  if ((atom->e_flag != 1) || (atom->rho_flag != 1))
    error->all(FLERR,
        "fix meso/concAporousLangmuir command requires atom_style with both energy and density, e.g. meso");

  if (narg != 5)
    error->all(FLERR,"Illegal number of arguments for fix meso/concAporousLangmuir command");

  // required args
  int m = 3;
  h = atof(arg[m++]);
  sAmax = atof(arg[m++]);

  time_integrate = 0;

  // find the aqueous mass fraction property
  int fxA;
  int ixA = atom->find_custom("xA", fxA);
  if (ixA < 0)
    error->all(FLERR,
        "Can't find property xA for fix meso/concAporousLangmuir");
  xA = atom->dvector[ixA];

  // find the change in aqueous mass fraction concentration property
  int fdxA;
  int idxA = atom->find_custom("dxA", fdxA);
  if (idxA < 0)
    error->all(FLERR,
        "Can't find property dxA for fix meso/concAporousLangmuir");
  dxA = atom->dvector[idxA];

  // find the absorbed mass fraction property
  int fyA;
  int iyA = atom->find_custom("yA", fyA);
  if (iyA < 0)
    error->all(FLERR,
        "Can't find property yA for fix meso/concAporousLangmuir");
  yA = atom->dvector[iyA];

  // find the change in the absorbed mass fraction
  int fdyA;
  int idyA = atom->find_custom("dyA", fdyA);
  if (idyA < 0)
    error->all(FLERR,
        "Can't find property dyA for fix meso/concAporousLangmuir");
  dyA = atom->dvector[idyA];

  // find the normalised concentration of adsorbed species
  int fthetaA;
  int ithetaA = atom->find_custom("thetaA", fthetaA);
  if (ithetaA < 0)
    error->all(FLERR,
        "Can't find property thetaA for fix meso/concAporousLangmuir");
  thetaA = atom->dvector[ithetaA];
}

/* ---------------------------------------------------------------------- */

int FixMesoConcAPorousLangmuir::setmask() {
  int mask = 0;
  mask |= FINAL_INTEGRATE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMesoConcAPorousLangmuir::init() {
  dtxA = update->dt;

  // need a full neighbor list, built whenever re-neighboring occurs
  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ---------------------------------------------------------------------- */

void FixMesoConcAPorousLangmuir::init_list(int, NeighList *ptr) {
  list = ptr;
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
      xA[i] += dtxA*dxA[i];
      if (type[i] == 2) // Only update mass fraction for solid particels
        yA[i] += dtxA*dyA[i];
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
  double yAmax;

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

      // Init value of yAmax to 0
      yAmax = 0.0;
      // Variable to check if there are fluid particles in the support
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
            yAmax = yAmax + (Nij*wfd*sAmax)/(ni*nj*imass);
          } // jtype fluid
        } // loop inside support kernel
      } // for loop jj
      // Calculate the absorbed concentration
      if (isfluidin)
	thetaA[i] = (yAmax == 0.0) ? 0.0 : (yA[i]/yAmax);
      else
	thetaA[i] = yA[i]/(sAmax/rho[i]);
      // printf("i %d thetaA[i] %f ymax %f \n", i, thetaA[i], yAmax);
    } // if i type is solid
  } // loop through i
}

/* ---------------------------------------------------------------------- */

void FixMesoConcAPorousLangmuir::reset_dt() {
  dtxA = update->dt;
}

/* ---------------------------------------------------------------------- */
