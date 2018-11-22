#include <math.h>
#include <stdlib.h>
#include "pair_sph_concAporousLangmuirsurfacereaction.h"
#include "sph_kernel_quintic.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "domain.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairSPHConcAPorousLangmuirSurfaceReaction::PairSPHConcAPorousLangmuirSurfaceReaction(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;

  // find the aqueous mass fraction property
  int fxA;
  int ixA = atom->find_custom("xA", fxA);
  if (ixA < 0)
    error->all(FLERR,
        "Can't find property xA for pair_style sph/concAporousLangmuirsurfacereaction");
  xA = atom->dvector[ixA];

  // find the change in aqueous mass fraction concentration property
  int fdxA;
  int idxA = atom->find_custom("dxA", fdxA);
  if (idxA < 0)
    error->all(FLERR,
        "Can't find property dxA for pair_style sph/concAporousLangmuirsurfacereaction");
  dxA = atom->dvector[idxA];

  // find the absorbed mass fraction property
  int fyA;
  int iyA = atom->find_custom("yA", fyA);
  if (iyA < 0)
    error->all(FLERR,
        "Can't find property yA for pair_style sph/concAporousLangmuirsurfacereaction");
  yA = atom->dvector[iyA];

  // find the change in the absorbed mass fraction
  int fdyA;
  int idyA = atom->find_custom("dyA", fdyA);
  if (idyA < 0)
    error->all(FLERR,
        "Can't find property dyA for pair_style sph/concAporousLangmuirsurfacereaction");
  dyA = atom->dvector[idyA];

  // find the normalised concentration of adsorbed species
  int fthetaA;
  int ithetaA = atom->find_custom("thetaA", fthetaA);
  if (ithetaA < 0)
    error->all(FLERR,
        "Can't find property thetaA for pair_style sph/concAporousLangmuirsurfacereaction");
  thetaA = atom->dvector[ithetaA];

  // find the local diffusivity constant
  int fDA;
  int iDA = atom->find_custom("DA", fDA);
  if (iDA < 0)
    error->all(FLERR,
        "Can't find property DA for pair_style sph/concAporousLangmuirsurfacereaction");
  DA = atom->dvector[iDA];

  // set comm size needed by this pair
  comm_forward = 3;
  comm_reverse = 2;
}

/* ---------------------------------------------------------------------- */

PairSPHConcAPorousLangmuirSurfaceReaction::~PairSPHConcAPorousLangmuirSurfaceReaction() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
  }
}

/* ---------------------------------------------------------------------- */
void PairSPHConcAPorousLangmuirSurfaceReaction::init_style()
{
  // need a full neighbour list
  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ---------------------------------------------------------------------- */

void PairSPHConcAPorousLangmuirSurfaceReaction::compute(int eflag, int vflag) {
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double delx, dely, delz;
  double xNij, yNij, zNij, Nij;

  int *ilist, *jlist, *numneigh, **firstneigh;
  double imass, jmass, h, ih, ihsq;
  double rsq, wf, wfd, D, K, deltaxA;

  double ni, nj;

  if (eflag || vflag)
    ev_setup(eflag, vflag);
  else
    evflag = vflag_fdotr = 0;

  double **x = atom->x;

  double *rmass = atom->rmass;
  double *rho = atom->rho;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double **cg = atom->colorgradient;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // Communicate the local xA to the ghost atoms
  comm->forward_comm_pair(this);

  // loop over neighbors of my atoms and do heat diffusion
  for (ii = 0; ii < inum; ii++)
    {
      dyA[ii] = 0.0;
      dxA[ii] = 0.0;
    }

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];

    // check if the i particles is within the domain
    if ((is_periodic==1) || (not (x[i][0] < domain->boxlo[0] || x[i][0] > domain->boxhi[0] ||
				  x[i][1] < domain->boxlo[1] || x[i][1] > domain->boxhi[1] ||
				  x[i][2] < domain->boxlo[2] || x[i][2] > domain->boxhi[2]))) {
      jlist = firstneigh[i];
      jnum = numneigh[i];

      imass = rmass[i];

      // variable to check if there is a fluid type inside solid
      bool isfluidin = false;
      for (jj = 0; jj < jnum; jj++) {
	j = jlist[jj];
	j &= NEIGHMASK;
	jtype = type[j];
	jmass = rmass[j];

	// check if the j particles is within the domain
	if ((is_periodic==1) || (not (x[j][0] < domain->boxlo[0] || x[j][0] > domain->boxhi[0] ||
				      x[j][1] < domain->boxlo[1] || x[j][1] > domain->boxhi[1] ||
				      x[j][2] < domain->boxlo[2] || x[j][2] > domain->boxhi[2]))) {
	  delx = x[j][0] - x[i][0];
	  dely = x[j][1] - x[i][1];
	  delz = x[j][2] - x[i][2];
	  rsq = delx * delx + dely * dely + delz * delz;

	  // Check if j is within the support kernel
	  if (rsq < cutsq[itype][jtype]) {
	    h = cut[itype][jtype];
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
	    // Diffusion now allowed inside solid as well
	    jmass = rmass[j];
	    // Calculating the particle exchange
	    // Reference: Tartakovsky(2007) - Simulations of reactive transport
	    // and precipitation with sph
	    // The constants give better results...
	    ni = rho[i] / imass;
	    nj = rho[j] / jmass;
	    deltaxA = (1.0/(imass*sqrt(rsq)))*
	      ((DA[i]*ni*imass + DA[j]*nj*jmass)/(ni*nj))*(xA[i] - xA[j])*wfd;
	    dxA[i] = dxA[i] + deltaxA;

	    if ((itype==1) && (jtype==2)) { // fluid-solid interaction
	      jmass = rmass[j];
	      // Calculate the normal vector of colour gradient
	      // Since the colour gradient is pointing AWAY from the surface, this needs to be modified
	      // to match the definition from Ryan's paper
	      xNij = -cg[i][0] + cg[j][0];
	      yNij = -cg[i][1] + cg[j][1];
	      zNij = -cg[i][2] + cg[j][2];
	      // Check if Nijsq is zero
	      double Nijsq = sqrt(xNij*xNij + yNij*yNij + zNij*zNij);
	      Nij = (Nijsq == 0.0) ? 0.0 : xNij*delx + yNij*dely + zNij*delz;
	      // Calculate the exchange in concentration
	      ni = rho[i] / imass;
	      nj = rho[j] / jmass;
	      deltaxA = (kaA*xA[i]*pow((1-thetaA[j]), lambda) - (kdA*pow(thetaA[j], lambda))/(ni*imass))*
		(2.0*(itype-jtype)*Nij*wfd)/(nj+ni);
	      dxA[i] = dxA[i] - deltaxA;
	    } // fluid-solid interaction
	    else if ((itype==2) && (jtype==1)) {
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
	      Nij = (Nijsq == 0.0) ? 0.0 : xNij*delx + yNij*dely + zNij*delz;
	      // Calculate the exchange in concentration
	      ni = rho[i] / imass;
	      nj = rho[j] / jmass;
	      deltaxA = (kaA*xA[j]*pow((1-thetaA[i]), lambda) - (kdA*pow(thetaA[i], lambda))/(ni*imass))*
		(2.0*(itype-jtype)*Nij*wfd)/(nj+ni);
	      dyA[i] = dyA[i] - deltaxA;
	    }
	  } // check i-j within support kernel
	} // check if j particle is inside
      } // jj loop
      // Extra reaction for Darcy's scale model inside porous solid
      if ((not isfluidin) && (itype==2)) {
	// Add the terms for the Darcy's scale model
	dxA[i] -= (As/Vp)*(kaA*xA[i]*pow((1-thetaA[i]), lambda) - (kdA*pow(thetaA[i], lambda))/(ni*imass));
	dyA[i] += (As/Vp)*(kaA*xA[i]*pow((1-thetaA[i]), lambda) - (kdA*pow(thetaA[i], lambda))/(ni*imass));
      }
    } // check i atom is inside domain
  } // ii loop
  // Communicate the ghost dxA and dmA to the locally owned atoms
  comm->reverse_comm_pair(this);
}

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairSPHConcAPorousLangmuirSurfaceReaction::allocate() {
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
  memory->create(cut, n + 1, n + 1, "pair:cut");
}

/* ----------------------------------------------------------------------
 global settings
 ------------------------------------------------------------------------- */

void PairSPHConcAPorousLangmuirSurfaceReaction::settings(int narg, char **arg) {
  if (narg != 0)
    error->all(FLERR,
        "Illegal number of setting arguments for pair_style sph/concAporousLangmuirsurfacereaction");
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairSPHConcAPorousLangmuirSurfaceReaction::coeff(int narg, char **arg) {
  if (narg != 9)
    error->all(FLERR,"Incorrect number of args for pair_style sph/concAporousLangmuirsurfacereaction coefficients");
  if (!allocated)
    allocate();

  int ilo, ihi, jlo, jhi;
  force->bounds(FLERR,arg[0], atom->ntypes, ilo, ihi);
  force->bounds(FLERR,arg[1], atom->ntypes, jlo, jhi);

  // Get kernel size
  double kernel_one = force->numeric(FLERR, arg[2]);

  // Get the parameters for Langmuir adsorption
  kaA = force->numeric(FLERR, arg[3]);
  kdA = force->numeric(FLERR, arg[4]);
  lambda = force->numeric(FLERR, arg[5]);
  As = force->numeric(FLERR, arg[6]);
  Vp = force->numeric(FLERR, arg[7]);

  // Check if periodicity for transport is allowed
  is_periodic = force->numeric(FLERR, arg[8]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut[i][j] = kernel_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0)
    error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
 init for one type pair i,j and corresponding j,i
 ------------------------------------------------------------------------- */

double PairSPHConcAPorousLangmuirSurfaceReaction::init_one(int i, int j) {

  if (setflag[i][j] == 0) {
    error->all(FLERR,"All pair sph/concAporousLangmuirsurfacereaction coeffs are not set");
  }

  cut[j][i] = cut[i][j];

  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

double PairSPHConcAPorousLangmuirSurfaceReaction::single(int i, int j, int itype, int jtype,
    double rsq, double factor_coul, double factor_lj, double &fforce) {
  fforce = 0.0;

  return 0.0;
}

/* ---------------------------------------------------------------------- */

int PairSPHConcAPorousLangmuirSurfaceReaction::pack_forward_comm(int n, int *list, double *buf,
					     int pbc_flag, int *pbc)
{
  int i, j, m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = xA[j];
    buf[m++] = yA[j];
    buf[m++] = thetaA[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairSPHConcAPorousLangmuirSurfaceReaction::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    xA[i] = buf[m++];
    yA[i] = buf[m++];
    thetaA[i] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int PairSPHConcAPorousLangmuirSurfaceReaction::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = dxA[i];
    buf[m++] = dyA[i];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairSPHConcAPorousLangmuirSurfaceReaction::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];

    dxA[j] += buf[m++];
    dyA[j] += buf[m++];
  }
}
