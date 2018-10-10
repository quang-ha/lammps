#include <math.h>
#include <stdlib.h>
#include "pair_sph_concALangmuirsurfacereaction.h"
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

PairSPHConcALangmuirSurfaceReaction::PairSPHConcALangmuirSurfaceReaction(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;

  // find the concentration property
  int fcA;
  int icA = atom->find_custom("cA", fcA);
  if (icA < 0)
    error->all(FLERR,
        "Can't find property cA for pair_style sph/concALangmuirsurfacereaction");
  cA = atom->dvector[icA];
  
  // find the local concentration property
  int fdcA;
  int idcA = atom->find_custom("dcA", fdcA);
  if (idcA < 0)
    error->all(FLERR,
        "Can't find property dcA for pair_style sph/concALangmuirsurfacereaction");
  dcA = atom->dvector[idcA];

  // find the mass fraction property
  int fyA;
  int iyA = atom->find_custom("yA", fyA);
  if (iyA < 0)
    error->all(FLERR,
        "Can't find property yA for pair_style sph/concALangmuirsurfacereaction");
  yA = atom->dvector[iyA];

  // find the change in the mass fraction
  int fdyA;
  int idyA = atom->find_custom("dyA", fdyA);
  if (idyA < 0)
    error->all(FLERR,
        "Can't find property dyA for pair_style sph/concALangmuirsurfacereaction");
  dyA = atom->dvector[idyA];
  
  // find the maximum mass fraction property
  int fyAmax;
  int iyAmax = atom->find_custom("yAmax", fyAmax);
  if (iyAmax < 0)
    error->all(FLERR,
        "Can't find property yAmax for pair_style sph/concALangmuirsurfacereaction");
  yAmax = atom->dvector[iyAmax];
  
  // find the local diffusivity constant
  int fDA;
  int iDA = atom->find_custom("DA", fDA);
  if (iDA < 0)
    error->all(FLERR,
        "Can't find property DA for pair_style sph/concALangmuirsurfacereaction");
  DA = atom->dvector[iDA];

  // find the adsorption rate
  int fkAa;
  int ikAa = atom->find_custom("kAa", fkAa);
  if (ikAa < 0)
    error->all(FLERR,
        "Can't find property kAa for pair_style sph/concALangmuirsurfacereaction");
  kAa = atom->dvector[ikAa];

  // find the desorption rate
  int fkAd;
  int ikAd = atom->find_custom("kAd", fkAd);
  if (ikAd < 0)
    error->all(FLERR,
        "Can't find property kAd for pair_style sph/concALangmuirsurfacereaction");
  kAd = atom->dvector[ikAd];
    
  // find the normalised concentration of adsorbed species
  int fthetaA;
  int ithetaA = atom->find_custom("thetaA", fthetaA);
  if (ithetaA < 0)
    error->all(FLERR,
        "Can't find property thetaA for pair_style sph/concALangmuirsurfacereaction");
  thetaA = atom->dvector[ithetaA];

  // Find the maximum surface concentration
  int fsA;
  int isA = atom->find_custom("sA", fsA);
  if (isA < 0)
    error->all(FLERR,
        "Can't find property sA for pair_style sph/concALangmuirsurfacereaction");
  sA = atom->dvector[isA];

  // set comm size needed by this pair
  comm_forward = 5;
  comm_reverse = 3;
}

/* ---------------------------------------------------------------------- */

PairSPHConcALangmuirSurfaceReaction::~PairSPHConcALangmuirSurfaceReaction() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
  }
}

/* ---------------------------------------------------------------------- */
void PairSPHConcALangmuirSurfaceReaction::init_style()
{
  // need a full neighbour list
  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ---------------------------------------------------------------------- */

void PairSPHConcALangmuirSurfaceReaction::compute(int eflag, int vflag) {
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz;
  double xNij, yNij, zNij, Nij;
  
  int *ilist, *jlist, *numneigh, **firstneigh;
  double imass, jmass, h, ih, ihsq;
  double rsq, wf, wfd, D, K, deltacA;

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
  
  // loop over neighbors of my atoms and do heat diffusion
  for (ii = 0; ii < inum; ii++)
    {
      dyA[ii] = 0.0;
      dcA[ii] = 0.0;
    }

  // Communicate the local cA to the ghost atoms
  comm->forward_comm_pair(this);
  
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];	

    // check if the i particles is within the domain
    if (not (xtmp < domain->boxlo[0] || xtmp > domain->boxhi[0] ||
	     ytmp < domain->boxlo[1] || ytmp > domain->boxhi[1] ||
	     ztmp < domain->boxlo[2] || ztmp > domain->boxhi[2]))
       { 
	 if (itype == 1) // only do diffusion for ifluid particle
	   {
	     jlist = firstneigh[i];
	     jnum = numneigh[i];
	     
	     imass = rmass[i];
	     
	     for (jj = 0; jj < jnum; jj++) {
	       j = jlist[jj];
	       j &= NEIGHMASK;
	       jtype = type[j];

	       // check if the j particles is within the domain
	       if (not (x[j][0] < domain->boxlo[0] || x[j][0] > domain->boxhi[0] ||
			x[j][1] < domain->boxlo[1] || x[j][1] > domain->boxhi[1] ||
			x[j][2] < domain->boxlo[2] || x[j][2] > domain->boxhi[2]))
		 {
		   delx = x[j][0] - xtmp;
		   dely = x[j][1] - ytmp;
		   delz = x[j][2] - ztmp;
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
		     if (jtype == 1) { // this is only for fluid-fluid interaction 			 
		       jmass = rmass[j];
		       
		       // Calculating the particle exchange
		       // Reference: Tartakovsky(2007) - Simulations of reactive transport
		       // and precipitation with sph
		       // The constants give better results...
		       ni = rho[i] / imass;
		       nj = rho[j] / jmass;
		       deltacA = (1.0/(imass*sqrt(rsq)))*
			 ((DA[i]*ni*imass + DA[j]*nj*jmass)/(ni*nj))*(cA[i] - cA[j])*wfd;
		       dcA[i] = dcA[i] + deltacA;
		     }
		     else { // if jtype is solid
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
		       deltacA = (kAa[i]*cA[i]*(1-thetaA[j])*(1-thetaA[j]) - (kAd[i]*thetaA[j]*thetaA[j])/(ni*imass))*
		         (2.0*(itype-jtype)*Nij*wfd)/(nj+ni);
		       dcA[i] = dcA[i] - deltacA;
		       dyA[j] = dyA[j] + deltacA;
		     } // jtype solid
		   } // check within support kernel
		 } // check if j particle is inside
	     } // jj loop
	   } //itype fluid
       } // check i atom is inside domain
  } // ii loop
  // Communicate the ghost dcA and dmA to the locally owned atoms
  comm->reverse_comm_pair(this);
}

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairSPHConcALangmuirSurfaceReaction::allocate() {
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

void PairSPHConcALangmuirSurfaceReaction::settings(int narg, char **arg) {
  if (narg != 0)
    error->all(FLERR,
        "Illegal number of setting arguments for pair_style sph/concALangmuirsurfacereaction");
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairSPHConcALangmuirSurfaceReaction::coeff(int narg, char **arg) {
  if (narg != 3)
    error->all(FLERR,"Incorrect number of args for pair_style sph/concALangmuirsurfacereaction coefficients");
  if (!allocated)
    allocate();

  int ilo, ihi, jlo, jhi;
  force->bounds(FLERR,arg[0], atom->ntypes, ilo, ihi);
  force->bounds(FLERR,arg[1], atom->ntypes, jlo, jhi);
  
  double cut_one = force->numeric(FLERR,arg[2]);
  
  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut[i][j] = cut_one;
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

double PairSPHConcALangmuirSurfaceReaction::init_one(int i, int j) {

  if (setflag[i][j] == 0) {
    error->all(FLERR,"All pair sph/concALangmuirsurfacereaction coeffs are not set");
  }

  cut[j][i] = cut[i][j];
  
  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

double PairSPHConcALangmuirSurfaceReaction::single(int i, int j, int itype, int jtype,
    double rsq, double factor_coul, double factor_lj, double &fforce) {
  fforce = 0.0;

  return 0.0;
}

/* ---------------------------------------------------------------------- */

int PairSPHConcALangmuirSurfaceReaction::pack_forward_comm(int n, int *list, double *buf,
					     int pbc_flag, int *pbc)
{
  int i, j, m;
  int *type = atom->type;
  
  m = 0;
  for (i = 0; i < n; i++)
    {
      j = list[i];
      buf[m++] = cA[j];
      buf[m++] = dcA[j];
      buf[m++] = yA[j];
      buf[m++] = dyA[j];
    }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairSPHConcALangmuirSurfaceReaction::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;
  int *type = atom->type;
  
  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    {
      cA[i] = buf[m++];
      dcA[i] = buf[m++];
      yA[i] = buf[m++];
      dyA[i] = buf[m++];
    }
}

/* ---------------------------------------------------------------------- */

int PairSPHConcALangmuirSurfaceReaction::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;
  int *type = atom->type;
  
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = dcA[i];
    buf[m++] = dyA[i];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairSPHConcALangmuirSurfaceReaction::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;
  int *type = atom->type;
  
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];

    dcA[j] += buf[m++];
    dyA[j] += buf[m++];
  }
}
