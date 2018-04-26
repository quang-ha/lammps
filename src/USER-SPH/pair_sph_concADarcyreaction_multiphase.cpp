#include <math.h>
#include <stdlib.h>
#include "pair_sph_concADarcyreaction_multiphase.h"
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

PairSPHConcADarcyReactionMultiPhase::PairSPHConcADarcyReactionMultiPhase(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;

  // find the concentration property
  int fcA;
  int icA = atom->find_custom("cA", fcA);
  if (icA < 0)
    error->all(FLERR,
        "Can't find property cA for pair_style sph/concADarcyreaction/multiphase");
  cA = atom->dvector[icA];

  // find the equlibrium concentration property
  int fcAeq;
  int icAeq = atom->find_custom("cAeq", fcAeq);
  if (icAeq < 0)
    error->all(FLERR,
        "Can't find property cAeq for pair_style sph/concADarcyreaction/multiphase");
  cAeq = atom->dvector[icAeq];
  
  // find the local concentration property
  int fdcA;
  int idcA = atom->find_custom("dcA", fdcA);
  if (idcA < 0)
    error->all(FLERR,
        "Can't find property dcA for pair_style sph/concADarcyreaction/multiphase");
  dcA = atom->dvector[idcA];

  // find the local effective Darcy-scale diffusivity constant
  int fDAeff;
  int iDAeff = atom->find_custom("DAeff", fDAeff);
  if (iDAeff < 0)
    error->all(FLERR,
        "Can't find property DAeff for pair_style sph/concADarcyreaction/multiphase");
  DAeff = atom->dvector[iDAeff];

  // find the effective Darcy-scale reaction constant
  int fkAeff;
  int ikAeff = atom->find_custom("kAeff", fkAeff);
  if (ikAeff < 0)
    error->all(FLERR,
        "Can't find property kAeff for pair_style sph/concADarcyreaction/multiphase");
  kAeff = atom->dvector[ikAeff];
  
  // set comm size needed by this pair
  comm_forward = 2;
  comm_reverse = 1;
}

/* ---------------------------------------------------------------------- */

PairSPHConcADarcyReactionMultiPhase::~PairSPHConcADarcyReactionMultiPhase() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(phasecut);
  }
}

/* ---------------------------------------------------------------------- */
void PairSPHConcADarcyReactionMultiPhase::init_style()
{
  // need a full neighbour list
  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ---------------------------------------------------------------------- */

void PairSPHConcADarcyReactionMultiPhase::compute(int eflag, int vflag) {
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz;

  int *ilist, *jlist, *numneigh, **firstneigh;
  double imass, jmass, h, d, ih, ihsq;
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

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
  // loop over neighbors of my atoms and do heat diffusion
  for (ii = 0; ii < inum; ii++)
    {
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
	 jlist = firstneigh[i];
	 jnum = numneigh[i];
	 
	 imass = rmass[i];
	 
	 for (jj = 0; jj < jnum; jj++) {
	   j = jlist[jj];
	   j &= NEIGHMASK;
	   
	   // check if the j particles is within the domain
	   if (not (x[j][0] < domain->boxlo[0] || x[j][0] > domain->boxhi[0] ||
		    x[j][1] < domain->boxlo[1] || x[j][1] > domain->boxhi[1] ||
		    x[j][2] < domain->boxlo[2] || x[j][2] > domain->boxhi[2]))
	     {
	       delx = xtmp - x[j][0];
	       dely = ytmp - x[j][1];
	       delz = ztmp - x[j][2];
	       rsq = delx * delx + dely * dely + delz * delz;
	       jtype = type[j];
	       
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
		 } // if kernel
		 
		 jmass = rmass[j];
		 
		 // Calculating the particle exchange
		 // Reference: Tartakovsky(2007) - Simulations of reactive transport
		 // and precipitation with sph
		 // The constants give better results...
		 ni = rho[i] / imass;
		 nj = rho[j] / jmass;
		 deltacA = (1.0/(sqrt(rsq)*imass))*((DAeff[i]*ni*imass + DAeff[j]*nj*jmass)/(ni*nj))*(cA[i] - cA[j])*wfd
		   - kAeff[i]*(cA[i] - cAeq[i]);;
		 dcA[i] = dcA[i] + deltacA;
	       } // check if j type is within cutoff distance
	     } // check if j particle is inside box
	 } // loop through j atoms
       } // check i atom is inside domain
  } // ii loop
  
  // Communicate the ghost dcA and dmA to the locally owned atoms
  comm->reverse_comm_pair(this);
}

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairSPHConcADarcyReactionMultiPhase::allocate() {
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
  memory->create(cut, n + 1, n + 1, "pair:cut");
  memory->create(phasecut, n + 1, n + 1, "pair:phasecut");
}

/* ----------------------------------------------------------------------
 global settings
 ------------------------------------------------------------------------- */

void PairSPHConcADarcyReactionMultiPhase::settings(int narg, char **arg) {
  if (narg != 0)
    error->all(FLERR,
        "Illegal number of setting arguments for pair_style sph/concprecipitation");
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairSPHConcADarcyReactionMultiPhase::coeff(int narg, char **arg) {
  if (narg != 4)
    error->all(FLERR,"Incorrect number of args for pair_style sph/concprecipitation coefficients");
  if (!allocated)
    allocate();

  int ilo, ihi, jlo, jhi;
  force->bounds(FLERR,arg[0], atom->ntypes, ilo, ihi);
  force->bounds(FLERR,arg[1], atom->ntypes, jlo, jhi);
  
  double cut_one = force->numeric(FLERR,arg[2]);
  double phasecut_one = force->numeric(FLERR,arg[3]);
  
  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut[i][j] = cut_one;
      phasecut[i][j] = phasecut_one;
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

double PairSPHConcADarcyReactionMultiPhase::init_one(int i, int j) {

  if (setflag[i][j] == 0) {
    error->all(FLERR,"All pair sph/concADarcyreaction/multiphase coeffs are not set");
  }

  cut[j][i] = cut[i][j];
  phasecut[j][i] = phasecut[i][j];
  
  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

double PairSPHConcADarcyReactionMultiPhase::single(int i, int j, int itype, int jtype,
    double rsq, double factor_coul, double factor_lj, double &fforce) {
  fforce = 0.0;

  return 0.0;
}

/* ---------------------------------------------------------------------- */

int PairSPHConcADarcyReactionMultiPhase::pack_forward_comm(int n, int *list, double *buf,
					     int pbc_flag, int *pbc)
{
  int i, j, m;

  m = 0;
  for (i = 0; i < n; i++)
    {
      j = list[i];
      buf[m++] = cA[j];
      buf[m++] = dcA[j];
    }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairSPHConcADarcyReactionMultiPhase::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    {
      cA[i] = buf[m++];
      dcA[i] = buf[m++];
    }
}

/* ---------------------------------------------------------------------- */

int PairSPHConcADarcyReactionMultiPhase::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = dcA[i];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairSPHConcADarcyReactionMultiPhase::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    dcA[j] += buf[m++];
  }
}
