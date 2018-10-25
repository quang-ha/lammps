#include <math.h>
#include <stdlib.h>
#include "pair_sph_concAsurfacereaction_multiphase.h"
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

PairSPHConcASurfaceReactionMultiPhase::PairSPHConcASurfaceReactionMultiPhase(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;

  // find the concentration property
  int fcA;
  int icA = atom->find_custom("cA", fcA);
  if (icA < 0)
    error->all(FLERR,
        "Can't find property cA for pair_style sph/concAsurfacereaction/multiphase");
  cA = atom->dvector[icA];

  // find the equlibrium concentration property
  int fcAeq;
  int icAeq = atom->find_custom("cAeq", fcAeq);
  if (icAeq < 0)
    error->all(FLERR,
        "Can't find property cAeq for pair_style sph/concAsurfacereaction/multiphase");
  cAeq = atom->dvector[icAeq];
  
  // find the local concentration property
  int fdcA;
  int idcA = atom->find_custom("dcA", fdcA);
  if (idcA < 0)
    error->all(FLERR,
        "Can't find property dcA for pair_style sph/concAsurfacereaction/multiphase");
  dcA = atom->dvector[idcA];

  // find the local diffusivity constant
  int fDA;
  int iDA = atom->find_custom("DA", fDA);
  if (iDA < 0)
    error->all(FLERR,
        "Can't find property DA for pair_style sph/concAsurfacereaction/multiphase");
  DA = atom->dvector[iDA];

  // find the solid-liquid interaction
  int fRA;
  int iRA = atom->find_custom("RA", fRA);
  if (iRA < 0)
    error->all(FLERR,
        "Can't find property RA for pair_style sph/concAsurfacereaction/multiphase");
  RA = atom->dvector[iRA];

  // find the change in mass of A property
  int fdmA;
  int idmA = atom->find_custom("dmA", fdmA);
  if (idmA < 0)
    error->all(FLERR,
        "Can't find property dmA for pair_style sph/concAsurfacereaction/multiphase");
  dmA = atom->dvector[idmA];

  // find the mass of A property
  int fmA;
  int imA = atom->find_custom("mA", fmA);
  if (imA < 0)
    error->all(FLERR,
        "Can't find property mA for pair_style sph/concAsurfacereaction/multiphase");
  mA = atom->dvector[imA];
  
  // find the mass threshold property
  int fmAthres;
  int imAthres = atom->find_custom("mAthres", fmAthres);
  if (imAthres < 0)
    error->all(FLERR,
        "Can't find property mAthres for pair_style sph/concAsurfacereaction/multiphase");
  mAthres = atom->dvector[imAthres];

  // set comm size needed by this pair
  comm_forward = 3;
  comm_reverse = 2;
}

/* ---------------------------------------------------------------------- */

PairSPHConcASurfaceReactionMultiPhase::~PairSPHConcASurfaceReactionMultiPhase() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(phase_support);
  }
}

/* ---------------------------------------------------------------------- */
void PairSPHConcASurfaceReactionMultiPhase::init_style()
{
  // need a full neighbour list
  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ---------------------------------------------------------------------- */

void PairSPHConcASurfaceReactionMultiPhase::compute(int eflag, int vflag) {
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double delx, dely, delz;

  int *ilist, *jlist, *numneigh, **firstneigh;
  double imass, jmass, h, ih, ihsq;
  double r, wf, wfd, D, K, deltacA;

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
  int nall = nlocal + atom->nghost;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // Communicate the local cA to the ghost atoms
  comm->forward_comm_pair(this);
  
  // loop over neighbors of my atoms and do heat diffusion
  for (ii = 0; ii < inum; ii++)
    {
      dmA[ii] = 0.0;
      dcA[ii] = 0.0;
    }
  
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    // check that we are only doing local and ghost atoms only
    itype = type[i];
    // check if the i particles is within the domain
    if (not (x[i][0] < domain->boxlo[0] || x[i][0] > domain->boxhi[0] ||
             x[i][1] < domain->boxlo[1] || x[i][1] > domain->boxhi[1] ||
             x[i][2] < domain->boxlo[2] || x[i][2] > domain->boxhi[2])) {
      jlist = firstneigh[i];
      jnum = numneigh[i];

      imass = rmass[i];

      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;
        // check that we are only doing local and ghost atoms only
        jtype = type[j];

        // check if the j particles is within the domain
        if (not (x[j][0] < domain->boxlo[0] || x[j][0] > domain->boxhi[0] ||
                 x[j][1] < domain->boxlo[1] || x[j][1] > domain->boxhi[1] ||
                 x[j][2] < domain->boxlo[2] || x[j][2] > domain->boxhi[2])) {
          delx = x[i][0] - x[j][0];
          dely = x[i][1] - x[j][1];
          delz = x[i][2] - x[j][2];
          r = sqrt(delx * delx + dely * dely + delz * delz);

          if (r < cut[itype][jtype]) {
            h = cut[itype][jtype];
            ih = 1.0/h;

            // kernel function
            if (domain->dimension == 3) {
              wfd = sph_dw_quintic3d(r*ih);
              wfd = wfd*ih*ih*ih*ih;
              wf = sph_kernel_quintic3d(r*ih)*ih*ih*ih;
            } else {
              wfd = sph_dw_quintic2d(r*ih);
              wfd = wfd*ih*ih*ih;
              wf = sph_kernel_quintic2d(r*ih)*ih*ih;
            }

            if ((itype==1) && (jtype==1)) { // fluid-fluid interaction
              jmass = rmass[j];
              // Calculating the particle exchange
              // Reference: Tartakovsky(2007) - Simulations of reactive transport
              // and precipitation with sph
              // The constants give better results...
              ni = rho[i] / imass;
              nj = rho[j] / jmass;
              deltacA = (1.0/(imass*r))*
                ((DA[i]*ni*imass + DA[j]*nj*jmass)/(ni*nj))*(cA[i] - cA[j])*wfd;
              dcA[i] = dcA[i] + deltacA;
            } // fluid-fluid interaction
            else if ((itype==1) && (jtype==2)) { // fluid-solid interaction
              if (r <= phase_support[itype][jtype]) {
                deltacA = 1.0*RA[i]*(cA[i] - cAeq[i]);
                dcA[i] = dcA[i] - deltacA;
              }
            } // fluid-solid interaction
            else if ((itype==2) && (jtype==1)) { // solid-fluid interaction
              if (r <= phase_support[itype][jtype]) {
                dmA[i] = dmA[i] + (imass + mA[i])*RA[i]*(cA[j] - cAeq[j]);
              }
            } // solid-fluid interaction
          } // check if j particle is inside kernel
        } // check if j particle is inside domain
      } // jj loop
    } // check i atom is inside domain
  } // ii loop

  // Communicate the ghost dcA and dmA to the locally owned atoms
  comm->reverse_comm_pair(this);
}

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairSPHConcASurfaceReactionMultiPhase::allocate() {
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
  memory->create(cut, n + 1, n + 1, "pair:cut");
  memory->create(phase_support, n + 1, n + 1, "pair:phase_support");
}

/* ----------------------------------------------------------------------
 global settings
 ------------------------------------------------------------------------- */

void PairSPHConcASurfaceReactionMultiPhase::settings(int narg, char **arg) {
  if (narg != 0)
    error->all(FLERR,
        "Illegal number of setting arguments for pair_style sph/concprecipitation");
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairSPHConcASurfaceReactionMultiPhase::coeff(int narg, char **arg) {
  if (narg != 4)
    error->all(FLERR,"Incorrect number of args for pair_style sph/concprecipitation coefficients");
  if (!allocated)
    allocate();

  int ilo, ihi, jlo, jhi;
  force->bounds(FLERR,arg[0], atom->ntypes, ilo, ihi);
  force->bounds(FLERR,arg[1], atom->ntypes, jlo, jhi);
  
  double kernel_one = force->numeric(FLERR,arg[2]);
  double phase_one = force->numeric(FLERR,arg[3]);
  
  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut[i][j] = kernel_one;
      phase_support[i][j] = phase_one;
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

double PairSPHConcASurfaceReactionMultiPhase::init_one(int i, int j) {

  if (setflag[i][j] == 0) {
    error->all(FLERR,"All pair sph/concAsurfacereaction/multiphase coeffs are not set");
  }

  cut[j][i] = cut[i][j];
  phase_support[j][i] = phase_support[i][j];
  
  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

double PairSPHConcASurfaceReactionMultiPhase::single(int i, int j, int itype, int jtype,
    double rsq, double factor_coul, double factor_lj, double &fforce) {
  fforce = 0.0;

  return 0.0;
}

/* ---------------------------------------------------------------------- */

int PairSPHConcASurfaceReactionMultiPhase::pack_forward_comm(int n, int *list, double *buf,
					     int pbc_flag, int *pbc)
{
  int i, j, m;
  
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = cA[j];
    buf[m++] = mA[j];
    buf[m++] = atom->type[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairSPHConcASurfaceReactionMultiPhase::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;
  
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    cA[i] = buf[m++];
    mA[i] = buf[m++];
    atom->type[i] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int PairSPHConcASurfaceReactionMultiPhase::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;
  int *type = atom->type;
  
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = dcA[i];
    buf[m++] = dmA[i];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairSPHConcASurfaceReactionMultiPhase::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;
  int *type = atom->type;
  
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];

    dcA[j] += buf[m++];
    dmA[j] += buf[m++];
  }
}
