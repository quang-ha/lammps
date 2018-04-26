#ifdef PAIR_CLASS

PairStyle(sph/concADarcyreaction/multiphase, PairSPHConcADarcyReactionMultiPhase)

#else

#ifndef LMP_PAIR_SPH_DARCYREACTION_MULTIPHASE_H
#define LMP_PAIR_SPH_DARCYREACTION_MULTIPHASE_H

#include "pair.h"

namespace LAMMPS_NS {
  class PairSPHConcADarcyReactionMultiPhase : public Pair {
  public:
    PairSPHConcADarcyReactionMultiPhase(class LAMMPS *);
    virtual ~PairSPHConcADarcyReactionMultiPhase();
    void init_style();
    virtual void compute(int, int);
    void settings(int, char **);
    void coeff(int, char **);
    virtual double init_one(int, int);
    virtual double single(int, int, int, int, double, double, double, double &);

    int pack_forward_comm(int, int *, double *, int, int *);
    void unpack_forward_comm(int, int, double *);
    
    int pack_reverse_comm(int, int, double *);
    void unpack_reverse_comm(int, int *, double *);
    
  protected:
    double **cut, **phasecut;
    double *cA, *dcA, *DAeff, *RA, *cAeq, *kAeff;
    void allocate();
  };
}

#endif
#endif
