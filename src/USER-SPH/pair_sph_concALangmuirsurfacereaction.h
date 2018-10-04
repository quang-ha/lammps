#ifdef PAIR_CLASS

PairStyle(sph/concALangmuirsurfacereaction, PairSPHConcALangmuirSurfaceReaction)

#else

#ifndef LMP_PAIR_SPH_LANGMUIR_SURFACEREACTION_H
#define LMP_PAIR_SPH_LANGMUIR_SURFACEREACTION_H

#include "pair.h"

namespace LAMMPS_NS {
  class PairSPHConcALangmuirSurfaceReaction : public Pair {
  public:
    PairSPHConcALangmuirSurfaceReaction(class LAMMPS *);
    virtual ~PairSPHConcALangmuirSurfaceReaction();
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
    double *cA, *dcA, *yA, *dyA, *yAmax, *DA, *kAa, *kAd, *thetaA, *sA;
    void allocate();
  };
}

#endif
#endif
