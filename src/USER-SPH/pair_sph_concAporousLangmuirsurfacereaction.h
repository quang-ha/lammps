#ifdef PAIR_CLASS

PairStyle(sph/concAporousLangmuirsurfacereaction, PairSPHConcAPorousLangmuirSurfaceReaction)

#else

#ifndef LMP_PAIR_SPH_POROUS_LANGMUIR_SURFACEREACTION_H
#define LMP_PAIR_SPH_POROUS_LANGMUIR_SURFACEREACTION_H

#include "pair.h"

namespace LAMMPS_NS {
  class PairSPHConcAPorousLangmuirSurfaceReaction : public Pair {
  public:
    PairSPHConcAPorousLangmuirSurfaceReaction(class LAMMPS *);
    virtual ~PairSPHConcAPorousLangmuirSurfaceReaction();
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
    double **cut;
    // Mass fraction for aqueous species and absrobed species
    double *xA, *dxA, *yA, *dyA;
    // Normalised mass fraction
    double *thetaA;
    // Binary diffusion coefficient
    double *DA;
    // Maximum absorbed concentration
    double yAmax;
    // Adsorption and desorption rate coefficient
    double kaA, kdA;
    // Number of adsorption sites needed
    int lambda;
    // Check if periodic property is on
    int is_periodic;
    // Surface area and micro-pore volume for porous solid
    double As, Vp;
    void allocate();
  };
}

#endif
#endif
