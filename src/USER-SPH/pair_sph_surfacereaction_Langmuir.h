#ifdef PAIR_CLASS

PairStyle(sph/surfacereaction/Langmuir, PairSPHSurfaceReactionLangmuir)

#else

#ifndef LMP_PAIR_SPH_SURFACEREACTION_LANGMUIR_H
#define LMP_PAIR_SPH_SURFACEREACTION_LANGMUIR_H

#include "pair.h"

namespace LAMMPS_NS {
  class PairSPHSurfaceReactionLangmuir : public Pair {
  public:
    PairSPHSurfaceReactionLangmuir(class LAMMPS *);
    virtual ~PairSPHSurfaceReactionLangmuir();
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
    // Normal vector for the surface reaction
    double *nx, *ny, *nz;
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
    void allocate();
  };
}

#endif
#endif
