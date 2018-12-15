/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(sph/surfacereaction/Langmuir/porous/simple/2, FixSPHSurfaceReactionLangmuirPorousSimple2)

#else

#ifndef LMP_FIX_SPH_SURFACEREACTION_LANGMUIR_POROUS_SIMPLE_2_H
#define LMP_FIX_SPH_SURFACEREACTION_LANGMUIR_POROUS_SIMPLE_2_H

#include "fix.h"

namespace LAMMPS_NS {

  class FixSPHSurfaceReactionLangmuirPorousSimple2 : public Fix {
  public:
    FixSPHSurfaceReactionLangmuirPorousSimple2(class LAMMPS *, int, char **);
    int setmask();
    virtual void init();
    virtual void final_integrate();
    void reset_dt();

  private:
    class NeighList *list;

  protected:
    // Time step for update
    double dtxA;
    // Change in aqueous species
    double *xA, *dxA;
    // Change in absorbed species
    double *yA, *dyA;

    class Pair *pair;
  };

}

#endif
#endif
