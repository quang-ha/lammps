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

FixStyle(meso/inletconcentrationA,FixMesoInletConcentrationA)

#else

#ifndef LMP_FIX_MESO_INLET_CONCENTRATION_A_H
#define LMP_FIX_MESO_INLET_CONCENTRATION_A_H

#include "fix.h"

namespace LAMMPS_NS {

class FixMesoInletConcentrationA : public Fix {
 public:
  FixMesoInletConcentrationA(class LAMMPS *, int, char **);
  int setmask();
  virtual void init();
  virtual void pre_exchange();
  
 private:
  class NeighList *list;
  // region to keep the inlet concentration
  int iregion, maxattempt, scaleflag;
  char *idregion;
  double xlo, xhi, ylo, yhi, zlo, zhi;
  // concentration to keep the inlet at
  double cAin;

  void options(int, char **);
  
 protected:
  double dtv,dtf;
  double *step_respa;
  double *cA, *dcA;
  int mass_require;

  class Pair *pair;
};

}

#endif
#endif
