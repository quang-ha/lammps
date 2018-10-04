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

FixStyle(meso/concALangmuir,FixMesoConcALangmuir)

#else

#ifndef LMP_FIX_MESO_CONCA_LANGMUIR_H
#define LMP_FIX_MESO_CONCA_LANGMUIR_H

#include "fix.h"

namespace LAMMPS_NS {

class FixMesoConcALangmuir : public Fix {
 public:
  FixMesoConcALangmuir(class LAMMPS *, int, char **);
  int setmask();
  virtual void init();
  virtual void init_list(int, class NeighList *);
  virtual void initial_integrate(int);
  virtual void final_integrate();
  virtual void end_of_step();
  void reset_dt();

 private:
  class NeighList *list;

 protected:
  double dtv,dtf;
  double h;
  double *step_respa;
  double *cA, *dcA, *yA, *dyA, *yAmax, *DA, *kAa, *kAd, *thetaA, *sA;
  int mass_require;

  class Pair *pair;
};

}

#endif
#endif
