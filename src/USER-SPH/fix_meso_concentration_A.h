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

FixStyle(meso/concentrationA,FixMesoConcentrationA)

#else

#ifndef LMP_FIX_MESO_CONCENTRATION_A_H
#define LMP_FIX_MESO_CONCENTRATION_A_H

#include "fix.h"

namespace LAMMPS_NS {

class FixMesoConcentrationA : public Fix {
 public:
  FixMesoConcentrationA(class LAMMPS *, int, char **);
  virtual ~FixMesoConcentrationA();
  int setmask();
  virtual void init();
  virtual void final_integrate();
  // virtual void end_of_step();
  void reset_dt();

  // int pack_forward_comm(int, int *, double *, int, int *);
  // void unpack_forward_comm(int, int, double *);

 private:
  class NeighList *list;
 protected:
  double dtcA;
  double *step_respa;
  double *cA, *cAeq, *dcA, *mA, *dmA, *RA, *mAthres, *rmass;
  int mass_require;

  class Pair *pair;
};

}

#endif
#endif
