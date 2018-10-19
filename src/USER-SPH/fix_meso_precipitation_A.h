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

FixStyle(meso/precipitationA,FixMesoPrecipitationA)

#else

#ifndef LMP_FIX_MESO_PRECIPITATION_A_H
#define LMP_FIX_MESO_PRECIPITATION_A_H

#include "fix.h"

namespace LAMMPS_NS {

class FixMesoPrecipitationA : public Fix {
 public:
  FixMesoPrecipitationA(class LAMMPS *, int, char **);
  virtual ~FixMesoPrecipitationA();
  int setmask();
  virtual void init();
  virtual void init_list(int, class NeighList *);
  virtual void initial_integrate(int);
  virtual void final_integrate();
  virtual void end_of_step();
  void reset_dt();

  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);

  void unpack_reverse_comm(int, int *, double *);
  int pack_reverse_comm(int, int, double *);
  
 private:
  class NeighList *list;
  
 protected:
  double dtv,dtf;
  double *step_respa;
  double *cA, *cAeq, *dcA, *mA, *dmA, *RA, *mAthres, *rmass;
  double *ischangecA;
  int mass_require;

  class Pair *pair;
};

}

#endif
#endif
