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

FixStyle(phase_change_threshold,FixPhaseChangeThreshold)

#ifndef LMP_FIX_PHASECHANGE_THRESHOLD_H
#define LMP_FIX_PHASECHANGE_THRESHOLD_H

#include <stdio.h>
#include "fix.h"

namespace LAMMPS_NS{

class FixPhaseChangeThreshold : public Fix {
 public:
  FixPhaseChangeThreshold(class LAMMPS*, int, char **);
  ~FixPhaseChangeThreshold();
  int setmask();
  void init();
  void init_list(int, class NeighList *);
  void pre_exchange();
  void write_restart(FILE *);
  void restart(char *);

 private:
  int from_type, to_type, nfreq;
  int iregion;
  char *idregion;
  double xlo, xhi, ylo, yhi, zlo, zhi;
  double cutoff;
  int nfirst;
  class NeighList *list;

  // threshold temperature for changing
  double Th;

  
  } 
}

#endif
#endif
