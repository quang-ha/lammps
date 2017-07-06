For the general installation instruction please refer to the
documentation of LAMMPS and to SPH-USER package:

- http://lammps.sandia.gov/doc/Section_start.html
- http://lammps.sandia.gov/doc/USER/sph/SPH_LAMMPS_userguide.pdf

The original code is from https://github.com/slitvinov/lammps-sph-multiphase where the remaining of the README is based on.

* Implementation

We add the following extension to USER-SPH package:

** atom_style meso/multiphase:

This is data structures which provides
- position
- velocity
- extrapolated velocity (`vest`)
- forces
- SPH density (`rho`)
- time derivative of SPH density (`drho`)
- internal energy per particle (`e`)
- time derivative of internal energy per particle (`de`)
- color gradient vector (`colorgradient`)
- per-particle heat capacity (`cv`)

This data structure can be activated by
```
atom_style meso/multiphase
```

** pair_sph_colorgradient

A [[http://lammps.sandia.gov/doc/pair_style.html][pair_style]] to calculate a color gradient

```
pair_style         sph/colorgradient
pair_coeff         I J     ${h} ${alpha}
```

Here, `I` and `J` are the types of SPH particles for which a color gradient is calculated, `alpha` is a surface tension coefficient, `h` is a cutoff.

** pair_sph_surfacetension

A [[http://lammps.sandia.gov/doc/pair_style.html][pair_style]] to calculate surface tension

```
pair_coeff         I J     sph/surfacetension ${h}
```

Here, `I` and `J` are the types of SPH particles for which a surface tension is calculated, `h` is a cutoff. Note that surface tension coefficient is included into color gradient.

** pair_sph_heatconduction_phasechange

A modified heat conduction equation to use for phase change model. Has to forms. Simple form is equivalent to the heat conduction equation from USER-SPH package.

```
pair_coeff         I J  sph/heatconduction/phasechange  ${D_heat_ld}
```

Here, `I` and `J` are the types of SPH particles which interact and
`D` is a heat diffusion coefficient.

Full form of the pair style is

```
pair_coeff         I J  sph/heatconduction/phasechange  ${D_heat_ld} TI TJ
```
where `TI` and `TJ` are temperatures for corresponding particles in
`I` and `J` interactions.

`NULL` can be used as a placeholder to indicate that normal temperate
should be used for corresponding particle

```
pair_coeff         I J  sph/heatconduction/phasechange  ${D_heat_ld} TI NULL
```

** fix_phase_change

Fix which adds a phase change

```
fix                fix_ID group_ID phase_change &
                   ${Tc} ${Tt} ${Hwv} ${dr} ${mass_v} &
		                      ${pcutoff} ${l_type} ${v_type} ${insert_every} 123456 ${prob} region
```

`fix_ID` and `group_ID` are described in LAMMPS documentation. `TC` is critical temperature of the phase change, =TT= is transition temperature for the algorithm (should be set above `TC`), `dr` a characteristic distance for a new particle position, `mass` a mass of a new particle, =h= cutoff of the interaction, `from_type` and `to_type` types of the particles involved in phase transition, `N` frequency of the check for phase transition algorithm, `seed` a seed for random number generator, `prob` probability of the phase transition if all criteria are met (`0<prob<1`), `region` a region where algorithm checks for potential phase transition.

* Examples
See [[file:examples/USER/sph/]]
