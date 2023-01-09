#!/usr/bin/env python

import ase.calculators.vasp
from ase.calculators.vasp import Vasp


def base_calc() -> ase.calculators.vasp.Vasp:
    return Vasp(kpts=(1, 1, 1),  # specifies the Bloch vectors (k points) used to sample the Brillouin zone
                potim=0.5,       # sets the time step in molecular dynamics or the step width in ionic relaxations
                encut=700,       # specifies the cutoff energy for the plane-wave-basis set in eV
                xc='PBE',        # XC functional
                ivdw=12,         # specifies a vdW dispersion correction, 12 implements D3 Becke-Johnson
                ediff=1E-6,      # specifies the global break condition for the electronic SC-loop, units in eV
                ediffg=-0.02,    # defines the break condition for the ionic relaxation loop
                isif=3,          # sets either relaxation or molecular dynamics runs
                ibrion=2,        # determines how the ions are updated and moved, 2 implements conjugate gradient algo
                nsw=100,         # sets the maximum number of ionic steps
                nelm=60,         # sets the maximum number of electronic self-consistency steps
                nelmin=4,        # specifies the minimum number of electronic self-consistency steps
                nelmdl=-4,       # specifies the number of non-self-consistent steps at the beginning, -ve means only in first ionic step
                ispin=2,         # spin-polarized calculations are performed, needs to be true for all radicals
                nupdown=-1,      # sets the difference between the number of electrons in the up and down spin components, -1 does a full relaxation
                prec='Normal',   # specifies the "precision"-mode
                algo='Fast',     # selects a fairly robust mixture of the Davidson and RMM-DIIS algorithms
                istart=1,        # determines whether or not to read the WAVECAR
                ismear=0,        # determines how the partial occupancies fnk are set for each orbital
                sigma=0.05,      # specifies the width of the smearing in eV
                isym=0,          # determines the way VASP treats symmetry
                lorbit=11,       # PAW projectors are used so that the projection is strictly within the PAW sphere
                nwrite=1,        # verbose level of OUTCAR
                icharg=2,        # determines how VASP constructs the initial charge density
                voskown=1,       # determines whether Vosko-Wilk-Nusair interpolation is used or not
                lplane=True,     # switches on the plane-wise data distribution in real space
                lreal='Auto',    # False for gas phase, determines whether the projection operators are evaluated in real-space or in reciprocal space
                lasph=True,      # include non-spherical contributions related to the gradient of the density in the PAW spheres
                lcharg=False,    # determines whether the charge densities files CHGCAR and CHG are written
                laechg=False,    # reconstruction the all-electron charge density
                lwave=False,     # determines whether the wavefunctions are written to the WAVECAR file at the end of a run
                npar=4,          # determines the number of bands that are treated in parallel, on HPC set it to sqrt of number of cores
                nsim=4,          # sets the number of bands that are optimized simultaneously by the RMM-DIIS algorithm
                )


def isif_2_calc() -> ase.calculators.vasp.Vasp:
    calc = base_calc()
    calc.set(isif=2, nsw=500, nelm=100)
    return calc


def isif_3_calc() -> ase.calculators.vasp.Vasp:
    calc = base_calc()
    calc.set(isif=3, nsw=500, nelm=100)
    return calc


def isif_4_calc() -> ase.calculators.vasp.Vasp:
    calc = base_calc()
    calc.set(isif=4, nsw=500, nelm=100)
    return calc


def nvt_calc() -> ase.calculators.vasp.Vasp:
    calc = base_calc()
    calc.set(encut=520, ibrion=0, isif=2, nsw=2000, potim=0.5, nelm=60)
    calc.set(mdalgo=2, tebeg=400.0, smass=1.0)  # Nose-Hoover implementation
    return calc


def npt_calc() -> ase.calculators.vasp.Vasp:
    calc = base_calc()
    calc.set(encut=520, ibrion=0, isif=3, nsw=2000, potim=0.5, nelm=60)
    calc.set(mdalgo=3, pstress=0.001, tebeg=400.0,  # Langevin implementation
             pmass=1000, langevin_gamma_l=0, langevin_gamma=(0, 0, 0, 0))
    return calc


def spe_calc() -> ase.calculators.vasp.Vasp:
    calc = base_calc()
    calc.set(encut=520, nsw=500, nelm=100)
    return calc


def vtst_neb_calc(my_images) -> ase.calculators.vasp.Vasp:
    calc = base_calc()
    calc.set(encut=400, ediffg=-0.05, isif=2, potim=0.2, nsw=500, nelm=100)
    return calc
