import json
import os
from datetime import datetime

from pyscf import gto, solvent
from pyscf.geomopt.geometric_solver import optimize
from pyscf.hessian import thermo
from pyscf.qsdopt.qsd_optimizer import QSD

os.environ["OMP_NUM_THREADS"] = "96"  # number of threads
os.environ["PYSCF_MAX_MEMORY"] = "192000"  # MB of memory

FUNCTIONAL = "B3LYP"


def get_molecule(fname, basis_set):
    return gto.M(
        atom=fname,
        basis=basis_set,
    )


def optfreq(mol, is_ts=False, use_cosmo=False):
    if mol.natm > 1:  # skip geometry optimization for monatomics
        mean_field = mol.KS(
            mol,
            xc=FUNCTIONAL,
        )
        if use_cosmo:
            mean_field = mean_field.ddCOSMO()
            # https://github.com/pyscf/pyscf/blob/14d88828cd1f18f1e5358da1445355bde55322a1/examples/solvent/21-tddft_geomopt.py
            mean_field.with_solvent.equilibrium_solvation = True

        optimizer_parameters = dict(
            prefix=mol.output,
            assert_convergence=True,
        )
        if is_ts:
            optimizer_parameters["transition"] = True
        mol = mean_field.Gradients().optimizer(solver="geomeTRIC").kernel(optimizer_parameters)

    mol_eq = mol

    print("~" * 50, "\n", " " * 15 + "OPTIMIZED GEOMETRY START", "\n", "~" * 50, flush=True)
    print(mol_eq.atom_coords(unit="ANG"))
    print("~" * 50, "\n", " " * 15 + " OPTIMIZED GEOMETRY END", "\n", "~" * 50, flush=True)

    # solve wavefunction - re-initalize with new geometry
    mean_field_eq = mol_eq.KS(
        mol_eq,
        xc=FUNCTIONAL,
    )
    if use_cosmo:
        mean_field_eq = mean_field_eq.ddCOSMO()
    mean_field_eq.kernel(dm0=None)
    return mean_field_eq


def therm(mean_field_eq, temperature_K=298.15, pressure_Pa=101325):
    hessian = mean_field_eq.Hessian().kernel()
    freq_info = thermo.harmonic_analysis(mean_field_eq.mol, hessian)
    thermo_info = thermo.thermo(
        mean_field_eq,
        freq_info["freq_au"],
        temperature_K,
        pressure_Pa,
    )
    return thermo_info


def run_all(data_dict, logpath="logfiles/", use_cosmo=False):
    for name, data in data_dict.items():
        try:
            print("~" * 50, "\n", " " * 15 + "Simulating:", name, "\n", "~" * 50, flush=True)
            data["pyscf_molecule"].output = logpath + name + "-" + datetime.now().strftime("%Y%m%d-%H%M%S") + ".log"
            data["pyscf_molecule"].verbose = 4
            data["pyscf_molecule"].build()  # make the above logging changes take affect
            mean_field_eq = optfreq(data["pyscf_molecule"], data["is_ts"], use_cosmo=use_cosmo)
            data["mean_field_eq"] = mean_field_eq
            thermoinfo = therm(mean_field_eq)
            data["therminfo"] = thermoinfo
            print(json.dumps(data, indent=4, default=lambda _: "<not serializble>"))
        except Exception as e:
            print(e)

    print("~" * 50, "\n", " " * 15 + "RESULTS START", "\n", "~" * 50, flush=True)
    print(json.dumps(data_dict, indent=4, default=lambda _: "<not serializble>"))
    print("~" * 50, "\n", " " * 15 + " RESULTS END", "\n", "~" * 50, flush=True)
    return data_dict
