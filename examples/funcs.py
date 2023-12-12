import json
import os
from datetime import datetime
import warnings

warnings.filterwarnings("ignore", message="Since PySCF.*")
                        
from pyscf import gto, solvent
from gpu4pyscf.dft import RKS
from pyscf.hessian import thermo

os.environ["OMP_NUM_THREADS"] = "96"  # number of threads
os.environ["PYSCF_MAX_MEMORY"] = "192000"  # MB of memory

FUNCTIONAL = "B3LYP"


def get_molecule(fname, basis_set, charge, spin):
    mol = gto.M(
        atom=fname,
        basis=basis_set,
    )
    mol.charge = charge
    mol.spin = spin
    return mol


def optfreq(mol, is_ts=False, use_cosmo=False, do_geomopt=True):
    # higher level of theory for single point energy
    if not do_geomopt or mol.natm == 1:
        mol.basis = "ano-rcc"
        mol.build()
    # returns a QM method object
    mean_field = RKS(
        mol,
        xc=FUNCTIONAL,
    )
    # inject COSMO into method object
    if use_cosmo:
        # mean_field = mean_field.ddCOSMO()
        mean_field = mean_field.PCM()
        # mf.with_solvent.lebedev_order = 29 # 302 Lebedev grids
        mean_field.with_solvent.method = 'COSMO'
        mean_field.with_solvent.eps = 78.3553
    # skip geometry optimization for monatomics
    if do_geomopt and mol.natm > 1:
        if use_cosmo:
            # https://github.com/pyscf/pyscf/blob/14d88828cd1f18f1e5358da1445355bde55322a1/examples/solvent/21-tddft_geomopt.py
            mean_field.with_solvent.equilibrium_solvation = True
        optimizer_parameters = dict(
            prefix=mol.output,
            assert_convergence=True,
        )
        if is_ts:
            optimizer_parameters["transition"] = True
        mol_eq = mean_field.Gradients().optimizer(solver="geomeTRIC").kernel(optimizer_parameters) # returns a pyscf molecule
        # solve again with a larger basis set for a more accurate single point energy
        return optfreq(mol_eq, is_ts, use_cosmo, do_geomopt=False)
    else:
        mean_field.kernel()
    return mean_field


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
        print("Simulating:", name, flush=True)
        data["pyscf_molecule"].output = logpath + name + "-" + datetime.now().strftime("%Y%m%d-%H%M%S") + ".log"
        data["pyscf_molecule"].verbose = 4
        data["pyscf_molecule"].build()  # make the above logging changes take affect
        mol_optim = optfreq(data["pyscf_molecule"], data["is_ts"], use_cosmo=use_cosmo)
        data["mol_optim"] = mol_optim
        thermoinfo = therm(mol_optim)
        data["therminfo"] = thermoinfo
        print(json.dumps(data, indent=4, default=lambda _: "<not serializble>"))

    print("RESULTS START")
    print(json.dumps(data_dict, indent=4, default=lambda _: "<not serializble>"))
    print("RESULTS END", flush=True)
    return data_dict


# testing code: python funcs.py
if __name__ == "__main__":
    water = get_molecule("xyzs/water.xyz", "sto-3g", 0, 0)

    data_dict = {
        "water": {"pyscf_molecule": water, "is_ts": False},
    }

    results = run_all(data_dict, use_cosmo=True)
