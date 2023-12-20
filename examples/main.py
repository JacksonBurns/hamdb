from datetime import datetime

start = datetime.now()

import numpy as np
from funcs import get_molecule, run_all

# used for geometry optimization
ORGANIC_BASIS_SET = "def2svpjkfit" # "6-31g"
METAL_BASIS_SET = "def2svpjkfit" # "weigend"

monomer = get_molecule("xyzs/monomer.xyz", ORGANIC_BASIS_SET, METAL_BASIS_SET, 0, 0)
cadmium = get_molecule("xyzs/cadmium.xyz", ORGANIC_BASIS_SET, METAL_BASIS_SET, 2, 0)
vdw_complex = get_molecule("xyzs/vdw_sm_planar.xyz", ORGANIC_BASIS_SET, METAL_BASIS_SET, 2, 0)
n_product = get_molecule("xyzs/n_product.xyz", ORGANIC_BASIS_SET, METAL_BASIS_SET, 2, 0)
n_product_ts = get_molecule("xyzs/n_product_ts.xyz", ORGANIC_BASIS_SET, METAL_BASIS_SET, 2, 0)
s_product = get_molecule("xyzs/s_product.xyz", ORGANIC_BASIS_SET, METAL_BASIS_SET, 2, 0)
s_product_ts = get_molecule("xyzs/s_product_ts.xyz", ORGANIC_BASIS_SET, METAL_BASIS_SET, 2, 0)

data_dict = {
    "monomer": {"pyscf_molecule": monomer, "is_ts": False},
    "cadmium": {"pyscf_molecule": cadmium, "is_ts": False},
    "vdw_complex": {"pyscf_molecule": vdw_complex, "is_ts": False},
    "n_product": {"pyscf_molecule": n_product, "is_ts": False},
    "n_product_ts": {"pyscf_molecule": n_product_ts, "is_ts": True},
    "s_product": {"pyscf_molecule": s_product, "is_ts": False},
    "s_product_ts": {"pyscf_molecule": s_product_ts, "is_ts": True},
}

results = run_all(data_dict, logpath="logfiles_mixed/", use_cosmo=True)


def get_G(results, name):
    try:
        return results[name]["therminfo"]["G_tot"][0]
    except:
        return 0.0


def get_E0(results, name):
    try:
        return results[name]["therminfo"]["E0"][0]
    except:
        return 0.0


def get_rates(results, product_name, ts_name):
    # constants
    T, P = 298.15, 101325
    R = 1.987 * 10**-3
    R_J = 8.3145
    ha_to_kcal = 627.509
    kb_over_h = 2.083661912 * 10**10

    del_G_overall = (get_G(results, product_name) - get_G(results, "vdw_complex")) * ha_to_kcal
    del_G_f = (get_G(results, ts_name) - get_G(results, "vdw_complex")) * ha_to_kcal
    del_G_r = (get_G(results, ts_name) - get_G(results, product_name)) * ha_to_kcal
    K_eq = np.exp(-del_G_overall / (R * T))
    K_c = K_eq * R_J * T / P
    k_f = kb_over_h * T * np.exp(-del_G_f / (R * T)) * (R_J * T / P)
    k_r = kb_over_h * T * np.exp(-del_G_r / (R * T)) * (R_J * T / P)
    del_E_overall = (get_E0(results, product_name) - get_E0(results, "vdw_complex")) * ha_to_kcal
    del_E_f = (get_E0(results, ts_name) - get_E0(results, "vdw_complex")) * ha_to_kcal
    del_E_r = (get_E0(results, ts_name) - get_E0(results, product_name)) * ha_to_kcal
    print(
        """
Results for {:s} from vdw:
 K_eq                  = {:.4e} 1
 K_c                   = {:.4e} 1
 Overall Energy Change = {: .4f}    kcal
 k_forward             = {:.4e} 1/s
 Forward Barrier       = {: .4f}    kcal
 k_reverse             = {:.4e} C/s
 Reverse Barrier       = {: .4f}    kcal
 k_forward/k_reverse   = {:.4e} 1
    """.format(
            product_name,
            K_eq,
            K_c,
            del_E_overall,
            k_f,
            del_E_f,
            k_r,
            del_E_r,
            k_f / k_r,
        )
    )

    del_G_overall = (get_G(results, product_name) - (get_G(results, "monomer") + get_G(results, "cadmium"))) * ha_to_kcal
    del_G_f = (get_G(results, ts_name) - (get_G(results, "monomer") + get_G(results, "cadmium"))) * ha_to_kcal
    del_G_r = (get_G(results, ts_name) - get_G(results, product_name)) * ha_to_kcal
    K_eq = np.exp(-del_G_overall / (R * T))
    K_c = K_eq * R_J * T / P
    k_f = kb_over_h * T * np.exp(-del_G_f / (R * T)) * (R_J * T / P)
    k_r = kb_over_h * T * np.exp(-del_G_r / (R * T)) * (R_J * T / P)
    del_E_overall = (get_E0(results, product_name) - (get_E0(results, "monomer") + get_E0(results, "cadmium"))) * ha_to_kcal
    del_E_f = (get_E0(results, ts_name) - (get_E0(results, "monomer") + get_E0(results, "cadmium"))) * ha_to_kcal
    del_E_r = (get_E0(results, ts_name) - get_E0(results, product_name)) * ha_to_kcal
    print(
        """
Results for {:s} from free starting materials:
 K_eq                  = {:.4e} 1
 K_c                   = {:.4e} 1
 Overall Energy Change = {: .4f}    kcal
 k_forward             = {:.4e} 1/s
 Forward Barrier       = {: .4f}    kcal
 k_reverse             = {:.4e} C/s
 Reverse Barrier       = {: .4f}    kcal
 k_forward/k_reverse   = {:.4e} 1
    """.format(
            product_name,
            K_eq,
            K_c,
            del_E_overall,
            k_f,
            del_E_f,
            k_r,
            del_E_r,
            k_f / k_r,
        )
    )

    del_G_overall = (get_G(results, "vdw_complex") - (get_G(results, "monomer") + get_G(results, "cadmium"))) * ha_to_kcal
    K_eq = np.exp(-del_G_overall / (R * T))
    K_c = K_eq
    del_E_overall = (get_G(results, "vdw_complex") - (get_E0(results, "monomer") + get_E0(results, "cadmium"))) * ha_to_kcal
    print(
        """
Results for vdW complexation:
 K_eq                  = {:.4e} 1
 K_c                   = {:.4e} 1
 Overall Energy Change = {: .4f}    kcal
    """.format(
            K_eq, K_c, del_E_overall
        )
    )


get_rates(results, "s_product", "s_product_ts")
get_rates(results, "n_product", "n_product_ts")

end = datetime.now()
print("\n", "   REST IN PEACE BANNEDBYGAUSSIAN.ORG AND LONG LIVE FOSS.")
print("        - decent people the world over", "\n")
print("TOTAL EXECUTION TIME: ", end - start)
