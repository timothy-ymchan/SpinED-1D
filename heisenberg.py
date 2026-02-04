from build_hamiltonian import build_nn_hamiltonian_nosym, build_nn_hamiltonian_charged_sector
from basis_construction import check_momentum, check_Zn_charge, get_sector_basis, construct_period_table_cache, refine_cache_period_table_by_charge, get_Zn_charge,get_cache_charged_momentum_state,clear_cache
from modular_arithmetic import to_base
from spin_ops import Id, Sx, Sy, Sz
from numpy import kron
import numpy as np
from scipy.linalg import eigh, expm, eig
from scipy.sparse.linalg import eigsh 
import json
import time

def main():
    for N in [3,4]:
        ED_nosym(N)




def Heisenberg_Hamiltonian():
    H_Heis = kron(Sx,Sx) + kron(Sy,Sy) + kron(Sz,Sz)
    return H_Heis

def ED_nosym(nsites):
    Hnn = Heisenberg_Hamiltonian()
    H_nosym = build_nn_hamiltonian_nosym(Hnn,nsites)
    hermitian_check = np.sum((H_nosym - H_nosym.conj().T)**2)
    print(hermitian_check)
    eigval = eigh(H_nosym,eigvals_only=True)
    eigval /= nsites
    result = {'eigval':eigval.tolist(),'nsites':nsites}
    
    # Dump results 
    timestamp = int(time.time())
    with open(f"./results/heisenberg/ED-heisenberg_nsites-{nsites}_{timestamp}.json","w") as file:
        json.dump(result,file)

main()