from build_hamiltonian import build_nn_hamiltonian_nosym, build_nn_hamiltonian_charged_sector
from basis_construction import check_momentum, check_Zn_charge, get_sector_basis,pprint_charge_str
from modular_arithmetic import to_base
from numpy import kron
import numpy as np
from scipy.linalg import eigh, expm, eig
from scipy.sparse.linalg import eigsh 
import time
import json 

X = 0.5*np.array([[0,1],[1,0]])
Z = 0.5*np.array([[1,0],[0,-1]])
Id = np.eye(2)

def main():
    for nsites in [2,4,6]:
        gs = np.linspace(0,1,51)
        for g in gs:
            ED_momentum_Z2(g,nsites)

def TFI_Hamiltonian(g):
    return kron(X,X) - g*(kron(Z,Id))

def ED_nosym(g,nsites):
    Hnn = TFI_Hamiltonian(g)
    H_nosym = build_nn_hamiltonian_nosym(Hnn,nsites)
    print('Block size ', H_nosym.shape)
    hermitian_check = np.sum((H_nosym - H_nosym.conj().T)**2)
    eigval = eigh(H_nosym,eigvals_only=True)
    eigval /= nsites
    result = {'eigval':eigval.tolist(),'g':g,'nsites':nsites}
    
    # Dump results 
    timestamp = int(time.time())
    with open(f"./TFI/nosym/ED_g-{g:.3f}_nsites-{nsites}_{timestamp}.json","w") as file:
        json.dump(result, file)

def ED_momentum(g,nsites):
    Hnn = TFI_Hamiltonian(g)
    result = {'eigval':{},'g':g,'nsites':nsites}
    for k in range(0,nsites):
        k_checker = lambda state, momentum: check_momentum(state,momentum,nsites,base=2) # Build basis on the fly to use less memory 
        sector_basis = get_sector_basis(nsites,(k,),(k_checker,),base=2)

        H = build_nn_hamiltonian_charged_sector(Hnn,nsites,sector_basis,k,base=2,charge=-1,charge_checker=lambda s,c:True) # Build sector 
        print(f'Sector size for {k}: {H.shape}')
        #eig_vals = eigh(H,eigvals_only=True)
        eig_vals = eigh(H,eigvals_only=True)
        eig_vals = np.sort(eig_vals)
        result['eigval'][k] = sorted((eig_vals/nsites).tolist())

    # Dump results 
    timestamp = int(time.time())
    with open(f"./TFI/momentum/ED_g-{g:.3f}_nsites-{nsites}_{timestamp}.json","w") as file:
        json.dump(result,file)

def ED_momentum_Z2(g,nsites):
    Hnn = TFI_Hamiltonian(g)
    result = {'eigval':{},'g':g,'nsites':nsites}
    charge_table = [0,1]
    Z2_checker = lambda state, charge : check_Zn_charge(state,charge,charge_table,mod_n=2,nsites=nsites,base=2) # Check Z3 charge
    k_checker = lambda state, momentum: check_momentum(state,momentum,nsites,base=2)
    result = {'eigval':{},'g':g,'nsites':nsites}

    for k in range(0,nsites):
        for charge in charge_table:
            print(f'k={k} c={charge}')
            sector_basis = get_sector_basis(nsites,(k,charge),(k_checker,Z2_checker),base=2)
            # for b in sector_basis: # To check the charges are right
            #     pprint_charge_str(b,charge_table,nsites)
            # print('Length: ',len(sector_basis))
  
            H = build_nn_hamiltonian_charged_sector(Hnn,nsites,sector_basis,k,charge=charge,charge_checker=Z2_checker,base=2) # Build sector 
            print(f'Sector size for k={k}, charge={charge}: {H.shape}')
            eig_vals = eigh(H,eigvals_only=True)
            result['eigval'][str((k,charge))]=sorted((eig_vals/nsites).tolist())

    # Dump results 
    timestamp = int(time.time())
    with open(f"./TFI/momentumZ2/ED_g-{g:.3f}_nsites-{nsites}_{timestamp}.json","w") as file:
       json.dump(result,file)
if __name__ == "__main__":
    main()