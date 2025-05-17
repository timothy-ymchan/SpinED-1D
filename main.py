from build_hamiltonian import build_nn_hamiltonian_nosym, build_nn_hamiltonian_charged_sector
from basis_construction import check_momentum, check_Zn_charge, get_sector_basis,pprint_charge_str
from modular_arithmetic import to_base
from spin_ops import Id, Sx, Sy, Sz
from numpy import kron
import numpy as np
from scipy.linalg import eigh, expm, eig
from scipy.sparse.linalg import eigsh 
import json
import time

def main():
    lamb = np.sqrt(3)/2
    mu = 2
    ED_momentum_Z3(lamb,mu,10)



def A4_Hamiltonian(lamb,mu):
    H_Heis = kron(Sx,Sx) + kron(Sy,Sy) + kron(Sz,Sz)
    H_q = kron(Sx,Sx)@kron(Sx,Sx) + kron(Sy,Sy)@kron(Sy,Sy) + kron(Sz,Sz)@kron(Sz,Sz)
    H_c = kron(Sx@Sy,Sz) + kron(Sz@Sx,Sy) + kron(Sy@Sz,Sx) + kron(Sy@Sx,Sz) + kron(Sx@Sz,Sy) + kron(Sz@Sy,Sx) + kron(Sx,Sy@Sz) + kron(Sz,Sx@Sy) + kron(Sy,Sz@Sx) + kron(Sx,Sz@Sy) + kron(Sz,Sy@Sx) + kron(Sy,Sx@Sz)
    
    return H_Heis + lamb*H_c + mu*H_q

def ED_nosym(lamb,mu,nsites):
    Hnn = A4_Hamiltonian(lamb,mu)
    H_nosym = build_nn_hamiltonian_nosym(Hnn,nsites)
    hermitian_check = np.sum((H_nosym - H_nosym.conj().T)**2)
    print(hermitian_check)
    eigval = eigh(H_nosym,eigvals_only=True)
    eigval /= nsites
    result = {'eigval':eigval.tolist(),'lamb':lamb,'mu':mu,'nsites':nsites}
    
    # Dump results 
    timestamp = int(time.time())
    with open(f"./results/nosym/ED_lamb-{lamb}_mu-{mu}_nsites-{nsites}_{timestamp}.json","w") as file:
        json.dump(result,file)

def ED_momentum(lamb,mu,nsites):
    Hnn = A4_Hamiltonian(lamb,mu)
    
    result = {'eigval':{},'lamb':lamb,'mu':mu,'nsites':nsites}
    for k in range(0,nsites):
        k_checker = lambda state, momentum: check_momentum(state,momentum,nsites) # Build basis on the fly to use less memory 
        sector_basis = get_sector_basis(nsites,(k),(k_checker,))

        H = build_nn_hamiltonian_charged_sector(Hnn,nsites,sector_basis,k) # Build sector 
        print(f'Sector size for {k}: {H.shape}')
        #eig_vals = eigh(H,eigvals_only=True)
        eig_vals = eigsh(H,k=20,which='LM',sigma=0,return_eigenvectors=False)
        eig_vals = np.sort(eig_vals)
        result['eigval'][k] = sorted((eig_vals/nsites).tolist())

    # Dump results 
    timestamp = int(time.time())
    with open(f"./results/momentum/ED_lamb-{lamb}_mu-{mu}_nsites-{nsites}_{timestamp}.json","w") as file:
        json.dump(result,file)

def ED_momentum_Z3(lamb,mu,nsites):
    # Get charged basis 
    U = expm((-1.j*2*np.pi/3)*(Sx+Sy+Sz)/np.sqrt(3))
    _, V = eig(U) # D = V^dagger U V ; Charge -1,0,1 (index: 0,1,2)
    charge_table = [-1,0,1]

    # Rotate Hamiltonian to charged basis
    Hnn = A4_Hamiltonian(lamb=lamb,mu=mu)
    S = kron(V,V)
    Hnn = S.conj().T@Hnn@S
    
    Z3_checker = lambda state, charge : check_Zn_charge(state,charge,charge_table,mod_n=3,nsites=nsites) # Check Z3 charge
    k_checker = lambda state, momentum: check_momentum(state,momentum,nsites)
    result = {'eigval':{},'lamb':lamb,'mu':mu,'nsites':nsites}
    for k in range(0,nsites):
        for charge in charge_table:
            print(f'k={k} c={charge}')
            sector_basis = get_sector_basis(nsites,(k,charge),(k_checker,Z3_checker))
            # for b in sector_basis: # To check the charges are right
            #     pprint_charge_str(b,charge_table,nsites)
            # print('Length: ',len(sector_basis))
  
            H = build_nn_hamiltonian_charged_sector(Hnn,nsites,sector_basis,k,charge=charge,charge_checker=Z3_checker) # Build sector 
            print(f'Sector size for k={k}, charge={charge}: {H.shape}')
            eig_vals = eigh(H,eigvals_only=True)
            result['eigval'][str((k,charge))]=sorted((eig_vals/nsites).tolist())
    # Dump results 
    timestamp = int(time.time())
    with open(f"./results/momentumZ3/ED_lamb-{lamb}_mu-{mu}_nsites-{nsites}_{timestamp}.json","w") as file:
       json.dump(result,file)

if __name__ == "__main__":
    main()
    

