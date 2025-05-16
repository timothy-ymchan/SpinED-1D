from build_hamiltonian import build_nn_hamiltonian_nosym, build_nn_hamiltonian_sector
from basis_construction import construct_momentum_table, construct_period_table
from spin_ops import Id, Sx, Sy, Sz
from numpy import kron
import numpy as np
from scipy.linalg import eigh
import json
import time

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
    pt = construct_period_table(nsites)
    mt = construct_momentum_table(pt)
    result = {'eigval':{},'lamb':lamb,'mu':mu,'nsites':nsites}
    for k in mt.keys():
        H = build_nn_hamiltonian_sector(Hnn,nsites,k,mt,pt)
        result['eigval'][k] = (eigh(H,eigvals_only=True)/nsites).tolist()

    # Dump results 
    timestamp = int(time.time())
    with open(f"./results/momentum/ED_lamb-{lamb}_mu-{mu}_nsites-{nsites}_{timestamp}.json","w") as file:
        json.dump(result,file)


if __name__ == "__main__":
    lambc = np.sqrt(3)/2
    mu = 2
    #for nsites in []:
    #nsites = 10
    #ED_momentum(lambc,mu,nsites)
    


    

