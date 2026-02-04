# Construct the 8 site hamiltonian directly
from spin_ops import Sx,Sy,Sz, Id
from numpy import kron 
import numpy as np 
from scipy.linalg import eigh
import time
import json

def A4_Hamiltonian(lamb,mu):
    H_Heis = kron(Sx,Sx) + kron(Sy,Sy) + kron(Sz,Sz)
    H_q = kron(Sx,Sx)@kron(Sx,Sx) + kron(Sy,Sy)@kron(Sy,Sy) + kron(Sz,Sz)@kron(Sz,Sz)
    H_c = kron(Sx@Sy,Sz) + kron(Sz@Sx,Sy) + kron(Sy@Sz,Sx) + kron(Sy@Sx,Sz) + kron(Sx@Sz,Sy) + kron(Sz@Sy,Sx) + kron(Sx,Sy@Sz) + kron(Sz,Sx@Sy) + kron(Sy,Sz@Sx) + kron(Sx,Sz@Sy) + kron(Sz,Sy@Sx) + kron(Sy,Sx@Sz)
    
    return H_Heis + lamb*H_c + mu*H_q

def Hamiltonian_7site(lamb,mu):
    Hnn = A4_Hamiltonian(lamb,mu)

    H = (
        kron(kron(kron(kron(kron(Hnn, Id), Id), Id), Id), Id) + 
        kron(kron(kron(kron(kron(Id, Hnn), Id), Id), Id), Id) +
        kron(kron(kron(kron(kron(Id, Id), Hnn), Id), Id), Id) +
        kron(kron(kron(kron(kron(Id, Id), Id), Hnn), Id), Id) +
        kron(kron(kron(kron(kron(Id, Id), Id), Id), Hnn), Id) +
        kron(kron(kron(kron(kron(Id, Id), Id), Id), Id), Hnn) 
    )

    return H

def Hamiltonian_6site(lamb,mu):
    Hnn = A4_Hamiltonian(lamb,mu)

    H = (
        kron(kron(kron(kron(Hnn, Id), Id), Id), Id) + 
        kron(kron(kron(kron(Id, Hnn), Id), Id), Id) +
        kron(kron(kron(kron(Id, Id), Hnn), Id), Id) +
        kron(kron(kron(kron(Id, Id), Id), Hnn), Id) +
        kron(kron(kron(kron(Id, Id), Id), Id), Hnn)
    )

    return H

def Hamiltonian_5site(lamb,mu):
    Hnn = A4_Hamiltonian(lamb,mu)

    H = (
        kron(kron(kron(Hnn, Id), Id), Id) + 
        kron(kron(kron(Id, Hnn), Id), Id) +
        kron(kron(kron(Id, Id), Hnn), Id) +
        kron(kron(kron(Id, Id), Id), Hnn)
    )

    return H

def main():
    l0 = 0.720
    l1 = 0.900
    dl = 0.001
    #lamb = np.sqrt(3)/2
    mu = 2
    hamiltonians = {
        5: Hamiltonian_5site,
        6: Hamiltonian_6site,
        7: Hamiltonian_7site,}
    
    for nsites, Hamiltonian in hamiltonians.items():
        print(f"Running for {nsites} sites")
        for lamb in np.linspace(l0, l1, int((l1 - l0) / dl) + 1):
            print("Running", lamb)
            H = Hamiltonian(lamb, mu)
            eigval = eigh(H, eigvals_only=True)
            eigval /= nsites
            result = {'eigval': np.sort(eigval).tolist(), 'lamb': lamb, 'mu': mu, 'nsites': nsites}
            
            # Dump results 
            timestamp = int(time.time())
            with open(f"./ed_explicit_check_open/ED_lamb-{lamb}_mu-{mu}_nsites-{nsites}_{timestamp}.json", "w") as file:
                json.dump(result, file)


if __name__ == "__main__":
    main()