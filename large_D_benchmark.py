# Construct the 8 site hamiltonian directly
from spin_ops import Sx,Sy,Sz, Id
from numpy import kron 
import numpy as np 
from scipy.linalg import eigh
import time
import json

# def A4_Hamiltonian(lamb,mu):
#     H_Heis = kron(Sx,Sx) + kron(Sy,Sy) + kron(Sz,Sz)
#     H_q = kron(Sx,Sx)@kron(Sx,Sx) + kron(Sy,Sy)@kron(Sy,Sy) + kron(Sz,Sz)@kron(Sz,Sz)
#     H_c = kron(Sx@Sy,Sz) + kron(Sz@Sx,Sy) + kron(Sy@Sz,Sx) + kron(Sy@Sx,Sz) + kron(Sx@Sz,Sy) + kron(Sz@Sy,Sx) + kron(Sx,Sy@Sz) + kron(Sz,Sx@Sy) + kron(Sy,Sz@Sx) + kron(Sx,Sz@Sy) + kron(Sz,Sy@Sx) + kron(Sy,Sx@Sz)
    
#     return H_Heis + lamb*H_c + mu*H_q

def LargeD(D):
    H_Heis = kron(Sx,Sx) + kron(Sy,Sy) + kron(Sz,Sz)
    return H_Heis

def Hamiltonian_7site(D):
    Hnn = LargeD(D)

    H = (
        kron(kron(kron(kron(kron(Hnn, Id), Id), Id), Id), Id) + 
        kron(kron(kron(kron(kron(Id, Hnn), Id), Id), Id), Id) +
        kron(kron(kron(kron(kron(Id, Id), Hnn), Id), Id), Id) +
        kron(kron(kron(kron(kron(Id, Id), Id), Hnn), Id), Id) +
        kron(kron(kron(kron(kron(Id, Id), Id), Id), Hnn), Id) +
        kron(kron(kron(kron(kron(Id, Id), Id), Id), Id), Hnn) +
        D*kron(kron(kron(kron(kron(kron(Sz@Sz, Id), Id), Id), Id), Id),Id) + 
        D*kron(kron(kron(kron(kron(kron(Id, Sz@Sz), Id), Id), Id), Id),Id) +
        D*kron(kron(kron(kron(kron(kron(Id, Id), Sz@Sz), Id), Id), Id),Id) +
        D*kron(kron(kron(kron(kron(kron(Id, Id), Id), Sz@Sz), Id), Id),Id) +
        D*kron(kron(kron(kron(kron(kron(Id, Id), Id), Id), Sz@Sz), Id),Id) +
        D*kron(kron(kron(kron(kron(kron(Id, Id), Id), Id), Id), Sz@Sz),Id) +
        D*kron(kron(kron(kron(kron(kron(Id, Id), Id), Id), Id), Id), Sz@Sz)
    )

    return H

def Hamiltonian_6site(D):
    Hnn = LargeD(D)

    H = (
        kron(kron(kron(kron(Hnn, Id), Id), Id), Id) + 
        kron(kron(kron(kron(Id, Hnn), Id), Id), Id) +
        kron(kron(kron(kron(Id, Id), Hnn), Id), Id) +
        kron(kron(kron(kron(Id, Id), Id), Hnn), Id) +
        kron(kron(kron(kron(Id, Id), Id), Id), Hnn) +
        D*kron(kron(kron(kron(kron(Sz@Sz, Id), Id), Id), Id),Id) +
        D*kron(kron(kron(kron(kron(Id, Sz@Sz), Id), Id), Id),Id) +
        D*kron(kron(kron(kron(kron(Id, Id), Sz@Sz), Id), Id),Id) +
        D*kron(kron(kron(kron(kron(Id, Id), Id), Sz@Sz), Id),Id) +
        D*kron(kron(kron(kron(kron(Id, Id), Id), Id), Sz@Sz),Id) +
        D*kron(kron(kron(kron(kron(Id, Id), Id), Id), Id), Sz@Sz)
    )

    return H


def main():
    D0 = -2.0
    D1 = 2.0
    dD = 0.1
    #lamb = np.sqrt(3)/2
    hamiltonians = {
        6: Hamiltonian_6site,
        7: Hamiltonian_7site,}
    
    for nsites, Hamiltonian in hamiltonians.items():
        print(f"Running for {nsites} sites")
        for D in np.linspace(D0, D1, int((D1 - D0) / dD) + 1):
            print("Running", D)
            H = Hamiltonian(D)
            eigval = eigh(H, eigvals_only=True)
            eigval /= nsites
            result = {'eigval': np.sort(eigval).tolist(), 'D': D, 'nsites': nsites}
            
            # Dump results 
            timestamp = int(time.time())
            with open(f"./largeD_benchmark/ED_D-{D}_nsites-{nsites}_{timestamp}.json", "w") as file:
                json.dump(result, file)


if __name__ == "__main__":
    #print(Sz@Sz)
    main()