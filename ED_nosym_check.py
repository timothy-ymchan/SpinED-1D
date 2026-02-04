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

def Hamiltonian_8site(lamb,mu,open=False):
    Hnn = A4_Hamiltonian(lamb,mu)
    k = 1 if not open else 0
    print("k is ",k)
    H = (
        kron(kron(kron(kron(kron(kron(Hnn, Id), Id), Id), Id), Id), Id) + 
        kron(kron(kron(kron(kron(kron(Id, Hnn), Id), Id), Id), Id), Id) +
        kron(kron(kron(kron(kron(kron(Id, Id), Hnn), Id), Id), Id), Id) +
        kron(kron(kron(kron(kron(kron(Id, Id), Id), Hnn), Id), Id), Id) +
        kron(kron(kron(kron(kron(kron(Id, Id), Id), Id), Hnn), Id), Id) +
        kron(kron(kron(kron(kron(kron(Id, Id), Id), Id), Id), Hnn), Id) +
        kron(kron(kron(kron(kron(kron(Id, Id), Id), Id), Id), Id), Hnn) + 

        k*kron(kron(kron(kron(kron(kron(kron(Sx,Id),Id),Id),Id),Id),Id),Sx) + 
        kron(kron(kron(kron(kron(kron(kron(Sy,Id),Id),Id),Id),Id),Id),Sy) + 
        kron(kron(kron(kron(kron(kron(kron(Sz,Id),Id),Id),Id),Id),Id),Sz) + 


        k*mu*(
            kron(kron(kron(kron(kron(kron(kron(Sx@Sx,Id),Id),Id),Id),Id),Id),Sx@Sx) + 
            kron(kron(kron(kron(kron(kron(kron(Sy@Sy,Id),Id),Id),Id),Id),Id),Sy@Sy) + 
            kron(kron(kron(kron(kron(kron(kron(Sz@Sz,Id),Id),Id),Id),Id),Id),Sz@Sz)
            ) +

        k*lamb * (
        kron(kron(kron(kron(kron(kron(kron(Sz, Id), Id), Id), Id), Id), Id), Sx @ Sy) + #
        kron(kron(kron(kron(kron(kron(kron(Sy, Id), Id), Id), Id), Id), Id), Sz @ Sx) + #
        kron(kron(kron(kron(kron(kron(kron(Sx, Id), Id), Id), Id), Id), Id), Sy @ Sz) + #
        kron(kron(kron(kron(kron(kron(kron(Sx, Id), Id), Id), Id), Id), Id), Sz @ Sy) + #
        kron(kron(kron(kron(kron(kron(kron(Sz, Id), Id), Id), Id), Id), Id), Sy @ Sx) + #
        kron(kron(kron(kron(kron(kron(kron(Sy, Id), Id), Id), Id), Id), Id), Sx @ Sz) + #

        kron(kron(kron(kron(kron(kron(kron(Sx @ Sy, Id), Id), Id), Id), Id), Id), Sz) + #
        kron(kron(kron(kron(kron(kron(kron(Sz @ Sx, Id), Id), Id), Id), Id), Id), Sy) + # 
        kron(kron(kron(kron(kron(kron(kron(Sy @ Sz, Id), Id), Id), Id), Id), Id), Sx) + # 
        kron(kron(kron(kron(kron(kron(kron(Sy @ Sx, Id), Id), Id), Id), Id), Id), Sz) + # 
        kron(kron(kron(kron(kron(kron(kron(Sx @ Sz, Id), Id), Id), Id), Id), Id), Sy) + # 
        kron(kron(kron(kron(kron(kron(kron(Sz @ Sy, Id), Id), Id), Id), Id), Id), Sx) # 
        )
    )
    return H

def main():
    l0 = 0.860
    l1 = 0.870
    dl = 0.0005
    #lamb = np.sqrt(3)/2
    mu = 2
    nsites = 8
    for lamb in np.linspace(l0,l1,int((l1-l0)/dl)+1):
        print("Running", lamb)
        H = Hamiltonian_8site(lamb,mu,open=False)
        eigval = eigh(H,eigvals_only=True)
        eigval /= nsites
        result = {'eigval':np.sort(eigval).tolist(),'lamb':lamb,'mu':mu,'nsites':nsites}
        
        # Dump results 
        timestamp = int(time.time())
        with open(f"./ed_check_periodic/ED_lamb-{lamb}_mu-{mu}_nsites-{nsites}_{timestamp}.json","w") as file:
            json.dump(result,file)

#kron(Hnn,Id,Id,Id,Id,Id,Id) + kron(Id,Hnn,Id,Id,Id,Id,Id) + kron(Id,Id,Hnn,Id,Id,Id,Id) + kron(Id,Id,Id,Hnn,Id,Id,Id) + kron(Id,Id,Id,Id,Hnn,Id,Id) + kron(Id,Id,Id,Id,Id,Hnn,Id) + kron(Id,Id,Id,Id,Id,Id,Hnn) + kron(Sx,Id,Id,Id,Id,Id,)
main()