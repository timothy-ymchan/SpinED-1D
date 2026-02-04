# Reproducing https://iopscience.iop.org/article/10.1088/0953-8984/2/26/010 for the sake of understanding SPT
# Construct the 8 site hamiltonian directly
from spin_ops import Sx,Sy,Sz, Id
from numpy import kron 
import numpy as np 
from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh 
import time
import json

def heisenberg_aklt(beta):
    SS = kron(Sx,Sx) + kron(Sy,Sy) + kron(Sz,Sz)
    return SS + beta * SS@SS 

def build_nn_hamiltonian_nosym_obc(H,nsites):
    # Build periodic hamiltonian from nearest neighbor `H` on `nsites`
    hib_onsite = round(np.sqrt(H.shape[0]))
    assert H.shape == (hib_onsite*hib_onsite, hib_onsite*hib_onsite), f"hib_onsite ({hib_onsite}) does not match that of the shape of H ({H.shape})"
    assert nsites >= 2, "Nearest Neighbour Hamiltonian needs at least 2 sites"

    # Build hamiltonian by direct product 
    H0 = H 
    for i in range(nsites-2):
        H0 = kron(H0,np.eye(hib_onsite))
    H0 = H0.reshape([hib_onsite]*(2*nsites)) # relabel the legs 
    #print('H0 has shape', H0.shape)

    # Cycle through all sites and sum the hamiltonian 
    idcs = np.arange(nsites)
    permute_axes = [np.concatenate(((idcs+i) % nsites,(idcs+i) % nsites+nsites)) for i in range(nsites-1)]
    H_open = sum([np.permute_dims(H0,ax) for ax in permute_axes])
    return H_open.reshape(hib_onsite**nsites,hib_onsite**nsites)


def main():
    beta = 0
    for nsites in [4]:
        print(f"ED for L = {nsites}")
        H_nn = heisenberg_aklt(beta)
        H = build_nn_hamiltonian_nosym_obc(H_nn,nsites)
        eigval= eigh(H,eigvals_only=True)#,k=20,which='LM',sigma=0,return_eigenvectors=False)
        result = {'eigval':np.sort(eigval).tolist(),'beta':beta,'nsites':nsites}
        # Dump results 
        timestamp = int(time.time())
        with open(f"./heisenberg_aklt/ED_beta{beta}_nsites-{nsites}_{timestamp}.json","w") as file:
            json.dump(result,file)

        # for lamb in np.linspace(l0,l1,num_pts):
        #     H_nn = A4_Hamiltonian(lamb,mu)
        #     H = build_nn_hamiltonian_nosym_obc(H_nn,nsites)
        #     print('lamb', lamb, 'sector size ', H.shape)
        #     #eigval = eigh(H,eigvals_only=True)
        #     eigval= eigsh(H,k=20,which='LM',sigma=0,return_eigenvectors=False)
        #     result = {'eigval':eigval.tolist(),'lamb':lamb,'mu':mu,'nsites':nsites}
        #     # Dump results 
        #     timestamp = int(time.time())
        #     with open(f"./results/obc-scans/ED_lamb-{lamb}_mu-{mu}_nsites-{nsites}_{timestamp}.json","w") as file:
        #         json.dump(result,file)

if __name__ == "__main__":
    main()