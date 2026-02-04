# Construct the 8 site hamiltonian directly
from spin_ops import Sx,Sy,Sz, Id
from numpy import kron 
import numpy as np 
from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh 
import time
import json

def A4_Hamiltonian(lamb,mu):
    H_Heis = kron(Sx,Sx) + kron(Sy,Sy) + kron(Sz,Sz)
    H_q = kron(Sx,Sx)@kron(Sx,Sx) + kron(Sy,Sy)@kron(Sy,Sy) + kron(Sz,Sz)@kron(Sz,Sz)
    H_c = kron(Sx@Sy,Sz) + kron(Sz@Sx,Sy) + kron(Sy@Sz,Sx) + kron(Sy@Sx,Sz) + kron(Sx@Sz,Sy) + kron(Sz@Sy,Sx) + kron(Sx,Sy@Sz) + kron(Sz,Sx@Sy) + kron(Sy,Sz@Sx) + kron(Sx,Sz@Sy) + kron(Sz,Sy@Sx) + kron(Sy,Sx@Sz)
    
    return H_Heis + lamb*H_c + mu*H_q

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
    l0 = 0.84
    l1 = 0.88
    dl = 0.001
    mu = 2
    num_pts = int((l1-l0)/dl)+1
    for nsites in [8,9]:
        print(f"ED for L = {nsites}")
        print(f'Scanning from {l0} to {l1} in {dl} steps (Total points = {num_pts})')
        for lamb in np.linspace(l0,l1,num_pts):
            H_nn = A4_Hamiltonian(lamb,mu)
            H = build_nn_hamiltonian_nosym_obc(H_nn,nsites)
            print('lamb', lamb, 'sector size ', H.shape)
            #eigval = eigh(H,eigvals_only=True)
            eigval= eigsh(H,k=20,which='LM',sigma=0,return_eigenvectors=False)
            result = {'eigval':eigval.tolist(),'lamb':lamb,'mu':mu,'nsites':nsites}
            # Dump results 
            timestamp = int(time.time())
            with open(f"./results/obc-scans/ED_lamb-{lamb}_mu-{mu}_nsites-{nsites}_{timestamp}.json","w") as file:
                json.dump(result,file)

if __name__ == "__main__":
    main()