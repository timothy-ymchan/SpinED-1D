
from numpy import kron 
import numpy as np
from scipy.linalg import eigh, expm, eig 
from scipy.sparse.linalg import eigsh
from build_hamiltonian import build_nn_hamiltonian_nosym, build_nn_hamiltonian_charged_sector
from basis_construction import check_momentum, check_Zn_charge, get_sector_basis, construct_period_table_cache, refine_cache_period_table_by_charge, get_Zn_charge,get_cache_charged_momentum_state,clear_cache
from modular_arithmetic import to_base
# from spin_ops import Id, Sx, Sy, Sz
import json
import time

def A4_operators():
    # Identity operator 
    Id = np.eye(3)

    # Op and Om
    Op = np.array([
        [1/2, 0, -1j*np.sqrt(3)/2],
        [0, -1, 0],
        [-1j*np.sqrt(3)/2, 0, 1/2]
    ], dtype=complex)

    Om = np.array([
        [1/2, 0,  1j*np.sqrt(3)/2],
        [0, -1, 0],
        [1j*np.sqrt(3)/2, 0, 1/2]
    ], dtype=complex)


    # T matrices
    T1 = np.array([
        [0, -1j, 0],
        [0,  0,  1j],
        [0,  0,  0]
    ], dtype=complex)

    T2 = np.array([
        [0,  0,  1j],
        [0,  0,  0],
        [-1j, 0,  0]
    ], dtype=complex)

    T3 = np.array([
        [0,  0,  0],
        [-1j, 0,  0],
        [0,  1j, 0]
    ], dtype=complex)


    # S matrices
    S1 = np.array([
        [0, 0, 0],
        [1, 0, 0],
        [0, 1, 0]
    ], dtype=float)

    S2 = np.array([
        [-1, 0, 0],
        [ 0, 0, 0],
        [ 0, 0, 1]
    ], dtype=float)

    S3 = np.array([
        [0, -1, 0],
        [0,  0, -1],
        [0,  0,  0]
    ], dtype=float)

    return Id, Op, Om, S1, S2, S3, T1, T2, T3

def op_dot(oplist1,oplist2):
    assert len(oplist1) == len(oplist2)
    return sum([kron(oplist1[i].conj().T,oplist2[i]) for i in range(len(oplist1))],start=np.zeros_like(kron(oplist1[0],oplist2[0])))

def A4_Hamiltonian(J0,J1,J3ss,J3tt,J3st):
    Id,Op,Om,S1,S2,S3,T1,T2,T3 = A4_operators()
    II = kron(Id,Id)
    OO = kron(Op,Om) + kron(Om,Op)
    SS = op_dot([S1,S2,S3],[S1,S2,S3])
    TT = op_dot([T1,T2,T3],[T1,T2,T3])
    ST = op_dot([S1,S2,S3],[T1,T2,T3])
    TS = op_dot([T1,T2,T3],[S1,S2,S3])
    return J0*II + J1*OO + J3ss*SS + J3tt*TT + J3st*(ST+TS)

def Haklt0():
    return A4_Hamiltonian(J0=4/9,J1=1/9,J3ss=5/6,J3tt=1/6,J3st=0)

def Haklt1p():
    return A4_Hamiltonian(J0=4/9,J1=1/9,J3ss=1/3,J3tt=2/3,J3st=np.sqrt(3)/6)

def Haklt1m():
    return A4_Hamiltonian(J0=4/9,J1=1/9,J3ss=1/3,J3tt=2/3,J3st=-np.sqrt(3)/6)

def aklt_triangle(x1,x2,x3):
    assert x1+x2+x3 == 1
    assert x1 >=0 and x2 >=0 and x3 >= 0
    return x1*Haklt0() + x2*Haklt1p() + x3*Haklt1m()

def ED_momentum(x1,x2,x3,nsites):
    Hnn = aklt_triangle(x1,x2,x3)
    
    result = {'eigval':{},'x1':x1,'x2':x2,'x3':x3,'nsites':nsites}
    for k in range(0,nsites):
        k_checker = lambda state, momentum: check_momentum(state,momentum,nsites) # Build basis on the fly to use less memory 
        sector_basis = get_sector_basis(nsites,(k,),(k_checker,))

        H = build_nn_hamiltonian_charged_sector(Hnn,nsites,sector_basis,k,1,lambda x,y:True) # Build sector 
        print(f'Sector size for {k}: {H.shape}')
        eig_vals = eigh(H,eigvals_only=True)
        # eig_vals = eigsh(H,k=20,which='LM',sigma=0,return_eigenvectors=False)
        eig_vals = np.sort(eig_vals)
        result['eigval'][k] = sorted((eig_vals/nsites).tolist())

    # Dump results 
    timestamp = int(time.time())
    with open(f"./aklt_triangle/ED_x1-{x1}_x2-{x2}_x3-{x3}_nsites-{nsites}_{timestamp}.json","w") as file:
        json.dump(result,file)

def main():
    for L in [3,4,5,6,7,8,9]:
        for x in np.linspace(0,1,51):
            ED_momentum(x,1-x,0,L)

if __name__ == "__main__":
    main()