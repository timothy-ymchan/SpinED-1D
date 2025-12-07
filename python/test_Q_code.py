import numpy as np
from numpy import kron
from scipy.linalg import eigh, expm, eig
from spin_ops import Id, Sx, Sy, Sz
from build_hamiltonian import build_nn_hamiltonian_nosym 
import json 

np.set_printoptions(precision=3)
# Obtain the basis
S = (Sx+Sy+Sz)/np.sqrt(3)
U = expm((-1.j*2*np.pi/3)*S)
_, V = eig(U) # D = V^dagger U V ; Charge -1,0,1 (index: 0,1,2)
charge_table = [-1,0,1]

# Check that the charge assigned is correct 
print(V.conj().T @ U @ V) # Should be diagonal
print(V.conj().T @ S @ V) 

# Define the Q operator on nsites 
def Q_op_Sz(nsites):
    Q0 = (Sx+Sy+Sz)/np.sqrt(3)
    Q = Q0
    for i in range(nsites-1):
        Q = kron(Q,Q0)
    return Q

def Q_op_Z3(nsites):
    Q0 = (Sx+Sy+Sz)/np.sqrt(3)
    Q0 = V.conj().T @ Q0 @ V
    Q = Q0
    for i in range(nsites-1):
        Q = kron(Q,Q0)
    return Q

# Hamiltonian
def A4_Hamiltonian_Sz(lamb,mu):
    H_Heis = kron(Sx,Sx) + kron(Sy,Sy) + kron(Sz,Sz)
    H_q = kron(Sx,Sx)@kron(Sx,Sx) + kron(Sy,Sy)@kron(Sy,Sy) + kron(Sz,Sz)@kron(Sz,Sz)
    H_c = kron(Sx@Sy,Sz) + kron(Sz@Sx,Sy) + kron(Sy@Sz,Sx) + kron(Sy@Sx,Sz) + kron(Sx@Sz,Sy) + kron(Sz@Sy,Sx) + kron(Sx,Sy@Sz) + kron(Sz,Sx@Sy) + kron(Sy,Sz@Sx) + kron(Sx,Sz@Sy) + kron(Sz,Sy@Sx) + kron(Sy,Sx@Sz)
    
    return H_Heis + lamb*H_c + mu*H_q

def A4_Hamiltonian_Z3(lamb,mu):
    H = A4_Hamiltonian_Sz(lamb,mu)
    VV = kron(V,V)
    return VV.conj().T @ H @ VV

def Hamiltonian_Sz(lamb,mu,nsites):
    Hnn = A4_Hamiltonian_Sz(lamb,mu)
    return build_nn_hamiltonian_nosym(Hnn,nsites)

def Hamiltonian_Z3(lamb,mu,nsites):
    Hnn = A4_Hamiltonian_Z3(lamb,mu)
    return build_nn_hamiltonian_nosym(Hnn,nsites)


# Test calculation 
nsites = 7
En = []
Q = []
# lamb_range = np.linspace(0.85,0.88,61)
# for lamb in lamb_range:
#     mu = 2 
#     H_Sz = Hamiltonian_Sz(lamb,mu,nsites)
#     Q_Sz = Q_op_Sz(nsites)
#     eigval, eigvec = eigh(H_Sz)
#     sorted_indices = np.argsort(eigval)
#     eigval = eigval[sorted_indices]
#     eigvec = eigvec[:,sorted_indices]
#     GS = eigvec[:,0]
#     # print(GS.shape)
#     GS_norm = np.einsum('i,i',GS.conj(),GS)
#     Q_expectation_Sz = np.einsum('i,ij,j',GS.conj(),Q_Sz,GS) / GS_norm
#     # print(f"Ground state energy (Sz basis): {eigval[0]:.6f}, <Q> = {Q_expectation_Sz:.6f}")
#     En.append(np.real(eigval)[:5].tolist())
#     Q.append(np.real(Q_expectation_Sz))

# with open(f'A4_Hamiltonian_Q_expectations_nsites-{nsites}.json','w') as f:
#     json.dump({'lambdas': lamb_range.tolist(), 'en': En, 'Q': Q}, f)

# Plot the results 
import matplotlib.pyplot as plt 
data = json.load(open(f'A4_Hamiltonian_Q_expectations_nsites-{nsites}.json','r'))
plt.figure(figsize=(8,6))
plt.plot(data['lambdas'],np.array(data['en'])[:,0],label='GS',marker='o')
plt.plot(data['lambdas'],np.array(data['en'])[:,1],marker='o')
plt.axvline(x=np.sqrt(3)/2)
plt.show()

plt.figure(figsize=(8,6))
plt.plot(data['lambdas'],data['Q'],marker='o')
plt.axvline(x=np.sqrt(3)/2)
plt.show()