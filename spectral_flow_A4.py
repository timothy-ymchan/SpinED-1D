import numpy as np 
import glob
import os
import json 
import matplotlib.pyplot as plt 
from ast import literal_eval as make_tuple

def remove_charge(eigvals):
    total = []
    for k,v in eigvals.items():
        total.append(np.array(v))
    total = np.concatenate(total)
    return np.sort(total)

# def get_next_smallest(eigval_dict):
#     # Get next smallest eigenvalue with charge label. Modify the original dictionary in-place 
#     charge_sectors = sorted(eigval_dict.keys())
#     min_sector = charge_sectors[0]
#     min_eig = eigval_dict[min_sector][0] # Get some arbitary sectors as guess
    
#     # Find the smallest eigenvalue and it's label 
#     for charge_sector in charge_sectors:
#         eig_vals = np.sort(eigval_dict[charge_sector])
#         if eig_vals[0] < min_eig:
#             min_eig = eig_vals[0]
#             min_sector = charge_sector
#         eigval_dict[charge_sector] = eig_vals # Enfore order in eigvals 
    
#     # Pop entries from eigval_dict
#     eigval_dict[min_sector] = eigval_dict[min_sector][1:] # Pop the first element
#     return min_sector, min_eig 
        

nsites = 9
sym_files = sorted(glob.glob(f'./results/period_three_study/*nsites-{nsites}*.json'))

K_max = 2
eig = {}
# eig = np.zeros((K_max,len(sym_files)))
# eig_charges = np.array([[None for j in range(len(sym_files))] for i in range(K_max)],dtype=np.dtypes.StringDType(na_object=None))
lambdas = np.zeros(len(sym_files))

# Loading the eigenvalues
for idx, fname in enumerate(sym_files):
    with open(fname, 'r') as f:
        data = json.load(f)
        eigval_dict = data['eigval']
        lamb = data["lamb"]
        lambdas[idx] = lamb
        for charge_sector in eigval_dict:
            if charge_sector not in eig.keys(): # Create sector list 
                eig[charge_sector] = [[] for i in range(K_max)]
            for i in range(K_max):
                eig[charge_sector][i].append(eigval_dict[charge_sector][i])

# Sort the lambdas and corresponding eigenvalues by the value of lambdas
sort_idx = np.argsort(lambdas)
lambdas_sorted = lambdas[sort_idx]
for charge_sector in eig.keys():
    for i in range(K_max):
        eig[charge_sector][i] = np.array(eig[charge_sector][i])[sort_idx]
    
#Plot the spectrum for each eigenvalue index
threshold_min  = 1.28
threshold_max = 1.3 
epsilon = 0
depsilon = 0.0005
for charge_sector in sorted(eig.keys()):
    for i in range(K_max):
        if np.min(eig[charge_sector][i]) < threshold_min and np.max(eig[charge_sector][i]) < threshold_max:
            k,q = make_tuple(charge_sector)
            k = k if k < nsites /2 else k-nsites
            plt.plot(lambdas_sorted,eig[charge_sector][i]+epsilon,label=f'k={k}, Z3 = {q}; eig[{i}]',ms=2,marker='o')

plt.title(f'nsites {nsites}')
plt.xlabel('Lambda')
plt.ylabel('Energy')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.tight_layout()
plt.grid()
plt.savefig(f'./spectral flows/period_three_study/nsite-{nsites}-flow.pdf')
plt.show()

#nosym_files = sorted(glob.glob(f'./results/nosym_scan/*nsites-{nsites}*.json'))
#nosym_files = sorted(glob.glob(f'./ed_check_open/*nsites-{nsites}*.json'))

# K_max = 2
# eig = np.zeros((K_max,len(nosym_files)))
# lambdas = np.zeros(len(nosym_files))

# for idx, fname in enumerate(nosym_files):
#     with open(fname, 'r') as f:
#         data = json.load(f)
#         eig_val = data['eigval']
#         eig_val = sorted(eig_val)
#         lamb = data["lamb"]
#         lambdas[idx] = lamb
#         for i in range(K_max):
#             eig[i, idx] = eig_val[i] if len(eig_val) > i else np.nan

# # Sort the lambdas and corresponding eigenvalues by the value of lambdas
# sort_idx = np.argsort(lambdas)
# lambdas_sorted = lambdas[sort_idx]
# eig_sorted = eig[:, sort_idx]

# # Plot the spectrum for each eigenvalue index
# for i in range(K_max):
#     plt.scatter(lambdas_sorted, eig_sorted[i], s=2, label=f"nosym eig[{i}]")
# plt.axvline(np.sqrt(3)/2,ls='--',color='grey')
# #plt.xlim(0.86575,0.86625)
# plt.xlabel('lambda')
# plt.ylabel('Spectrum')
# plt.title(f'nsites: {nsites}')
# plt.legend()
# plt.tight_layout()
# plt.savefig(f"./spectral flows/nsite-{nsites}-flow.pdf")
# plt.show()

# # Compute and plot the second derivative of the ground state (eig[0]) with respect to lambda
# d2_eig0 = np.gradient(np.gradient(eig_sorted[0], lambdas_sorted), lambdas_sorted)
# plt.plot(lambdas_sorted, -d2_eig0, label="d² eig[0] / d lambda²", marker='o')
# plt.xlabel('lambda')

# Repeat the analysis for different nsites
# nsites_list = [4,6,8, 10,12,14]  # Adjust as needed
# K_max = 20
# for nsites in nsites_list:
#     files = sorted(glob.glob(f'./TFI/momentumZ2/*nsites-{nsites}*.json'))
#     if not files:
#         continue
#     eig = np.zeros((K_max, len(files)))
#     lambdas = np.zeros(len(files))
#     for idx, fname in enumerate(files):
#         with open(fname, 'r') as f:
#             data = json.load(f)
#             eig_val = remove_charge(data['eigval'])
#             eig_val = sorted(eig_val)
#             lamb = data["g"]
#             lambdas[idx] = lamb
#             for i in range(K_max):
#                 eig[i, idx] = eig_val[i] if len(eig_val) > i else np.nan
#     sort_idx = np.argsort(lambdas)
#     lambdas_sorted = lambdas[sort_idx]
#     eig_sorted = eig[:, sort_idx]
#     d2_eig0 = np.gradient(np.gradient(eig_sorted[0], lambdas_sorted), lambdas_sorted)
#     plt.plot(lambdas_sorted, -d2_eig0, label=f"nsites={nsites}",marker='o')

# plt.grid()
# plt.title('Exact diagonalization TFI')
# plt.xlabel('lambda')
# plt.ylabel('-d²E/dlambda²')
# plt.legend()
# plt.tight_layout()
# plt.show()


# Now gaps_dict contains the gaps for each file, keyed by filename
