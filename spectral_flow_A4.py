import numpy as np 
import glob
import os
import json 
import matplotlib.pyplot as plt 

def remove_charge(eigvals):
    total = []
    for k,v in eigvals.items():
        total.append(np.array(v))
    total = np.concatenate(total)
    return np.sort(total)
nsites = 3
files = sorted(glob.glob(f'./results/momentumZ3_scan/*nsites-{nsites}*.json'))

K_max = 2
eig = np.zeros((K_max,len(files)))
lambdas = np.zeros(len(files))

for idx, fname in enumerate(files):
    with open(fname, 'r') as f:
        data = json.load(f)
        eig_val = remove_charge(data['eigval'])
        eig_val = sorted(eig_val)
        lamb = data["lamb"]
        lambdas[idx] = lamb
        for i in range(K_max):
            eig[i, idx] = eig_val[i] if len(eig_val) > i else np.nan

# Sort the lambdas and corresponding eigenvalues by the value of lambdas
sort_idx = np.argsort(lambdas)
lambdas_sorted = lambdas[sort_idx]
eig_sorted = eig[:, sort_idx]

# Plot the spectrum for each eigenvalue index
for i in range(K_max):
    plt.scatter(lambdas_sorted, eig_sorted[i], s=2, label=f"eig[{i}]")
#plt.xlim(0.86575,0.86625)
plt.xlabel('lambda')
plt.ylabel('Spectrum')
plt.title(f'nsites: {nsites}')
plt.legend()
plt.tight_layout()
plt.show()
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
