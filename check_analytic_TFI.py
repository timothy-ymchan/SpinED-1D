import numpy as np 
import glob
import os
import json 
import matplotlib.pyplot as plt 
from math import floor

def remove_charge(eigvals):
    total = []
    for k,v in eigvals.items():
        total.append(np.array(v))
    total = np.concatenate(total)
    return np.sort(total)

def analytic_energy(J,g,L):
    #print(floor(L/2))
    if L % 2 == 0:
        k = np.array([(2*n-1)*np.pi/L for n in range(1,floor(L/2)+1)])
        return -2*J*np.sum(np.sqrt((np.cos(k)+g)**2 + np.sin(k)**2)) 
    else:
        print('Odd')
        k_PBC = np.array([(2*n)*np.pi/L for n in range(1,floor(L/2))])
        k_ABC = np.array([(2*n-1)*np.pi/L for n in range(1,floor(L/2)+1)])
        E_PBC = -2*J*np.sum(np.sqrt((np.cos(k_PBC)+g)**2 + np.sin(k_PBC)**2)) + 2*J  + g*J
        E_ABC = -2*J*np.sum(np.sqrt((np.cos(k_ABC)+g)**2 + np.sin(k_ABC)**2)) 
        return  E_PBC,E_ABC#np.minimum(E_PBC,E_ABC)
    #print(len(k))
    


nsites_list = [7]  # Adjust as needed
K_max = 20
for nsites in nsites_list:
    # if nsites % 2 == 1:
    #     continue
    files = sorted(glob.glob(f'./TFI/momentumZ2/*nsites-{nsites}*.json'))
    if not files:
        continue
    eig = np.zeros((K_max, len(files)))
    lambdas = np.zeros(len(files))
    for idx, fname in enumerate(files):
        with open(fname, 'r') as f:
            data = json.load(f)
            eig_val = remove_charge(data['eigval'])
            eig_val = sorted(eig_val)
            lamb = data["g"]
            lambdas[idx] = lamb
            for i in range(K_max):
                eig[i, idx] = eig_val[i] if len(eig_val) > i else np.nan
    sort_idx = np.argsort(lambdas)
    lambdas_sorted = lambdas[sort_idx]
    eig_sorted = eig[:, sort_idx]
    if nsites % 2 == 0:
        En_anal = np.array([analytic_energy(.25,-2*g,nsites) for g in lambdas_sorted])
        plt.plot(lambdas_sorted,En_anal,label='Analytical')
    else:
        E_PBC= np.array([analytic_energy(.25,-2*g,nsites)[0] for g in lambdas_sorted])
        E_ABC= np.array([analytic_energy(.25,-2*g,nsites)[1] for g in lambdas_sorted])
        plt.plot(lambdas_sorted,E_PBC)
        plt.plot(lambdas_sorted,E_ABC)
    #print(En_anal[:2])
    
    plt.plot(lambdas_sorted,eig_sorted[0]*nsites,label='Numerical',ls='',marker='x')
    plt.plot(lambdas_sorted,eig_sorted[1]*nsites,label='Numerical',ls='',marker='x')
plt.legend()
plt.show()