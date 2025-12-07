import numpy as np 
import glob
import os
import json 
import matplotlib.pyplot as plt 

def compute_gap(result,sector=-1):
    k0 = min(result['eigval'].keys())
    E_min = min(result['eigval'][k0]) # Get ground state energy

    Gap = np.inf
    if sector < 0:
        for k in result['eigval'].keys():
            if k != k0:
                E_k = min(result['eigval'][k])
                if E_k-E_min < Gap:
                    Gap = E_min
    else:
        E_k = np.minimum(result['eigval'][sector])
        Gap = E_k - E_min
    return Gap 

if __name__ == "__main__":
    nsites = []
    gaps = []
    for filename in glob.glob('./results/momentum/*.json'):
        with open(filename,"r") as file:
            result = json.load(file)
            nsites.append(result['nsites'])
            gaps.append(compute_gap(result))
    nsites = np.array(nsites)
    gaps = np.array(gaps)
    print(nsites)
    print(gaps)
    plt.plot(1/nsites,gaps,marker='o')
    plt.xlabel('1/nsites')
    plt.ylabel('Gap')
    plt.show()