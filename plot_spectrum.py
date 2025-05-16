import matplotlib.pyplot as plt 
import json

import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path



def pbc_keys(keys, nsites):
    # Keys are string keys of '0','1',...'nsites-1'
    int_keys = [int(k) for k in keys]
    return np.array([k if k < nsites/2 else k-nsites for k in int_keys])

def get_momentum(result):
    nsites = result['nsites']
    momentum_keys = np.array([k for k in result['eigval'].keys()])
    momentum_vals = pbc_keys(momentum_keys,nsites)
    idc_sort = np.argsort(momentum_vals)
    #print(momentum_vals[idc_sort])
    return momentum_keys[idc_sort], momentum_vals[idc_sort]

def get_xpos(eig_val,k):
    length = len(eig_val)
    #return [k+np.sign(k+0.001)*0.95*l/length for l in range(length)]
    if k != 0:
        x_pos = [k+np.sign(k)*0.95*l/length for l in range(length)]
    else:
        x_pos = []
        dx = 1/length
        in_multiplet = -1
        multiplet_pos = 0
        for i in range(len(eig_val)-1):
            if abs(eig_val[i+1] - eig_val[i]) < 1e-7: # Degeneracy conditions
                in_multiplet += 1
                if in_multiplet % 2 == 0:
                    multiplet_pos += dx
                    x_pos.append(multiplet_pos)
                else:
                    multiplet_pos = -multiplet_pos
                    x_pos.append(multiplet_pos)
            else:
                if in_multiplet != -1 and in_multiplet %2 == 0:
                    multiplet_pos = -x_pos[-1]
                    x_pos.append(multiplet_pos)
                else:
                    x_pos.append(0)
                in_multiplet = -1
                multiplet_pos = 0
                
        x_pos.append(0)
    return x_pos

def get_xpos_tower(eig_val,dx=1/20):
    return [dx*i for i in range(len(eig_val))]
   
# Filename
filename = "./results/nosym/ED_lamb-0.8660254037844386_mu-2_nsites-6_1747406556.json"
# # Plotting the results with symmetry
# with open(filename,"r") as file:
#     result = json.load(file)
#     plt.grid()
#     momentum_keys,momentum_val = get_momentum(result)
#     for i,momentum in enumerate(momentum_keys):
#         all_eigvals = np.array(result['eigval'][momentum])
#         idc_filter = np.where(all_eigvals < 2)
#         eig_vals = all_eigvals[idc_filter]
#         k = momentum_val[i]
#         x_pos = get_xpos(eig_vals,k)
#         plt.scatter(x_pos,eig_vals,s=5)
#     x_labels = [f'k={k}' for k in momentum_val]
#     plt.xticks(momentum_val, x_labels)
#     # Show plot
#     basename = Path(filename).stem
#     plt.savefig(f'./figs/spectrum_{basename}.pdf')


with open(filename,"r") as file:
    result = json.load(file)
    plt.grid()
    all_eigvals = np.array(result['eigval'])
    idc_filter = np.where(all_eigvals < 1.6)
    eig_vals = all_eigvals[idc_filter]
    #print(eig_vals)
    x_pos = get_xpos_tower(eig_vals,dx=1/20)
    plt.scatter(x_pos,eig_vals,s=5)
    x_labels = [f'nosym']
    plt.xticks([0], x_labels)
    basename = Path(filename).stem
    plt.title(basename)
    plt.savefig(f'./figs/nosym/spectrum_{basename}.pdf')