import numpy as np 
import pandas as pd 
import glob 
import re 
import json 
import matplotlib.pyplot as plt 

def exact_simulation_params(file_name):
    pattern = r"ED_lamb-([0-9.]+)_mu-([0-9.]+)_nsites-([0-9]+)_([0-9]+)\.json"
    match = re.search(pattern, file_name)
    if match:
        lamb = float(match.group(1))
        mu = float(match.group(2))
        nsites = int(match.group(3))
        timestamp = int(match.group(4))
        return lamb, mu, nsites, timestamp
    else:
        raise ValueError("Filename does not match expected pattern.")

def load_simulation_data(folder_path):
    path_table = []
    for path in glob.glob(f"{folder_path}/ED_lamb-*.json"):
        lamb, mu, nsites, timestamp = exact_simulation_params(path)
        path_table.append({
            'path': path,
            'lamb': lamb,
            'mu': mu,
            'nsites': nsites,
            'timestamp': timestamp
        })
    return pd.DataFrame(path_table)

def format_momentum(momentum, nsites):
    return  momentum - nsites if momentum > nsites//2 else momentum

# def read_spectrum(file_path):
#     result = []
#     with open(file_path,'r') as file:
#         data = json.load(file) 
#         return data['eigval']
    
def read_GS(paths):
    evs = []
    for path in paths:
        with open(path,'r') as file:
            data = json.load(file) 
            eigvals = data['eigval'][0]
            evs.append(eigvals)
    return np.array(evs)

# Load all the paths and create a dataframe with paths and simulation params 
def main():
    df = load_simulation_data("./python/results/nosym")
    nsites_list = [5,6,7]
    df = df[df['nsites'].isin(nsites_list)]
    for nsites, nsites_group in df.groupby('nsites'):
        nsites_group = nsites_group.sort_values(by='lamb')
        nsites_group['GS'] = read_GS(nsites_group['path'])
        plt.plot(nsites_group['lamb'], nsites_group['GS'], marker='o', label=f'nsites={nsites}')
        plt.show()    

if __name__ == "__main__":
    main()