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

def read_spectrum(file_path):
    result = []
    with open(file_path,'r') as file:
        data = json.load(file) 
        for entry in data['eigval']:
            qn = entry['qn']
            Z3_charge = qn['charge']
            momentum = format_momentum(qn['k'], data['nsites'])
            for en in sorted(entry['en'])[:5]:
                result.append({'energy': en,'momentum': momentum, 'charge': Z3_charge})
    return pd.DataFrame(result)
# Load all the paths and create a dataframe with paths and simulation params 
def main():
    df = load_simulation_data("./python/results/momentumZ3")
    lamb_list = [0.840,0.8660254037844386,0.880]
    nsites_list = [7,8,9]
    df = df[df['lamb'].isin(lamb_list)] # Filter out specific values of lambdas
    df = df[df['nsites'].isin(nsites_list)] # Filter out specific values of nsites
    # print(df)
    # Group by lamb 
    for lamb, lamb_group in df.groupby('lamb'):
        # Group by nsites
        print(lamb, len(lamb_group), "\t", len(nsites_list))

        mu = lamb_group['mu'].iloc[0]
        filename = f"./python/results/low lying spectrum/spectrum_lamb-{lamb}_mu-{mu}.json"
        assert len(lamb_group) == len(nsites_list), "Mismatched in nsites"

        fig,ax = plt.subplots(nrows=1, ncols=len(lamb_group), figsize=(5*len(lamb_group),4),sharey=True) # I will decide whethet I should share y 

        for i, (nsites, nsites_group) in enumerate(lamb_group.groupby('nsites')):
            ax[i].grid(True)
            ax[i].set_title(f'L={nsites} lambda={lamb:.3f} mu={mu}')
            assert len(nsites_group) == 1, "There should be only one entry per nsites"
            path = nsites_group['path'].iloc[0]
            spectrum = read_spectrum(path)
            gs_energy = spectrum['energy'].min() # Get the minimum energy 
            spectrum['energy'] = spectrum['energy'] - gs_energy # Shift the energy
            
            spectrum = spectrum[np.abs(spectrum['momentum']) <= 2]
            
            charge_colormap = { -1: 'tab:blue', 0: 'tab:orange', 1: 'tab:green' }
            charge_markers = { -1: 's', 0: 'o', 1: '^' }
            for Z3_charge, charge_group in spectrum.groupby('charge'):
                color = charge_colormap[Z3_charge]
                marker = charge_markers[Z3_charge]
                ax[i].scatter(charge_group['momentum'], charge_group['energy'], label=f'Z3={Z3_charge}', color=color,marker=marker)
            #fig, ax = plt.subplots(figsize=())
            ax[i].set_xlabel('k')
            ax[i].set_ylabel('E - E0')
            ax[i].legend()
            
        plt.savefig("./python/results/low lying spectrum/" + f"spectrum_lamb-{lamb}_mu-{mu}.pdf", bbox_inches='tight')

if __name__ == "__main__":
    main()