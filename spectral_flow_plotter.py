import numpy as np 
import glob
import os
import json 
import matplotlib.pyplot as plt 
from ast import literal_eval as make_tuple
import argparse
import importlib
# Sample command:
# spectral_flow_plotter.py --in-pattern "./results/period_three_study/*_nsites-6*.json" --num-eigs 2 --scan-param "lamb" --out-dir "./" --charged --eigs-range "(-1,1.44)"


def remove_charge(eigvals):
    total = []
    for k,v in eigvals.items():
        total.append(np.array(v))
    total = np.concatenate(total)
    return np.sort(total)


def main():
    parser = argparse.ArgumentParser(description="Plotting spectral flows from JSON files")
    parser.add_argument('--in-pattern', type=str, required=True, help='Glob pattern for input files')
    parser.add_argument('--charged', action='store_true', help='Include charge sectors')
    parser.add_argument('--out-dir', type=str, required=True, help='Output directory for plots')
    parser.add_argument('--num-eigs', type=int, default=20, help='Number of eigenvalues to read from each file')
    parser.add_argument('--scan-param', type=str, default='lamb', help='Parameter to scan (e.g., lamb or mu)')
    parser.add_argument('--charge-resolver', type=str, default=None, help='Charge resolver for custom charge sectors labeling')
    parser.add_argument('--eigs-range', type=make_tuple, default=None, help='Range of eigenvalues to plot (min,max)')
    parser.add_argument('--per-site', action='store_true', help='Divide energies by nsites (per-site energies)')

    # If a charge resolver is specified, import it from this file
    charge_resolver_func = None
    if parser.parse_known_args()[0].charge_resolver is not None:
        charge_resolver_func = globals()[parser.parse_known_args()[0].charge_resolver]

    args = parser.parse_args()

    files = sorted(glob.glob(args.in_pattern))
    if not files:
        print(f"No files matched pattern: {args.in_pattern}")
        return

    os.makedirs(args.out_dir, exist_ok=True)

    # Alert users about charged flag 
    if args.charged:
        print("Info: --charged flag is set. Charge sectors will be included in the plot.")
    else:
        print("Info: --charged flag is not set. Charge sectors will be removed and eigenvalues will be combined.")

    eigs = {}
    scan_values = []
    nsites = None 
    for fname in files:
        with open(fname, 'r') as f:
            content = json.load(f)
        eigval_dict = content['eigval'] 

        # Check nsites consistency
        if nsites is None:
            nsites = content['nsites']
        else:
            assert nsites == content['nsites'], f"nsites mismatch in file {fname}: expected {nsites}, got {content['nsites']}"

        # Extract the scan parameter value and append to scan_values
        if args.scan_param not in content:
            raise ValueError(f"Scan parameter '{args.scan_param}' not found in file {fname}")
        scan_value = content[args.scan_param]
        scan_values.append(scan_value)

        # If eigval_dict is not a dict, convert it to one with a default key
        if not isinstance(eigval_dict, dict):
            eigval_dict = {'uncharged': np.array(eigval_dict)} 

        # If --charged flag is not set, remove charge sectors
        if not args.charged:
            eigval_dict = {'uncharged': remove_charge(eigval_dict)}

        # Put eigenvalues into eigs dictionary named by charge sector and filtered by num_eigs
        for charge_sector in eigval_dict:
            if charge_sector not in eigs:
                eigs[charge_sector] = [[] for _ in range(args.num_eigs)]
            eigval_dict[charge_sector] = np.sort(eigval_dict[charge_sector]) # Enforce order in eigvals
            for i in range(args.num_eigs):
                evals = eigval_dict[charge_sector][i]
                if args.per_site:
                    evals /= nsites
                eigs[charge_sector][i].append(evals)

    # Print the sizes of each charge sector for debugging
    for charge_sector in eigs:
        print(f"Charge sector {charge_sector} has {len(eigs[charge_sector][0])} eigenvalues")

    # Sort the scan values and corresponding eigenvalues by the value of scan_values
    sort_idx = np.argsort(scan_values)
    scan_values_sorted = np.array(scan_values)[sort_idx]
    for charge_sector in eigs:
        for i in range(args.num_eigs):
            eigs[charge_sector][i] = np.array(eigs[charge_sector][i])[sort_idx]
    
    # Plotting
    # Assign a color to each charge sector using the default color cycle
    color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
    charge_sectors = sorted(eigs.keys())
    sector_colors = {sector: color_cycle[i % len(color_cycle)] for i, sector in enumerate(charge_sectors)}

    marker_styles = ['o', 's', '^', 'v', 'D', 'P', '*', 'X', '<', '>']

    plt.figure(figsize=(10, 6))
    for i in range(args.num_eigs):
        for j, charge_sector in enumerate(charge_sectors):
            # Check if the eigenvalues are within the specified range
            eig_min, eig_max = -np.inf, np.inf
            if args.eigs_range is not None:
                eig_min, eig_max = args.eigs_range
            sector_min = np.min(eigs[charge_sector][i])
            sector_max = np.max(eigs[charge_sector][i])

            if args.eigs_range is None or (sector_min >= eig_min and sector_max <= eig_max):
                plt.plot(
                    scan_values_sorted,
                    eigs[charge_sector][i],
                    label=f'charge {charge_sector} [eig {i+1}]',
                    marker=marker_styles[i % len(marker_styles)],
                    markersize=8,
                    linestyle='-',
                    markerfacecolor='none',  # hollow marker
                    markeredgewidth=1.5,
                    color=sector_colors[charge_sector]
                )
            
    plt.xlabel(args.scan_param)
    plt.ylabel('Eigenvalue')
    plt.title(f'Spectral Flow for Eigenvalue {i+1} (N={nsites})')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.grid(True)
    plt.savefig(os.path.join(args.out_dir, f'largeD_obc_ed_explicit_spectral_flow_nsites_{nsites}.pdf'))
    #plt.savefig(os.path.join(args.out_dir, f
    #plt.show()
    plt.close()


if __name__ == "__main__":
    main()