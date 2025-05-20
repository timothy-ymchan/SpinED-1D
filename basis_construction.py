import logging 
from modular_arithmetic import circ_shift, to_base, num_to_digits
from math import log
import os
import shutil
import glob
import json

def get_momentum_projection(s,nsites,momentum,base=3):
    si = s
    min_s = s
    R_to_min = 0
    for R in range(1, nsites+1):
        si = circ_shift(si, shift=1, base=base, ndit=nsites)
        if si < min_s:
            min_s = si 
            R_to_min = R
        if si == s:
            break 
    
    # Check momentum compatibility
    if momentum*R % nsites != 0:
        return -1,-1,-1
    return min_s, R, R_to_min

def get_Zn_charge(state,charge_table,mod_n,nsites,base=3):
    # charge_table is a look up table
    hib_dim = len(charge_table) 
    assert hib_dim == base, "State base should match the length of charge_table"
    tot_charge = 0
    loop_count = 0
    while state > 0:
        tot_charge += charge_table[state % base] 
        state //= base
        loop_count += 1
    tot_charge += (nsites-loop_count)*charge_table[0] # The rest are zeros
    return tot_charge % mod_n


def get_period(s,nsites,base=3):
    R = 1
    si = s
    logging.debug(f'[get_period] s0 {to_base(si,base=base,ndigit=nsites)}')
    while R < nsites:
        sip1 = circ_shift(si,shift=1,base=base,ndit=nsites)
        logging.debug(f'[get_period] s{R} {to_base(sip1,base=base,ndigit=nsites)}')
        if sip1 < s:
            return -1
        elif sip1 == s:
            return R 
        si = sip1 
        R += 1
    return R

def pprint_period_table(period_table,base=3):
    nsites = int(log(max(period_table[1])+1)/log(base))
    print('Number of sites: ',nsites)
    for P in sorted(period_table.keys()):
        reps = period_table[P]
        print(f"Period {P}:")
        for rep in reps:
            print(f"  {to_base(rep,base=base,ndigit=nsites)}")


def construct_period_table(nsites,base=3):
    # Reference: https://physics.bu.edu/~sandvik/perimeter/l07.pdf
    basis = {} # (period, [rep1, rep2, ...]) where repi are representatives of the basis 
    for s in range(0,base**nsites):
        P = get_period(s,nsites)
        if P > 0:
            if P in basis:
                basis[P].append(s)
            else:
                basis[P] = [s]
    assert sum([P*len(B) for P,B in basis.items()]) == base**nsites, f"The total number of basis does not add up to 2^nsites = {2**nsites}"
    
    # Debug message to check the construction is correct 
    if logging.getLogger().isEnabledFor(logging.DEBUG):
        pprint_period_table(period_table=basis)
    return basis 

def clear_cache(cache_dir):
    print(f"Cleaning up {cache_dir}. Removing all contents.")
    for filename in os.listdir(cache_dir):
        file_path = os.path.join(cache_dir, filename)
        if os.path.isfile(file_path) or os.path.islink(file_path):
            os.unlink(file_path)
        elif os.path.isdir(file_path):
            shutil.rmtree(file_path)

def construct_period_table_cache(nsites,cache_dir,base=3):

    # Manage the cache_dir
    if not os.path.exists(cache_dir):
        print(f"Directory does not exist. Creating {cache_dir}")
        os.makedirs(cache_dir)
    else:
        print(f"Directory {cache_dir} exists. Removing all contents.")
        for filename in os.listdir(cache_dir):
            file_path = os.path.join(cache_dir, filename)
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
    
    # Start building basis
    periods = []
    for s in range(0, base**nsites):
        P = get_period(s, nsites)
        if P > 0:
            period_file = os.path.join(cache_dir, f"period_{P}.state")
            mode = 'a' if os.path.exists(period_file) else 'w'
            with open(period_file, mode) as f:
                f.write(f"{s}\n")
            if not (P in periods):
                periods.append(P)
    
    # Write a config.json 
    config_path = os.path.join(cache_dir, "sector_info.json")
    with open(config_path, "w") as f:
        json.dump({"periods":sorted(periods)}, f, indent=2)
    # Print complete message
    print(f"Completed constructing period table cache in {cache_dir}")

def refine_cache_period_table_by_charge(nsites,charge_getter,cache_dir,base=3):
    # Check if cache_dir exists and contains *.state files
    if not os.path.exists(cache_dir):
        raise FileNotFoundError(f"Cache directory '{cache_dir}' does not exist.")
    state_files = glob.glob(os.path.join(cache_dir, "*.state"))
    if not state_files:
        raise FileNotFoundError(f"No '.state' files found in cache directory '{cache_dir}'.")
    
    
    # Read periods from sector_info.json
    config_path = os.path.join(cache_dir, "sector_info.json")
    with open(config_path, "r") as f:
        config = json.load(f)
    periods = config.get("periods", [])

    period_charge = []
    for P in periods:
        state_file = f"{cache_dir}/period_{P}.state"
        charge_file_path = f"{cache_dir}/period_{P}"
        with open(state_file, 'r') as f:
            states = [int(line.strip()) for line in f if line.strip()]
        for state in states:
            state_charge = charge_getter(state,nsites,base)
            charge_file = f"{charge_file_path}-charge_{state_charge}.state"
            mode = 'a' if os.path.exists(charge_file) else 'w'
            with open(charge_file, mode) as cf:
                cf.write(f"{state}\n")
            if (P,state_charge) not in period_charge:
                period_charge.append((P,state_charge))
    
    # Add period_charge to sector_info.json
    with open(config_path, "r") as f:
        config = json.load(f)
    config["period_charge"] = period_charge
    with open(config_path, "w") as f:
        json.dump(config, f, indent=2)
    
    print(f"Finished refining cache period table by charge in {cache_dir}")

def get_cache_charged_momentum_state(cache_dir,momentum,charge,nsites):
    config_path = os.path.join(cache_dir, "sector_info.json")
    with open(config_path, "r") as f:
        config = json.load(f)
    period_charge = config.get("period_charge", [])
    compatible_period = [p for (p,c) in sorted(period_charge) if c == charge and p*momentum % nsites == 0]
    states = []
    for P in compatible_period:
        state_file = os.path.join(cache_dir, f"period_{P}-charge_{charge}.state")
        if os.path.exists(state_file):
            with open(state_file, "r") as f:
                states.extend([int(line.strip()) for line in f if line.strip()])
        else:
            raise FileNotFoundError(f"State file '{state_file}' does not exist.")
    return sorted(states)

    
def check_momentum(state,momentum,nsites,base=3):
    # Check momentum of basis, only return the smallest weight state
    P = get_period(state,nsites,base)
    if P < 0:
        return False
    return momentum*P % nsites == 0 # Check if it is divisible

def check_Zn_charge(state,charge,charge_table,mod_n,nsites,base=3):
    # charge_table is a look up table
    return (get_Zn_charge(state,charge_table,mod_n,nsites,base) - charge) % mod_n == 0

def get_sector_basis(nsites,qnums,qnum_checker,base=3):
    state_list = []
    for state in range(base**nsites):
        if all([qc(state,qn) for (qn,qc) in zip(qnums,qnum_checker)]):
            state_list.append(state)
    return sorted(state_list)

def pprint_charge_str(state,charge_table,nsites,base=3):
    # Pretty print charge explitcitly for debug purpose
    digits = to_base(state,base,nsites)
    out_str = ""
    for d in digits:
        out_str = out_str + f"({charge_table[int(d)]})"
    print(out_str)

# def construct_momentum_basis(nsites,momentum,base=3):
#      # Construct basis set directly to save computation time

# def construct_momentum_table(period_table,base=3):
#     periods = sorted(period_table.keys())
#     nsites = int(log(max(period_table[1])+1)/log(base))
#     print('Number of sites: ',nsites)
#     # construct the period table
#     momentum_table = {} 
#     for P in periods:
#         for m in range(0,P):
#             assert nsites % P == 0, f"Period {P} is not a divisor of nsites {nsites}"
#             k = int(m*nsites/P)
#             if k in momentum_table.keys():
#                 momentum_table[k][P] = period_table[P]
#             else:
#                 momentum_table[k] = {P:period_table[P]}

#     # Check hilbert space size make sense
#     hib_size = sum([sum([len(momentum_table[k][P]) for P in momentum_table[k].keys()]) for k in momentum_table.keys()])
#     assert hib_size == base**nsites,  f"The total number of basis does not add up to base^nsites = {base**nsites}"
#     return momentum_table

def construct_momentum_table(period_table,base=3):
    periods = sorted(period_table.keys())
    nsites = round(log(max(period_table[1])+1)/log(base))
    logging.debug('[construct_momentum_table] number of sites: ',nsites)

    # construct the period table
    momentum_table = {} 
    for P in periods:
        for m in range(0,P):
            assert nsites % P == 0, f"Period {P} is not a divisor of nsites {nsites}"
            k = int(m*nsites/P)
            if k in momentum_table.keys():
                momentum_table[k].append(P)
            else:
                momentum_table[k] = [P]
    
    # Check that the number of basis make sense
    hib_size = sum([sum([len(period_table[P]) for P in momentum_table[k]]) for k in momentum_table.keys()])  
    assert hib_size == base**nsites,  f"The total number of basis does not add up to base^nsites = {base**nsites}"
    return momentum_table


if __name__ == "__main__":
    #construct_period_table_cache(12,"./momentum_basis_N12")
    charge_table = [-1,0,1]
    Z3_charge_getter = lambda s,nsite,base: get_Zn_charge(s,charge_table,3,nsite,base)
    states = get_cache_charged_momentum_state("./momentum_basis_N12",momentum=1,charge=0,nsites=12)
    print(states)
    print(len(states))
    #refine_cache_period_table_by_charge(12,Z3_charge_getter,"./momentum_basis_N12")
    #get_cache_charged_momentum_state(1,0,1,2)
    # # Get charges for different states
    # charge_table = [0,1,2]
    # nsites = 3
    # # for state in range(3**nsites):
    # #     print(to_base(state,3,nsites),': ',get_Zn_charge(state,charge_table,mod_n=2,base=3))
    # Zn_checker = lambda state, charge : check_Zn_charge(state,charge,charge_table,mod_n=3)
    # k_checker = lambda state, momentum: check_momentum(state,momentum,nsites)
    # sector_basis = get_sector_basis(nsites,(2,),(k_checker,))
    # for state in sector_basis:
    #     print(to_base(state,3,ndigit=nsites))
    #  # # Compute period table 
    # pt = construct_period_table(nsites=nsites)
    # pprint_period_table(pt)
    # mt = construct_momentum_table(pt)
    # print(mt)
    # # Get period of various numbers 
    # num = 5
    # print(f'The period of  {to_base(num,3,ndigit=5)} is',get_period(num,5))
    # num = 5
    # print(f'The period of {to_base(num,3,ndigit=6)} is ',get_period(num,6))
    # num =  9 + 1
    # print(f'The period of {to_base(num,3,ndigit=4)} is ',get_period(num,4))

    # num = (9 + 1)*3
    # print(f'The period of {to_base(num,3,ndigit=4)} is ',get_period(num,4))


    # # Compute period table 
    # pt = construct_period_table(nsites=6)
    # pprint_period_table(pt)
    # mt = construct_momentum_table(pt)
    # print(mt)


    # Momentum projection
    # base = 3
    # nsites = 6
    # momentum = 4
    # print(f"States with momentum={momentum}")
    # for s in range(0,base**nsites):
    #     s_rep, R, R_steps = get_momentum_projection(s,nsites,momentum=momentum)
    #     if s_rep >= 0 and R == 6:
    #         print(f"{to_base(s,3,nsites)}; s_rep={to_base(s_rep,3,nsites)}; R={R}; R_steps={R_steps}")