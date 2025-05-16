import logging 
from modular_arithmetic import circ_shift, to_base, change_digit
from math import log

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
    # Get period of various numbers 
    num = 5
    print(f'The period of  {to_base(num,3,ndigit=5)} is',get_period(num,5))
    num = 5
    print(f'The period of {to_base(num,3,ndigit=6)} is ',get_period(num,6))
    num =  9 + 1
    print(f'The period of {to_base(num,3,ndigit=4)} is ',get_period(num,4))

    num = (9 + 1)*3
    print(f'The period of {to_base(num,3,ndigit=4)} is ',get_period(num,4))


    # Compute period table 
    pt = construct_period_table(nsites=6)
    pprint_period_table(pt)
    mt = construct_momentum_table(pt)
    print(mt)


    # Momentum projection
    # base = 3
    # nsites = 6
    # momentum = 4
    # print(f"States with momentum={momentum}")
    # for s in range(0,base**nsites):
    #     s_rep, R, R_steps = get_momentum_projection(s,nsites,momentum=momentum)
    #     if s_rep >= 0 and R == 6:
    #         print(f"{to_base(s,3,nsites)}; s_rep={to_base(s_rep,3,nsites)}; R={R}; R_steps={R_steps}")