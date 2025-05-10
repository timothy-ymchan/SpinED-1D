# Build momentum basis 
# We want build a lookup table that does the following:
# It is a dictionary, with keys equal to the period 
# The values are the smallest value of the equivalency class
import logging 
#logging.basicConfig(level=logging.DEBUG)

def circ_shift(num,shift,nbit):
    maxval = 2**nbit
    assert num < maxval, f"Overflow Error: {num} >= {maxval} in circ_shift"
    shift = shift % nbit 
    mask = 2**nbit - 1 # String of 1s of length nbit 
    return ((num << shift) | (num >> (nbit - shift))) & mask

def get_period(s,nsites):
    R = 1
    si = s
    logging.debug(f'[get_period] s0 {format(si, f'0{nsites}b')}')
    while R < nsites:
        sip1 = circ_shift(si,1,nsites)
        logging.debug(f'[get_period] s{R} {format(sip1, f'0{nsites}b')}')
        if sip1 < s:
            return -1
        elif sip1 == s:
            return R 
        si = sip1 
        R += 1
    return R

def pprint_period_table(period_table):
    nsites = int(max(period_table[1]).bit_length())
    for P in sorted(period_table.keys()):
        reps = period_table[P]
        print(f"Period {P}:")
        for rep in reps:
            print(f"  {format(rep, f'0{nsites}b')}")


def construct_period_table(nsites):
    # Reference: https://physics.bu.edu/~sandvik/perimeter/l07.pdf
    basis = {} # (period, [rep1, rep2, ...]) where repi are representatives of the basis 
    for s in range(0,2**nsites):
        P = get_period(s,nsites)
        if P > 0:
            if P in basis:
                basis[P].append(s)
            else:
                basis[P] = [s]
    assert sum([P*len(B) for P,B in basis.items()]) == 2**nsites, f"The total number of basis does not add up to 2^nsites = {2**nsites}"
    
    # Debug message to check the construction is correct 
    if logging.getLogger().isEnabledFor(logging.DEBUG):
        pprint_period_table(period_table=basis)
    return basis 

def construct_momentum_table(period_table):
    periods = sorted(period_table.keys())
    nsites = int(max(period_table[1]).bit_length())

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
    assert hib_size == 2**nsites,  f"The total number of basis does not add up to 2^nsites = {2**nsites}"
    return momentum_table

def build_momentum_sector(Hi,momentum,momentum_table,period_table):
    # Build dense matrix for a momentum sector 
    periods = sorted(momentum_table[momentum])
    M = sum([len(period_table[P]) for P in periods]) 
    for P in periods:
        for state in period_table[P]:
            pass 

if __name__ == "__main__":
    nsites = 5 

    # Try to get period 
    print('The period of  00001 is',get_period(1,5))
    print('The period of 000001 is ',get_period(1,6))

    print('The period of 001001 is ',get_period(0b001001,6))
    print('The period of 010010 is ',get_period(0b010010,6))

    print('The period of 00101 is ',get_period(0b00101,5))
    print('The period of 01001 is ',get_period(0b01001,5))


    # # Constructing the nsite basis
    
    print('Constructing period for nsites=5')
    bt = construct_period_table(5)
    pprint_period_table(bt)
    print('Momentum table ',construct_momentum_table(bt))


    print('Constructing period basis for nsites=6')
    bt = construct_period_table(6)
    pprint_period_table(bt)
    print('Momentum table ', construct_momentum_table(bt))