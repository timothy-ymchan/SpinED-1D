from modular_arithmetic import circ_shift, to_base
from basis_construction import get_momentum_projection, get_period, check_momentum, get_sector_basis, construct_period_table, construct_momentum_table, get_Zn_charge
import numpy as np
from numpy import kron
from scipy.linalg import expm, eigh
from spin_ops import Sx,Sy,Sz, Id
import logging 

def find_index(num, num_list):
    # Binary search for num in sorted num_list, return index or -1 if not found
    left, right = 0, len(num_list) - 1
    while left <= right:
        mid = (left + right) // 2
        if num_list[mid] == num:
            return mid
        elif num_list[mid] < num:
            left = mid + 1
        else:
            right = mid - 1
    return None

def build_nn_hamiltonian_nosym(H,nsites):
    # Build periodic hamiltonian from nearest neighbor `H` on `nsites`
    hib_onsite = round(np.sqrt(H.shape[0]))
    assert H.shape == (hib_onsite*hib_onsite, hib_onsite*hib_onsite), f"hib_onsite ({hib_onsite}) does not match that of the shape of H ({H.shape})"
    assert nsites >= 2, "Nearest Neighbour Hamiltonian needs at least 2 sites"

    # Build hamiltonian by direct product 
    H0 = H 
    for i in range(nsites-2):
        H0 = kron(H0,np.eye(hib_onsite))
    H0 = H0.reshape([hib_onsite]*(2*nsites)) # relabel the legs 
    print('H0 has shape', H0.shape)

    # Cycle through all sites and sum the hamiltonian 
    idcs = np.arange(nsites)
    #print(idcs)
    permute_axes = [np.concatenate(((idcs+i) % nsites,(idcs+i) % nsites+nsites)) for i in range(nsites)]
    H_periodic = sum([np.permute_dims(H0,ax) for ax in permute_axes])
    return H_periodic.reshape(hib_onsite**nsites,hib_onsite**nsites)


def apply_nn_hamiltonian(H,m,a,nsites,base=3):
    # Apply hamiltonian `H(m+1,m)`` to the state `a` periodically, with a number of sites `nsites`
    # The site indicies goes from `0` to `nsites-1`, which has a physical layout from right to left
    hib_size = base**2 
    assert H.shape == (hib_size,hib_size), f"Hamiltonian shape {H.shape} does not match that of the hilbert space {(hib_size, hib_size)}"

    # Get the shifted_states, it's head, and it's tail 
    state_shifted = circ_shift(a,-m,base,nsites) # When m = 0, we don't need to shift the states
    state_tail = state_shifted % base**2 # Get the two qudit tail of the state
    state_head = state_shifted - state_tail

    # Get the Hamiltonian coefficients 
    coeff = H[:,state_tail]

    # Get the shifted states 
    b = [circ_shift(state_head + i,m,base,nsites) for i in range(hib_size)]

    return coeff, b # the coefficients, the states

def build_nn_hamiltonian_charged_sector(H,nsites,charged_states,momentum,charge,charge_checker,base=3):
    # Float type 
    float_type = np.complex128
    # Compute crystal momentum 
    k = 2*np.pi *momentum/nsites

    # Build the matrix 
    sector_size = len(charged_states)
    logging.debug(f"Building k={momentum} sector with size {sector_size}x{sector_size}")
    H_sector = np.zeros(shape=(sector_size,sector_size),dtype=float_type)
    
    for a in charged_states:
        Pa = get_period(a,nsites,base)
        for m in range(nsites):
            coeffs, bs = apply_nn_hamiltonian(H,m,a,nsites,base)
            for coeff, b in zip(coeffs,bs):
                b_rep, Pb, b_r = get_momentum_projection(b,nsites,momentum,base)

                if b_rep >= 0 and charge_checker(b,charge): # Only if compatible 
                    assert b_rep in charged_states, f"The state {b_rep} is lot in the list of states {charged_states}"
                    phase = np.exp(-1.j*k*b_r)
                    norm = np.sqrt(Pa/Pb)
                    logging.debug(f"<{to_base(b_rep,base,nsites)}|H_{m}|{to_base(a,base,nsites)}> += {norm*phase*coeff}")
                    idc_out, idc_in = find_index(b_rep,charged_states), find_index(a,charged_states)
                    H_sector[idc_out,idc_in] += norm*phase*coeff
    
    atol = 100*np.finfo(float_type).eps
    if not np.allclose(H_sector, H_sector.conj().T, atol=atol):
        max_diff = np.max(np.abs(H_sector - H_sector.conj().T))
        raise AssertionError(f"H_sector is not Hermitian up to 100 times machine precision {atol}. Max difference: {max_diff}")
    return H_sector

# def build_nn_hamiltonian_sector(H,momentum,momentum_table)
# logging.basicConfig(level=logging.DEBUG)
if __name__ == "__main__":
    nsites = 4
    
    # Period table 
    pt = construct_period_table(nsites)
    mt = construct_momentum_table(pt)
    print(pt)
    print(mt)

    # Build random matrices to test the code
    O1 = np.random.rand(3,3) + 1j * np.random.rand(3,3)
    O1 = 0.5 * (O1 + O1.conj().T)
    O2 = np.random.rand(3,3) + 1j * np.random.rand(3,3)
    O2 = 0.5 * (O2 + O2.conj().T)
    Hnn = kron(O1,O2)#kron(Sx,Sx)+kron(Sy,Sy)

    for momentum in range(0,nsites):
        k_checker = lambda state, momentum: check_momentum(state,momentum,nsites)
        sector_basis = get_sector_basis(nsites,(momentum,),(k_checker,))
        print(f'Sector basis for k={momentum}',sector_basis)
        H = build_nn_hamiltonian_charged_sector(Hnn,nsites,sector_basis,momentum,charge=-1,charge_checker=lambda s,c:True)
        #H = build_nn_hamiltonian_sector(Hnn, nsites, momentum, momentum_table, period_table)
        print("Sector size: ",H.shape)
        eigval = eigh(H, eigvals_only=True)
        eigval_str = ", ".join(f"{v:8.3f}" for v in eigval)
        print(f"k={momentum}:\n\t{eigval_str}")

    # Exact diagonalization without symmetry 
    H = build_nn_hamiltonian_nosym(Hnn,nsites)
    print("Sector size: ", H.shape)
    eigval = eigh(H,eigvals_only=True)
    eigval_str = ", ".join(f"{v:8.3f}" for v in eigval)
    print(f"No symmetry:\n{eigval_str}")


    # print("Checking that the state enumeration works")
    # state_str = "012"
    # nsites = len(state_str)
    # a = digits_str_to_num(state_str,base=3)
    # Rx = expm(1.j*np.pi*Sx)
    # H = kron(Id,Rx)
    # print('Rx matrix:')
    # for row in Rx:
    #     print('\t'.join(f'{elem.real:+.3f}{elem.imag:+.3f}j' for elem in row))
    # print("Input states: ",state_str)
    # coeff, out_states = apply_nn_hamiltonian(H=H,m=-1,a=a,nsites=nsites)
    # print("Possible out states are when applying to m=-1:")
    # for c,b in zip(coeff,out_states):
    #     print(to_base(b,base=3,ndigit=nsites), 'coeff:', c)        

    # print()
    # coeff, out_states = apply_nn_hamiltonian(H=H,m=0,a=a,nsites=nsites)
    # print("Possible out states are when applying to m=0:")
    # for c,b in zip(coeff,out_states):
    #     print(to_base(b,base=3,ndigit=nsites), 'coeff:', c) 

    # print()
    # coeff, out_states = apply_nn_hamiltonian(H=H,m=1,a=a,nsites=nsites)
    # print("Possible out states are when applying to m=1:")
    # for c,b in zip(coeff,out_states):
    #     print(to_base(b,base=3,ndigit=nsites), 'coeff:', c) 





