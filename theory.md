# Translation symmetry

## Definitions
- Let $\mathcal B_N$ be the set of binary bit strings of length $N$
- Let $a\in \mathcal B_N$
- Let $\text{Orb}(a)$ be the orbit of $a$ under translation group action acting on the left 
- Let $R_a = |\text{Orb}(a)|$
- Let $a \sim b$ if $b\in \text{Orb}(a)$ 

## Momentum eigenstates
- Let $a\in \mathcal B_N$, we define 
  - $$\ket{a(k)} = \frac{\sqrt{R_a}}{N} \sum_{r=0}^{N-1} e^{-ikr} T^r \ket{a}$$
  - where $k=2\pi m/N$ for $m=0,1,\cdots, N-1$

- The state has the following properties:

1. $T\ket{a(k)} = e^{ik}\ket{a(k)}$ 
2. If $\ket{b} = T^m \ket{a} \in \text{Orb}(a)$ then $\ket{b(k)} = T^m \ket{a(k)} = e^{ikm} \ket{a(k)}$ (i.e. states in the same orbit differ by at phase)
3. If $b\not \in \text{Orb}(a)$, then $\braket{b|T^r|a} = 0$ for all $r$, so $\braket{b(k)|a(k)} = 0$
4. In particular, if $T^{R_a}\ket{a} = \ket{a}$, we have $\ket{a(k)}=e^{ik R_a}\ket{a(k)}$. Therefore $\ket{a(k)} \neq 0$ only when $kR_a = 0 \mod 2\pi$, i.e. $m R_a = 0 \mod N$
5. Conversely, $\ket{a(k)}$ is non-zero if $k R_a = 0 \mod 2\pi$. This is because $T^r \ket{a}$ are all linearly independent for $0\leq r < R_a$. Therefore, the resulting sum can only be zero if the coefficients of each $T^r\ket{a}$ cancels individually.
6. $\braket{a(k)|a(k')} = \delta_{k,k'}$ 

- Property 2 tell us that to define a basis, we need to set a basis order. We will make the following convention: Define the a total order $\leq$ in the set $\mathcal B_n$, then the smallest state in an orbit will always have phase $=1$

- We now discuss the conversion from $\ket{a(k)}$ back to $\ket{a}$
1. There are exactly $R_a$ solutions to the equation $kR_a = 0 \mod 2\pi$. More precisely, if $k=2\pi m/N$, we have $m=0, N/R_a, 2N/R_a,\cdots (R_a-1) N/R_a$ 
2. To invert the relation, we have $T^m \ket{a} = \frac{1}{\sqrt R_a} \sum_{k=0}^{N-1} e^{ikm} \ket{a(k)}$

## Local Hamiltonian
- A Hamiltonian $H$ is **local** if it is of the form $H = \sum_i H_i$ where $H_i$ is supported on a finite number sites
- A hamiltonian is **translational invariant** $[T,H]$
- More often the Hamiltonian will be in a stronger form, and there exist a decomposition $H=\sum_i H_i$, so that $T^j H_i T^{-j} = H_{i+j}$ and these terms are periodic $H_{i+N} = H_i$
- Consider a prototypical example
  - In general, the product $H_i\ket{a}$ will contain a linear combination of local terms $\alpha_1 \ket{b_1} + \alpha_2 \ket{b_2} + \cdots$
  - Then the symmetrization operator in front will act on these vectors and make them $\ket{b_i(k)}$
  - But given the basis defined according to our convention above, these $\ket{b_i (k)}$ might differ from our basis ket by a phase. So we need to correct for that. 
- We would like to project these Hamiltonian of the stronger form to 


## Basis construction
- We want to group the bit strings according to the equivalence class defined above
- We use the same idea in the lecture notes https://physics.bu.edu/~sandvik/perimeter/l07.pdf
- The basic strategy is the following:
  - Sweep from $s=0$ to $s=2^{N}$ 
  - Get the period of $s$ by a method called $\text{get_period}(s,N)$
  - The period is negative if the following holds:
    - If $s$ is not the smallest element in the orbit, there must exist some $r < R$ such that $T^r s < s$. If that is the case, we return -1
    - If not, we return the period of $s$ as $R$
  - The basis table is a dictionary, in which the keys are the period $P$, and the element is a list of all orbit representatives of with orbit length $P$. A orbit representative is defined to be the smallest element in the orbit.


## Build block sectors 
- Consider the Hamiltonian $H$ that commutes with translation operation 
- We have the following 
$$H\ket{a(k)} = \frac{1}{\sqrt{N_a}} \sum_{m=0}^{N-1}\sum_{r=0}^{N-1} e^{-ikr} T^r H_m\ket{a}$$
- Suppose $H_m$ is two-site operator, i.e. 
$$H_m\ket{a} = \sum_{(i_{m+1},i_m)} \beta(i_{m+1},i_m) \ket{a_{N-1},\cdots,i_{m+1},i_m,\cdots,a_1,a_0}$$
- Let $\omega(i_{m+1},i_m)$ be the shift needed to bring the configuration back to the it's representative states in the table $\text{Rep}(\ket{a_{N-1},\cdots,i_{m+1},i_m,\cdots,a_1,a_0})$, i.e. 
$$T^{\omega(i_{m+1},i_m)}\ket{a_{N-1},\cdots,i_{m+1},i_m,\cdots,a_1,a_0} = \text{Rep}(\ket{a_{N-1},\cdots,i_{m+1},i_m,\cdots,a_1,a_0})$$
- Then we see that we have the following equation 

#### Algorithm for building block sectors 
