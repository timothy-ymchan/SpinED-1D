# Operator storage convention
- Consider a Hamiltonian generated from kronecker product
- The definition of `kron` is `
- `kron(a,b)[k0,k1,...,kN] = a[i0,i1,...,iN] * b[j0,j1,...,jN]`
- where 
- `kt = it * st + jt,  t = 0,...,N`, 
- where `a.shape = (r0,r1,..,rN)` and `b.shape = (s0,s1,...,sN)`

- For 2D, we have a matrix of the following form 
  
$$
A \otimes B = 
\begin{bmatrix}
a_{11}B & a_{12}B & \cdots & a_{1n}B \\
a_{21}B & a_{22}B & \cdots & a_{2n}B \\
\vdots & \vdots & \ddots & \vdots \\
a_{m1}B & a_{m2}B & \cdots & a_{mn}B
\end{bmatrix}
$$

- Therefore, consider the spin index $s_i$ and $s_{i+1}$, then we have the following relationship 
- $\braket{s'_i,s'_{i+1}|H|s_i,s_{i+1}} = H(s'_id+s'_{i+1},s_id+s_{i+1})$ where $d$ is the 1-site hilbert space dimension

## Hamiltonian construction
- Input: 
1. `H(m+1,m)` acting on site `m+1` and `m`
2. `a`, which is the lowest weight state of a specific crystal momentum `k` 

- Step:
- Apply `H(m+1,m)` to state `a` and compute all the weights
  - Suppose `a=[in in-1 ... im+1 im ... i3 i2 i1]`
  - The action `H(m+1,m) = sum jm,jm+1 beta(jm+1,jm)|in.... jm+1 jm ... i1> <in in-1 ... im+1 im ... i3 i2 i1|`
  - Return the `beta(jm+1,jm)` and the list of out states `b(jm+1,jm)`
- For each out state `b(jm+1,jm)` project onto crystal momentum basis with momentum `k`
  - Either the state is incompatible, and we get `0`
  - Either the state is compatible, and we get a phase depending on `R(jm+1,jm)`, i.e. the steps needed to go to the minimal state

- About the phase
  - Suppose `T^r b = b_rep`, then `b(k) = exp(-ikr)b_rep(k)` 
- Input:
1. `S=[s1,s2,...,sn]` are lowest representative states with momentum `k`
2. `H(m+1,m)` is the nearest-neighbor hamiltonian acting on site `m+1` and `m`

- Steps:
- For each state `si` in `S`
  - 