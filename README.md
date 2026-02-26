## Spin ED 1D
- Exact diagonalization code for 1D spin chain with nearest-neighbor interactions in Python and in Julia
- Julia version supports Lanczos methods
- Implemented periodic boundary conditions and crystal momentum conservation
- Implemented total charge conservation (any abelian charge should work, although I only tried $\mathbb Z_2$ and $\mathbb Z_3$)

## Next step
- semi momentum conservation (low priority, my use case does not have parity)
- Use block Lanczos in the Julia version to eliminate fragile behavior when the spectrum is degenerate
- Anti-periodic boundary conditions

## Transverse field ising model demonstration demo 

![tfi_ising](https://github.com/user-attachments/assets/22195e9b-e821-406b-aca9-63c1bb865f18)
