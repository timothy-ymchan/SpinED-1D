import numpy as np

# Define identity operator
Id = np.eye(3)

# Define spin-1 operators
Sx = (1 / np.sqrt(2)) * np.array([[0, 1, 0],
                                  [1, 0, 1],
                                  [0, 1, 0]])

Sy = (1 / np.sqrt(2)) * np.array([[0, -1j, 0],
                                  [1j, 0, -1j],
                                  [0, 1j, 0]])

Sz = np.array([[1, 0, 0],
               [0, 0, 0],
               [0, 0, -1]])

# Define raising and lowering operators
Sp = Sx + 1j * Sy 
Sm = Sx - 1j * Sy
