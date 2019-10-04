from sympy.matrices import eye
from SymbolicCollisions.core.cm_symbols import e_D3Q27, moments_dict
from SymbolicCollisions.core.printers import print_as_vector, get_print_symbols_in_indx_notation
from sympy import Matrix

import re
import os
import pwd

from SymbolicCollisions.core.cm_symbols import omega_ade, omega_b, omega_v, m00
from SymbolicCollisions.core.cm_symbols import Force_str as F_str
from SymbolicCollisions.core.cm_symbols import dynamic_import, moments_dict
from SymbolicCollisions.core.DiscreteCMTransforms import get_m00
from SymbolicCollisions.core.printers import print_as_vector, get_print_symbols_in_m_notation
from SymbolicCollisions.core.printers import print_u2, print_sigma_cht
from SymbolicCollisions.core.MatrixGenerator import get_m_order_as_in_r, get_e_as_in_r, MatrixGenerator
from sympy.matrices import Matrix
import numpy as np
import pandas as pd

m_seed = [0, 1, 2]
rmoments_order = get_m_order_as_in_r(m_seed, m_seed, m_seed)
q, d = rmoments_order.shape

moments_order = np.array(moments_dict[f'D{d}Q{q}_tclb'])
print(f"order of moments | rmoments: \n "
      f"{pd.concat([pd.DataFrame.from_records(moments_order),pd.DataFrame.from_records(rmoments_order)], axis=1)}")

# e_seed = [0, 1, -1]
# ex_D3Q27new, ey_D3Q27new, ez_D3Q27new, e_D3Q27new = get_e_as_in_r(e_seed, e_seed, e_seed)

# hardcode 'e' vectors as in WMRT matrix: ex=Mraw[1,:], ey = Mraw[2,:], ez=Mraw[3,:]
ex_D3Q27new = Matrix([0 , 1,  -1,  0,  1,  -1,  0,  1,  -1,  0,  1,  -1,  0,  1,  -1,  0,  1,  -1,  0,  1,  -1,  0, 1,  -1,  0,  1,  -1])
ey_D3Q27new = Matrix([0,  0,  0,  1,  1,  1,  -1,  -1,  -1,  0,  0,  0,  1,  1,  1,  -1,  -1,  -1,  0,  0,  0,  1,  1,  1,  -1,  -1,  -1])
ez_D3Q27new = Matrix([0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1])

e_D3Q27new = ex_D3Q27new.col_insert(1, ey_D3Q27new)
e_D3Q27new = e_D3Q27new.col_insert(2, ez_D3Q27new)

print(f"lattice velocities - e: \n {np.array(e_D3Q27new)}")


def print_matrix_to_file(matrix, path):
    print(f'writing to {path}')
    with open(path, "w") as f:
        rows = matrix.tolist()
        for row in rows:
            # row = [str(round(i, 1)) for i in row]
            # row = [str(i).rstrip('0') for i in row]
            # row = [str(i).rstrip('.') for i in row]
            f.write(str(row) + '\n')

my_dir = os.path.dirname(os.path.abspath(__file__))

matrixGenerator = MatrixGenerator(ex_D3Q27new, ey_D3Q27new, ez_D3Q27new, moments_order)
Mraw = matrixGenerator.get_raw_moments_matrix()
Nraw = matrixGenerator.get_shift_matrix(Mraw.inv())

print_matrix_to_file(Mraw, os.path.join(my_dir, "Mraw.txt"))
print_matrix_to_file(Mraw.inv(), os.path.join(my_dir, "Mraw_inv.txt"))
print_matrix_to_file(Nraw, os.path.join(my_dir, "Nraw.txt"))
print_matrix_to_file(Nraw.inv(), os.path.join(my_dir, "Nraw_inv.txt"))


wM = np.genfromtxt("WMRT_matrix_v2.csv", delimiter='\t')

wMinv = np.linalg.inv(wM)
#wMinv2 = Matrix(wM).inv() #TODO: Why do they differ?! np.linalg.inv(wM)

print_matrix_to_file(wM, os.path.join(my_dir, "wM.txt"))
print_matrix_to_file(wMinv, os.path.join(my_dir, "wMinv.txt"))

wN = matrixGenerator.get_shift_matrix(wMinv)
# wN2 = matrixGenerator.get_shift_matrix2(wMinv)
# wN2inv = np.linalg.inv(wN2)
print_matrix_to_file(wN, os.path.join(my_dir, "wN.txt"))
print_matrix_to_file(wN.inv(), os.path.join(my_dir, "wNinv.txt"))

print("BYE")
