from SymbolicCollisions.core.printers import print_as_vector
from SymbolicCollisions.core.ContinuousCMTransforms import ContinuousCMTransforms, get_mom_vector_from_continuous_def
from SymbolicCollisions.core.cm_symbols import \
    F3D, dzeta3D, u3D, rho, \
    F2D, dzeta2D, u2D, rho,\
    cs2_thermal
from SymbolicCollisions.core.cm_symbols import moments_dict
from sympy.matrices import Matrix
from sympy import Symbol

import warnings
import time

start = time.process_time()

lattice = 'D2Q9'
ccmt = ContinuousCMTransforms(dzeta3D, u3D, F3D, rho, cs2=cs2_thermal)
# ccmt = ContinousCMTransforms(dzeta3D, u3D, F3D, rho)
# ccmt = ContinousCMTransforms(dzeta2D, u2D, F2D, rho, cs2=cs2_thermal)

print('\n\n// === continous cm === \n ')
print("\n--- EQUILIBRIA ---")

print("\n----------------------- calculate particular moment -------------------------------")

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    # row = moments_dict['D3Q15'][14]
    row = moments_dict['D2Q9'][8]
    moment = ccmt.get_cm(row, ccmt.get_Maxwellian_DF)
    print_as_vector(Matrix([moment]), 'particular_moment')


print("\n----------------------- calculate all moments -------------------------------")
cm_eq = get_mom_vector_from_continuous_def(ccmt.get_Maxwellian_DF,
                                           continuous_transformation=ccmt.get_cm,
                                           moments_order=moments_dict[lattice])
print_as_vector(cm_eq, 'cm_eq')

cm_eq_tot_e = get_mom_vector_from_continuous_def(ccmt.get_total_energy_Maxwellian_DF,
                                           continuous_transformation=ccmt.get_cm,
                                           moments_order=moments_dict[lattice])
print_as_vector(cm_eq_tot_e, 'cm_eq_tot_e')

cm_eq_internal_e = get_mom_vector_from_continuous_def(ccmt.get_internal_energy_Maxwellian_DF,
                                           continuous_transformation=ccmt.get_cm,
                                           moments_order=moments_dict[lattice])
print_as_vector(cm_eq_internal_e, 'cm_eq_internal_e')

print("\n----------------------- without regex -------------------------------")
print_as_vector(cm_eq, 'cm_eq', raw_output=True)

print('\n\n Done in %s [s].'
      % str(time.process_time() - start))
