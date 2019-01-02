from sympy.matrices import eye

from SymbolicCollisions.core.cm_symbols import sv, sb, Mraw_D2Q9, NrawD2Q9, S_relax_D2Q9, S_relax_ADE_D2Q9
from SymbolicCollisions.core.sym_col_fun import get_DF, get_m00, \
    get_mom_vector_from_continuous_def, get_mom_vector_from_discrete_def, get_continuous_force_He_MB, get_discrete_force_Guo
from SymbolicCollisions.core.printers import print_u2, print_as_vector, print_ccode
from SymbolicCollisions.core.hardcoded_results import hardcoded_cm_pf_eq, hardcoded_F_cm_pf

print("\n\n=== PRETTY CODE: relax and collide ===\n\n")

pop_in_str = 'f_in'  # symbol defining populations
temp_pop_str = 'temp'  # symbol defining populations
cm_eq_pop_str = 'cm_eq'  # symbol defining populations
F_cm_str = 'F_phi_cm'

# "Modelling incompressible thermal flows using a central-moments-based lattice Boltzmann method" L. Fei et al. 2017
# eq8 : (eye(9)-S)*cm + S*cm_eq + (eye(9)-S/2.)*force_in_cm_space

print("CudaDeviceFunction void relax_and_collide_CM_phase_field("
      "real_t %s[9], "
      "real_t tau, "
      "vector_t u, "
      "vector_t F_phi"
      ") \n{"
      % pop_in_str)

print_u2()
print("real_t %s = 1./tau;" % sv)
# print("real_t bulk_visc = 1./6. ;")
# print("real_t %s = 1./(3*bulk_visc + 0.5);" % sb)
print("real_t %s = omega_bulk;" % sb)  # s_b = 1./(3*bulk_visc + 0.5)
print("")

print_ccode(get_m00(pop_in_str), assign_to='real_t m00')

print("\nreal_t %s[9]; real_t %s[9]; real_t %s[9];\n" % (temp_pop_str, cm_eq_pop_str, F_cm_str))
print("for (int i = 0; i < 9; i++) {\n\t"
      "%s[i] = %s[i];}" % (temp_pop_str, pop_in_str))

populations = get_DF(pop_in_str)
temp_populations = get_DF(temp_pop_str)
cm_eq = get_DF(cm_eq_pop_str)
F_cm = get_DF(F_cm_str)
m = Mraw_D2Q9 * temp_populations

print("\n//raw moments from density-probability functions")
print("//[m00, m10, m01, m20, m02, m11, m21, m12, m22]")
print_as_vector(m, print_symbol=pop_in_str)

print("\n//central moments from raw moments")
cm = NrawD2Q9 * populations
print_as_vector(cm, print_symbol=temp_pop_str, regex=True)

print("\n//collision in central moments space")
print("//calculate equilibrium distributions in cm space")
# print_as_vector_re(get_cm_vector_from_discrete_def(lambda i: m00 * get_gamma(i)), cm_eq_pop_str)
print_as_vector(hardcoded_cm_pf_eq, cm_eq_pop_str, regex=True)  # save time

print("//calculate forces in cm space")
# print_as_vector(get_cm_vector_from_discrete_def(get_discrete_force_Guo_second_order), F_cm_str, regex=True)
# print_as_vector(get_cm_vector_from_continuous_def(get_continuous_force_He_MB), F_cm_str, regex=True)
print_as_vector(hardcoded_F_cm_pf, F_cm_str, regex=True)  # save time
print("//collide")
# cm_after_collision = (eye(9) - S_relax) * temp_populations + S_relax * cm_eq + (eye(9) - S_relax / 2) * F_cm  # eq 8
# Relax 2nd moments, FOI
# cm_after_collision = (eye(9) - S_relax) * temp_populations + S_relax * hardcoded_cm_pf_eq +  hardcoded_F_cm_pf

# Relax 2nd moments, SOI
cm_after_collision = (eye(9) - S_relax_D2Q9) * temp_populations + S_relax_D2Q9 * hardcoded_cm_pf_eq + (eye(9) - S_relax_D2Q9 / 2) * hardcoded_F_cm_pf

# Relax 1st moments, SOI
#cm_after_collision = (eye(9) - S_relax_phi) * temp_populations + S_relax_phi * hardcoded_cm_pf_eq + (eye(9) - S_relax_phi / 2) * hardcoded_F_cm_pf  # eq 8

print_as_vector(cm_after_collision, print_symbol=pop_in_str, regex=True)

print("\n//back to raw moments")
m = NrawD2Q9.inv() * populations
print_as_vector(m, print_symbol=temp_pop_str, regex=True)

print("\n//back to density-probability functions")
populations = Mraw_D2Q9.inv() * temp_populations
print_as_vector(populations, print_symbol=pop_in_str, regex=True)

print("\n}\n")
