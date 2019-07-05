from Benchmarks.ADE.steady_two_layer_cylinder_analytical_2D import PipeWithinPipeNeumann
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm as colormap
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import os
import pwd
from DataIO.VTIFile import VTIFile

# -------- numerical solution ---------------
wd = os.getcwd()
wd = os.path.dirname(wd)  # go level up

reference_lattice_size = 32
gauge = 4
lattice_size = int(gauge * reference_lattice_size)

k = f"0.1666666"
home = pwd.getpwuid(os.getuid()).pw_dir
# k_0.1666666_size_128lu

filename_vtk = f'neumann_bc_k_{k}_size_{int(gauge * reference_lattice_size)}lu_VTK_P00_00010000.vti'
main_folder = os.path.join(home, 'DATA_FOR_PLOTS', 'batch_ruraWrurze_NeumannBC', f'k_{k}_size_{lattice_size}lu_bb')
filepath_vtk = os.path.join(main_folder, filename_vtk)
vti_reader = VTIFile(filepath_vtk)
T_num = vti_reader.get("T")
T_num_slice = T_num[:, :, 1]

# ----------------------- calc dimensions -----------------------

ySIZE, xSIZE = T_num_slice.shape
assert ySIZE == xSIZE == int(gauge * reference_lattice_size)

r0 = gauge * (8 / 2)  # inner radius
r2 = gauge * (30 / 2)  # outer radius

# J0 = 0.22  # heat flux (dT/dr) for r = r0
J0 = 1  # heat flux (dT/dr) for r = r0
T2 = 0  # temperature for r = r2

# ----------------------- compute anal solution ---------------------------
x0 = gauge * (reference_lattice_size / 2)  # center of the pipe
y0 = gauge * (reference_lattice_size / 2)

pwp = PipeWithinPipeNeumann(r0, r2, J0, T2)


class Plate1DNeumann:
    def __init__(self, J, T0):
        """
        :param J: heat flux
        :param T0: temperature at x = 0
        """
        self.J = J
        self.T0 = T0
        
    def get_temperature(self, x):
        return self.J*x + self.T0


x_grid = np.linspace(0, xSIZE, xSIZE, endpoint=False) + 0.5
y_grid = np.linspace(0, ySIZE, ySIZE, endpoint=False) + 0.5
xx, yy = np.meshgrid(x_grid, y_grid)
T_anal = np.zeros((ySIZE, xSIZE))

for i in range(ySIZE):
    for j in range(xSIZE):
        r = pwp.get_r_from_xy(xx[i][j], yy[i][j], x0, y0)
        T_anal[i][j] = pwp.get_temperature_r(r)


T_err_field_eq = T_anal - T_num_slice
# T_L2 = np.sum(abs(T_err_field)

# -------- error norm hacks ---------------
# T_err_field = (T_anal - T_num) / T_anal
# T_err_field[np.isnan(T_err_field)]=0
# T_err_field[np.isinf(T_err_field)]=0
# T_err_field = np.clip(T_err_field, -1, 1)
# np.isnat()

T_mse_eq = np.sum((T_anal - T_num_slice) * (T_anal - T_num_slice)) / len(T_anal)
T_L2_eq = np.sqrt(
    np.sum((T_anal - T_num_slice) * (T_anal - T_num_slice))
    / np.sum(T_anal * T_anal))  # Eq. 4.57

x_slice = np.arange(0, int(xSIZE / 2), 1) + 0.5
T_num_slice = T_num_slice[:, int(xSIZE / 2)]  # take Y slice
T_num_slice = T_num_slice[int(xSIZE / 2):]  # half of it

T_anal = T_anal[:, int(xSIZE / 2)]  # take Y slice
T_anal = T_anal[int(xSIZE / 2):]  # half of it

x = x_grid[:int(xSIZE / 2)]  # half of it

# step = 0.01
# r = np.arange(r0, r2, step) + 0.5
mask = (x > r0) & (x < r2)
r = x[mask]
T_r_anal = np.array([pwp.get_temperature_r(r_) for r_ in r])

###################################################################################################################

fig_name = f'pipe_within_pipe_Neumann_anal_vs_num_J0{J0}_T2{T2}_{xSIZE}lu.png'

# -------------------- make dummy plot --------------------
plt.rcParams.update({'font.size': 14})
plt.figure(figsize=(14, 8))

axes = plt.gca()

plt.plot(r, T_r_anal,
         color="black", marker="", markevery=5, markersize=5, linestyle="-", linewidth=2,
         label='analytical solution')

plt.plot(x, T_num_slice,
         color="black", marker="", markevery=1, markersize=15, linestyle=":", linewidth=2,
         label='current model')


# ------ format y axis ------ #
yll = T_r_anal.min()
# yhl = T_r_anal.max()
axes.set_ylim([yll, 1.05*T2])
# axes.set_yticks(np.linspace(yll, yhl, 5))
# axes.set_yticks(np.arange(yll, yhl, 1E-2))
# axes.set_yticks([1E-4, 1E-6, 1E-8, 1E-10, 1E-12])
# axes.yaxis.set_major_formatter(xfmt)

# plt.yscale('log')
# ------ format x axis ------ #
plt.xlim(0, int(xSIZE / 2))
# plt.xlim(int(xSIZE / 2), xSIZE)

# plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
# plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))  # scilimits=(-8, 8)


plt.title(f'T across pipes,  \n '
          f'{xSIZE}x{xSIZE} [lu]'
          f'\t' r'$T_{MSE}$=' + f'{T_mse_eq:.2e}'
          )

plt.xlabel(r'$r$')
plt.ylabel(r'$Temperature$')
plt.legend()
plt.grid()

fig = plt.gcf()  # get current figure
fig.savefig(fig_name, bbox_inches='tight')
plt.show()

# plt.close(fig)  # close the figure
