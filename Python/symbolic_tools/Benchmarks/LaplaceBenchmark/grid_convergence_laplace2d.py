from sympy.abc import x
from DataIO.VTIFile import VTIFile
import os
import pwd
import pandas as pd
import matplotlib.pyplot as plt

import numpy as np

# -------- numerical solution ---------------
wd = os.getcwd()
wd = os.path.dirname(wd)  # go level up

x_size = np.array([32, 64, 128, 256])
fig_name = f'LaplaceBenchmark_grid_convergence_from_{x_size[0]}_to_{x_size[-1]}.png'

home = pwd.getpwuid(os.getuid()).pw_dir
main_folder = os.path.join(home, 'DATA_FOR_PLOTS', 'LaplaceBenchmark')


def get_t_mse(folder, n=len(x_size)):
    T_mse = np.zeros(n)
    for i in range(n):

        filename_vtk = f'laplace_template_nx_{x_size[i]}_ny_{x_size[i] + 2}_VTK_P00_00250000.vti'
        filepath_vtk = os.path.join(main_folder, folder, filename_vtk)
        vti_reader = VTIFile(filepath_vtk)

        filename_txt = f'laplace_template_nx_{x_size[i]}_ny_{x_size[i] + 2}_TXT_P00_00250000_T.txt'
        filepath_txt = os.path.join(main_folder, folder, filename_txt)

        T_num_txt = pd.read_csv(filepath_txt, delimiter=" ",  header=None)
        T_num = vti_reader.get("T")
        U = vti_reader.get("U", is_vector=True)

        T_num = np.delete(T_num, 0, axis=0)  # delete first row - extra bc (stops periodicity)

        n_rows, n_columns = T_num.shape
        T_num = np.delete(T_num, (n_rows - 1), axis=0)  # delete last row - extra bc (stops periodicity)

        # -------- analytical solution ---------------
        ySIZE, xSIZE = T_num.shape
        step = 1
        my_fun = -4 * x * (x - xSIZE) / (xSIZE * xSIZE)
        n_fourier = 25
        anal_input = InputForLaplace2DAnalytical(xSIZE, ySIZE, step, my_fun, n_fourier)

        dump_fname = os.path.join(main_folder, f'n_fourier{n_fourier}', f'T_anal_x{xSIZE}y{ySIZE}.npy')

        if os.path.isfile(dump_fname):
            print(f'{dump_fname} found, loading results from disc')
            T_anal = np.load(dump_fname)
            x_grid = np.linspace(0, xSIZE, xSIZE, endpoint=False) + 0.5
            y_grid = np.linspace(0, ySIZE, ySIZE, endpoint=False) + 0.5
            xx, yy = np.meshgrid(x_grid, y_grid)
        else:
            print(f'{dump_fname} not found, starting calculations')
            xx, yy, T_anal = analytical_laplace_2d(anal_input)
            np.save(dump_fname, T_anal)

        T_mse[i] = np.sum((T_anal - T_num) * (T_anal - T_num))/len(T_anal)
        # print(f"T_mse={T_mse} for grid {x_size[i]} x {x_size[i]} [lu]")

    return T_mse


T_err_EQ = get_t_mse('eq_scheme_laplace_template')
T_err_ABB = get_t_mse('abb_laplace_template')

# T_mse_EQ = get_t_mse('eq_scheme_laplace_template_without_wall')
# T_mse_ABB = get_t_mse('abb_laplace_template_without_wall')


print("---------- PLOTTING -------------")


plt.rcParams.update({'font.size': 14})
plt.figure(figsize=(14, 8))

axes = plt.gca()
plt.plot(x_size, T_err_EQ,
         color="black", marker="o", markevery=1, markersize=5, linestyle="", linewidth=2,
         label='Equilibrium scheme')

plt.plot(x_size, T_err_ABB,
         color="black", marker=">", markevery=1, markersize=5, linestyle="", linewidth=2,
         label='Anti-Bounce-Back scheme')


# --------- format y axis --------- #
# yll = 1 * 1E-7
# yhl = 6 * 1E-7
# axes.yaxis.set_major_formatter(xfmt)

# plt.yscale('log')
axes.set_ylim(0, 0.004)
axes.set_yticks(np.linspace(0, 0.004, 5))

# plt.ylim(0, 1.05*max(T_mse_ABB[0], T_mse_EQ[0]))
# plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))  # scilimits=(-8, 8)
# --------- format x axis --------- #
axes.set_xticks(np.insert(x_size, 1, 0))

# plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
# plt.xlim(x_size.min(), x_size.max())


plt.title(f'Laplace Benchmark - Grid Convergence Study\n '
          # r'$x_{range}$' + f'={range}'
          # f'; \t'
          # r'$x_{step}$' + f'={step:.4f}'
          )
plt.xlabel(r'lattice size [lu]')
plt.ylabel(r'$T_{mse}$')
plt.legend()
plt.grid()

fig = plt.gcf()  # get current figure
fig.savefig(fig_name, bbox_inches='tight')
plt.show()
# plt.close(fig)  # close the figure


print("bye")
