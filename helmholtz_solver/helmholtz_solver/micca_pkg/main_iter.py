import numpy as np
import os

from micca_pkg.micca_flame import geom_1
from micca_pkg import os_utils
from micca_pkg import params

from helmholtz_pkg.flame_transfer_function import state_space
from helmholtz_pkg.mshr import MeshXDMF
from helmholtz_pkg.passive_flame import PassiveFlame
from helmholtz_pkg.active_flame import ActiveFlame
from helmholtz_pkg.eigensolvers import fixed_point_iteration_eps
from helmholtz_pkg.eigenvectors import normalize_eigenvector, normalize_adjoint
from helmholtz_pkg.shape_derivatives import shape_derivatives

# ________________________________________________________________________________


def helmholtz_solver(k, target, **kwargs):

    (mesh_subdir,
     eigenvectors_subdir) = os_utils.create_iter_subdirs(mesh_dir, eigenvectors_dir, i=k)

    mesh_filename = os.path.join(mesh_subdir, 'micca_flame')
    geom_1(mesh_filename, fltk=False, **kwargs)

    geometry = MeshXDMF(mesh_filename, write_xdmf_file=True)
    geometry()
    mesh = geometry.mesh
    subdomains = geometry.subdomains
    boundaries = geometry.boundaries

    # ________________________________________________________________________________

    matrices = PassiveFlame(mesh, boundaries, boundary_conditions,
                            c=params.c,
                            degree=degree)
    matrices.assemble_A()
    matrices.assemble_C()
    A = matrices.A
    C = matrices.C

    D = ActiveFlame(mesh, subdomains,
                    params.x_r, params.rho_amb, params.Q_tot, params.U_bulk, FTF,
                    degree=degree)

    D.assemble_submatrices('direct')

    E = fixed_point_iteration_eps(A, C, D, target**2, i=0, tol=1e-4)

    omega_1, p_1 = normalize_eigenvector(mesh, E, i=0, degree=degree)
    omega_2, p_2 = normalize_eigenvector(mesh, E, i=2, degree=degree)

    # ________________________________________________________________________________

    D.assemble_submatrices('adjoint')

    E_adj = fixed_point_iteration_eps(A, C, D, target**2, i=1, tol=1e-4, problem_type='adjoint')

    omega_adj_1, p_adj_1 = normalize_eigenvector(mesh, E_adj, i=1, degree=degree)
    omega_adj_2, p_adj_2 = normalize_eigenvector(mesh, E_adj, i=3, degree=degree)

    p_adj_norm_1 = normalize_adjoint(omega_1, p_1, p_adj_1, matrices, D)
    p_adj_norm_2 = normalize_adjoint(omega_2, p_2, p_adj_2, matrices, D)

    # ________________________________________________________________________________

    omega = (omega_1 + omega_2) / 2

    results = shape_derivatives(geometry, boundary_conditions,
                                omega, (p_1, p_2), (p_adj_norm_1, p_adj_norm_2), params.c,
                                local=True)

    # ________________________________________________________________________________

    # Save eigenvalues, eigenvectors and shape derivatives

    os_utils.save_eigenvector(eigenvectors_subdir, 'p_1', p_1)
    os_utils.save_eigenvector(eigenvectors_subdir, 'p_2', p_2)

    os_utils.save_eigenvector(eigenvectors_subdir, 'p_adj_1', p_adj_norm_1)
    os_utils.save_eigenvector(eigenvectors_subdir, 'p_adj_2', p_adj_norm_2)

    eigs = {'omega_1': omega_1,
            'omega_2': omega_2,
            'omega_adj_1': omega_adj_1,
            'omega_adj_2': omega_adj_2}

    os_utils.pickle_dump(pickle_dir, 'eigs', eigs, i=k)
    os_utils.save_eigs_as_text(pickle_dir, eigs, i=k)

    os_utils.pickle_dump(pickle_dir, 'shape_derivatives', results, i=k)

    os_utils.pickle_dump(pickle_dir, 'msh_params', kwargs, i=k)

    return omega, results

# ________________________________________________________________________________


(my_dir,
 mesh_dir,
 results_dir,
 eigenvectors_dir,
 pickle_dir) = os_utils.create_dirs()
# shutil.copy('MeshDir/MICCA_flame.geo', mesh_dir)

msh_params = {'R_in_p': .14,
              'R_out_p': .21,
              'l_p': .07,
              'h_b': .0165,
              'l_b': .014,
              'h_pp': .00945,
              'l_pp': .006,
              'h_f': .018,
              'l_f': .006,
              'R_in_cc': .15,
              'R_out_cc': .2,
              'l_cc': .2,
              'l_ec': 0.041,
              'lc_1': 5e-2,
              'lc_2': 1e-2
              }

boundary_conditions = {1: 'Neumann',
                       2: 'Neumann',
                       3: 'Neumann',
                       4: 'Neumann',
                       5: 'Neumann',
                       6: 'Neumann',
                       7: 'Neumann',
                       8: 'Neumann',
                       9: 'Neumann',
                       10: 'Neumann',
                       11: 'Dirichlet'}

degree = 2

# FTF = n_tau(params.N3, params.tau)
FTF = state_space(params.S1, params.s2, params.s3, params.s4)

foo = {'pl_rear': 1,
       'pl_outer': 2,
       'pl_inner': 3,
       'pl_front': 4,
       'b_lateral': 5,
       'b_front': 6,
       'pp_lateral': 7,
       'cc_rear': 8,
       'cc_outer': 9,
       'cc_inner': 10,
       'cc_front': 11
       }

# ________________________________________________________________________________

target_ = 3e3

for k_ in range(5):

    eig, shape_der2 = helmholtz_solver(k_, target_, **msh_params)
    target_ = np.abs(eig)

    # ________________________________________________________________________________

    shape_der2 = dict(zip(list(foo.keys()), list(shape_der2.values())))
    shape_der1 = dict()

    shape_der1['R_in_p'] = - shape_der2['pl_inner'][0].imag
    shape_der1['R_out_p'] = shape_der2['pl_outer'][0].imag
    shape_der1['l_p'] = shape_der2['pl_rear'][0].imag

    # shape_der1['h_b'] = shape_der2['b_lateral'][0].imag
    # shape_der1['h_pp'] = shape_der2['pp_lateral'][0].imag

    shape_der1['R_in_cc'] = - shape_der2['cc_inner'][0].imag
    shape_der1['R_out_cc'] = shape_der2['cc_outer'][0].imag
    shape_der1['l_cc'] = shape_der2['cc_front'][0].imag

    max_key = max(shape_der1, key=lambda x: abs(shape_der1[x]))
    max_value = abs(shape_der1[max_key])

    for key, value in shape_der1.items():
        msh_params[key] -= 2.5e-3 * value / max_value
