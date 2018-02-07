#!/usr/bin/python

from binning import Binning
import numpy as np
import sys

zbins = np.arange(0, 3.001, 0.01)
abins = np.arange(0.01, 1.001, 5e-3)

params = {
    'Omega_Lambda': 0,
    'Omega_fld': 0,
    'Omega_smg': -1}


def common_params():
    params.update({'h': np.random.uniform(0.6, 0.8),
                   'Omega_cdm': np.random.uniform(0.15, 0.35)
                  })

if sys.argv[1] == '1':  # Quintessence monomial
    def params_func():
        parameters_smg = [
            np.random.uniform(1, 7),  # "N"
            1., #10**np.random.uniform(-1, 1),  # "V0" tunned
            1e-100,  # phi_prime_ini,
            np.random.uniform(1, 7)  # "phi_ini"
        ]

        params.update({"parameters_smg": str(parameters_smg).strip('[]'),
                       "gravity_model": "quintessence_monomial"})
        common_params()
        return params

    binning = Binning(zbins, abins, 'quintessence_monomial-prueba', sys.argv[2])
    #binning.compute_bins('./quintessence_monomial-w0_wa-201702071315.dat')

elif sys.argv[1] == '2':
    gravity_model = "quintessence_modulus"

    def params_func():
        parameters_smg = [
            1.e-100,  # phi_prime_ini
            np.random.uniform(low=-1, high=1),  # phi_ini
            1.,  # 10 ** np.random.uniform(low=-1, high=1) # V0, tunned
            10 ** np.random.uniform(low=-3, high=-1),  # E_D
            np.random.randint(low=1, high=6),  # p_D, The highest number is 5.
            np.random.randint(low=10, high=21),  # n_max, The highest number is 20.
            np.random.uniform(low=0, high=1)  # alpha
        ]
        parameters_2_smg = list(np.random.standard_normal(1 + parameters_smg[-2]))

        params.update({"parameters_smg": str(parameters_smg).strip('[]'),
                       "parameters_2_smg": str(parameters_2_smg).strip('[]'),
                       "gravity_model": gravity_model})
        common_params()
        return params

    binning = Binning(zbins, abins, gravity_model, sys.argv[2])
    #binning.compute_bins('./quintessence_modulus-w0_wa-201702071819.dat')

elif sys.argv[1] == '3':
    gravity_model = "quintessence_axion"

    def params_func():
        E_F = 10 ** np.random.uniform(low=-3, high=-1)  # E_F

        parameters_smg = [
            1e-100,  # phi_prime_ini
            np.random.uniform(low=-np.pi/E_F, high=np.pi/E_F),  # phi_ini
            1.,  # 10 ** np.random.uniform(low=-1, high=1),  # V0, tunned
            E_F,  # E_F
            10 ** np.random.uniform(low=-3, high=-1),  # E_NP
            np.random.randint(low=10, high=21)  # n_max, The highest number is 20.
        ]

        parameters_2_smg = list(np.random.standard_normal(parameters_smg[-1] - 1))

        params.update({"parameters_smg": str(parameters_smg).strip('[]'),
                       "parameters_2_smg": str(parameters_2_smg).strip('[]'),
                       "gravity_model": gravity_model})
        common_params()
        return params

    binning = Binning(zbins, abins, gravity_model, sys.argv[2])

    #binning.compute_bins('./quintessence_axion-w0_wa-201702071506.dat')

elif sys.argv[1] == '4':
    gravity_model = 'quintessence_eft'

    def params_func():
        E_F = 10 ** np.random.uniform(low=-3, high=-1)

        parameters_smg = [
            1e-100,  # phi_prime_ini
            np.random.uniform(low=-1/E_F, high=1/E_F),  # phi_ini
            1.,  # 10 ** np.random.uniform(low=-1, high=1),  # V0, tunned
            E_F,
            np.random.randint(low=5, high=11),  # n_min, The highest number is 10.
            np.random.randint(low=5, high=11),  # n_Q, The highest number is 10.
        ]

        parameters_2_smg = list(np.random.standard_normal(2 + parameters_smg[-1]))

        params.update({"parameters_smg": str(parameters_smg).strip('[]'),
                       "parameters_2_smg": str(parameters_2_smg).strip('[]'),
                       "gravity_model": gravity_model})
        common_params()
        return params

    binning = Binning(zbins, abins, gravity_model, sys.argv[2])
    #binning.compute_bins('./quintessence_eft-w0_wa-201702081639.dat')

binning.compute_bins_from_params(params_func, int(sys.argv[3]))
