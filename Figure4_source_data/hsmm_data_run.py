#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 11:19:22 2019

@author: lcaopcn
"""


import numpy as np
import matplotlib.pyplot as plt

from pybasicbayes.util.text import progprint_xrange
from pyhsmm.basic.distributions import PoissonDuration

import hdf5storage
import pyhsmm_spiketrains.models
reload(pyhsmm_spiketrains.models)

N_iter = 500

# load matlab data
run_obs = hdf5storage.loadmat('pyhsmm_data_run.mat')
rfpa = run_obs['rH'].copy(order='C')
Nc = rfpa.shape[1]
rspk = run_obs['rspk'].copy(order='C')
Nn = rspk.shape[1]
rpos = run_obs['rpos'].copy(order='C')


# create model for fpa
fpa_hsmm = pyhsmm_spiketrains.models.PoissonHDPHMM(N=Nc, K_max=64)
fpa_hsmm.add_data(rfpa)


fpa_lls = []
for itr in progprint_xrange(N_iter):
    fpa_hsmm.resample_model()

    # Collect the log likelihood and predictive log likelihood
    fpa_lls.append(fpa_hsmm.log_likelihood(rfpa))

# Get the inferred state sequence
fpa_hsmm.relabel_by_usage()
fpa_train_inf = fpa_hsmm.stateseqs[0]
fpa_N_used_inf = len(list(fpa_hsmm.used_states))

# Plot the log likelihood over time
plt.figure()
plt.plot(fpa_lls, 'b')
plt.xlabel("Iteration")
plt.ylabel("lfp Log Likelihood")

# Visualize the inferred transition matrices
plt.figure()
plt.imshow(fpa_hsmm.A[:fpa_N_used_inf, :fpa_N_used_inf], interpolation="none", cmap="Greys", vmin=0, vmax=1)
plt.title("lfp Inf. Transitions")


# create model for spk
spk_hsmm = pyhsmm_spiketrains.models.PoissonHDPHMM(N=Nn, K_max=64)
spk_hsmm.add_data(rspk)

spk_lls = []
for itr in progprint_xrange(N_iter):
    spk_hsmm.resample_model()

    # Collect the log likelihood and predictive log likelihood
    spk_lls.append(spk_hsmm.log_likelihood(rspk))

# Get the inferred state sequence
spk_hsmm.relabel_by_usage()
spk_train_inf = spk_hsmm.stateseqs[0]
spk_N_used_inf = len(list(spk_hsmm.used_states))

# Plot the log likelihood over time
plt.figure()
plt.plot(spk_lls, 'b')
plt.xlabel("Iteration")
plt.ylabel("spike Log Likelihood")

# Visualize the inferred transition matrices
plt.figure()
plt.imshow(spk_hsmm.A[:spk_N_used_inf, :spk_N_used_inf], interpolation="none", cmap="Greys", vmin=0, vmax=1)
plt.title("spike Inf. Transitions")

results = {
          "N_iter": N_iter,
          "fpa_lls": fpa_lls,
          "fpa_N_used_inf": fpa_N_used_inf,
          "fpa_train_inf": fpa_train_inf,
          "fpa_Transitions_mat": fpa_hsmm.A,
          "fpa_rates_mat":fpa_hsmm.rates,
          "spk_lls": spk_lls,
          "spk_N_used_inf": spk_N_used_inf,
          "spk_train_inf": spk_train_inf,
          "spk_Transitions_mat": spk_hsmm.A,
          "spk_rates_mat":spk_hsmm.rates,
        }
hdf5storage.savemat('pyhsmm_data_run_result.mat',results)