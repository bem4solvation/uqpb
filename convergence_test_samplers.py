# Convergence of QMC versus MC
# A simple script to compare the convergence of the mean value
# We run the same experiments n_runs
# For a given run, we generate coeffs \sim U[0,1]^natom

import numpy as np
from sampler import sample
from matplotlib import pyplot as plt

# Create a number of runs
n_runs = 2
n_test = 5000
n_atom = 4

M = n_test

# LHS not

# samplers = ['pseudo', 'LHS', 'Halton', 'Hammersley', 'Sobol']

# Sobol up to dimension 40
# LHS not so efficient and SLOW
# For very high dimension, QCM not interesting
# For Sobol: We added randomness -> not the same results  for each run.

samplers = ["pseudo", "LHS", "Halton", "Hammersley", "Sobol"]

means_sampler = []
for sampler in samplers:
    print(sampler)
    means = []
    for n_run in range(n_runs):
        coeffs = sample(n_test, n_atom, sampler=sampler)

        mean = []
        for m in range(M):
            mean.append(coeffs[:m].mean())

        mean = np.array(mean)
        means.append(mean)

    means = np.array(means)
    means_sampler.append(means)


# Plto Mean v/s M
fig, axs = plt.subplots(3, 2, figsize=(10, 10), sharex=True, sharey=True)

for i in range(len(samplers)):
    ax = axs.ravel()[i]
    for n_run in range(n_runs):
        ax.plot(means_sampler[i][n_run])

    ax.set_title(samplers[i])
    ax.legend()

plt.savefig("MeanM.pdf")


# Plot the error vs/ M

fig, axs = plt.subplots(3, 2, figsize=(10, 10), sharex=True, sharey=True)

for i in range(len(samplers)):
    means = means_sampler[i]
    ax = axs.ravel()[i]
    for n_run in range(n_runs):
        y = np.abs(means[n_run] - 0.5)
        x = np.arange(M).astype(float)

        ax.loglog(x, y, color="k")
        ax.loglog(np.abs(np.mean(means, axis=0) - 0.5), color="g")
        ax.loglog(x, 2 * x**-1, color="grey")
        ax.loglog(x, 2 * x ** (-1 / 2), color="grey")

    ax.set_title(samplers[i])
    ax.legend()
plt.savefig("ErrorM.pdf")