### WHAM output

All of the processed simulation results can be found in .txt files in the `wham/` folder. Each of the running estimated properties has corresponding errors that are calculated with a block-averaging procedure.

* [runav_rate.txt](#the-transition-rate)
  The running estimate of the transition rate, meaning an estimate of the rate after each accepted path. errRATE.txt contains the error estimates.

* [Pcross.txt](#the-crossing-probability)
  The estimate of the total crossing probability as a function of the order parameter.

* [ploc_unscaled.txt](#the-crossing-probability)
  The local crossing probabilities in each of the plus ensembles ([0+], [1+], ...,).

* runav_flux.txt
  The running estimate of the flux through the first interfaces $\lambda_0$. errFLUX.txt contains the errors.

* runav_Ptot.txt
  The running estimate of the total crossing probability. This is the last number of contained in Pcross.txt, estimated after each accepted path. errPtot.txt contains the errors.

* runav_ploc.txt
  The running average of the local crossing probabilities in each ensemble. These are the values where the curves from ploc_unscaled.txt reach the neighboring interface. errploc.txt contains the errors.
* pathdistr.txt
  A distribution of the weighted path lengths in each ensemble.
* pathlengths.txt
  The average path length in each ensemble.
* runav_L0.txt
  The running averages of paths in either [0-] or [0+], which are used for the flux calculation. errL0.txt contains the errors.
* ploc_pointmatch.txt
  The same as ploc_unscaled.txt, except that the values are scaled by point matching, which is used to construct the total crossing probabilities in Pcross.txt.
* ploc_WHAM.txt
  Same as above, but scaled with the WHAM weights, which is more accurate than point-matching.

The first line of these files tells you what the different columns are.

#### The transition rate

The running *estimate* of the transition rate can be found in the file `wham/runav_rate.txt`. The last column in this file should be the most accurate. If you plot this, you get the estimate of the transition rate after each accepted path in your simulation. The last value of this running estimate is the best estimate, as that point uses all of the available data.

This rate is in infretis units, meaning the rate is the number of transitions per unit of infretis time.

$$\Delta t = \text{subcycles }\times \text{timestep}$$

So the rate in correct time units is

$$r = r' / \Delta t$$

where $r'$ is the rate in `runav_rate.txt`.

Below is a plot of the running estimate of this rate, and the corresponding error estimates that are calculated in `errRATE.txt`. If you inspect the first row of this .txt file, you can also see a value of the averaged relative error estimate. This is the horizontal line in the right plot below.

![transition-rate](https://github.com/user-attachments/assets/24a55fca-4e04-4ecb-a9dc-0d5c8ef97152)

#### The crossing probability
Another interesting property we can inspect is the crossing probability in `Pcross.txt` (right figure below), which is the probability of reaching a given value of the order parameter, given that you start in state A. The crossing probability is constructed from the local crossing probabilities in `ploc_unscaled.txt` using either point matching or WHAM. The local crossing probability is the probability that you cross interface $i+1$ given that you are part of the [i+] ensemble. The local crossing probabilities are a useful check that the interfaces are placed adequately. In an ideal case, all interfaces have the same crossing probability. 

![probabilites](https://github.com/user-attachments/assets/7cfff00e-960d-433d-8c15-0732ba721d32)

### Error Estimates

Each property calculated with the WHAM has a corresponding error estimate, given by **the relative standard error of the mean**. For a property $A$ (for example the local crossing probability in ensemble [0+]), the relative standard error of the mean is given by

$\sigma_{\langle A \rangle} = \frac{\sigma_A}{\sqrt{N}\langle A \rangle}$

where $\sigma_A$ is the standard deviation of property A, $\langle A \rangle$ is the average of property A, and $N$ is the number of samples.

For the properties that are true running *averages* (the local crossing proabilities and the path lengths), these error estimates are identical to those from a regular block averaging procedure [More on this can be found here.](https://doi.org/10.1002/jcc.27319)



Suppose you want to sample some property $A$ during an infretis simulation. That could for example be the local crossing probability in ensemble [0+], which is calculated in infretis by sampling a set of paths using a Monte Carlo (MC) procedure. The property $A$ can fluctuate a bit, and since we have a finite number of samples, we only have an estimate of the true underlying average $\bar{A}$. If we now estimate the average $\langle A \rangle$ from our finite number of samples, we would likely be a bit off due to random fluctuations (not taking into account systematic errors, more on this later).



We can say something about how far we are off from the true mean $\bar{A}$ by invoking the central-limit theorem which, under certain weak assumptions and a resonable number of samples, implies that an estimate of the average $\langle A \rangle$ is distributed normally around the true average $\bar{A}$ with variance $\sigma^2_{\langle A \rangle}$ (the variance of the mean). Taking the square root then results in the standard error of the mean, which says something about how far off we may be from the true mean due to the *random sampling process*.



It should be stressed that the error estimates are based on random fluctuations in the variable $A$, and says something about the random sampling process. It does not say anything about systematic errors; when we have long-lasting correlations with lifetimes longer than what we have simulated. For example, if we do not sample a certain reaction channel during the simulation at all, our error estimate (and probably also the mean) would be off.

New paths are generated by modifying an existing path using, in most cases, molecular dynamics (MD). Because we use MD the new paths will usually be rather similar (correlated) to the old path we started from, but as we continue modifying a path, they start to decorrelate. This is an important point. Suppose we want to calculate the variance of the mean of property $A$ given paths 1, 2, ..., M:

$\sigma^2_{\langle A \rangle} = \frac{1}{M}\sigma^2_A$

Since the $M$ paths are correlated, we would underestimate the variance of the mean since the effective number of independent measurements is actually less than $M$. We could take into account the autocorrelation function of property $A$ to correct for this. However, a simpler solution is to use a **block averaging** procedure, which allows us to calculate a more "correct" value for the standard error of the mean taking into account these correlations. But again, it can't correct for correlations with lifetimes longer than what we have sampled.



##### Examples

Below, we apply the error estimation on the running average of a random uniform variable with $N$ samples to see how far off our estimate of the mean is from the true mean of $0.5$.

```python
from inftools.analysis.rec_error import rec_block_errors
import numpy as np
import matplotlib.pyplot as plt

N = 10000
x = np.random.rand(N)
runav_x = np.cumsum(x)/np.arange(1, N+1)

half_average_error, statistical_ineficiency, relative_error = rec_block_errors(runav_x, N/200)
f, (a0, a1) = plt.subplots(2,1)
a0.plot(runav_x, label = "running average")
a0.legend()
a1.plot(relative_error)
a1.axhline(half_average_error,c="C1", label="half_average_relative_error")
a1.axhline(relative_error[0], c="C2", label="relative_error_[0]")
# relative error so divide by the average
a1.axhline(np.std(x)/(0.5*N**runav_x[-1]),c="C3", label = "std(x)/(sqrt(N)*runav_x[-1])")
a1.set(ylim = (half_average_error - half_average_error*0.5,half_average_error + half_average_error*0.5))
a1.legend()
plt.tight_layout()
plt.show()
```

Next, we consider the case where the values are correlated, where we underestimate the relative standard error of the mean compared to the block average error.

```python
corr = 0.8
N = 10000

# from https://stackoverflow.com/a/33904277
mu = 0.5
sigma = 1.0
c = mu * (1 - corr)
sigma_e = np.sqrt((sigma ** 2) * (1 - corr ** 2))
x = [np.random.normal(mu, sigma)]
for _ in range(1, N):
	x.append(c + corr * x[-1] + np.random.normal(0, sigma_e))

runav_x = np.cumsum(x)/np.arange(1, N+1)

half_average_error, statistical_ineficiency, relative_error = rec_block_errors(runav_x, N/200)
plt.plot(relative_error)
plt.axhline(half_average_error,c="C1", label="half_average_relative_error")
plt.axhline(np.std(x)/(N**0.5*runav_x[-1]),c="C3", label = "std(x)/(sqrt(N)*runav_x[-1])")
plt.ylim(0, np.max(relative_error))
plt.legend()
plt.show()
```

