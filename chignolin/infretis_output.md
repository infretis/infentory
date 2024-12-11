### WHAM output

All of the processed simulation results can be found in .txt files in the `wham/` folder. Each of the running estimated properties has corresponding errors that are calculated with a block-averaging procedure.

* runav_rate.txt
  The running estimate of the transition rate, meaning an estimate of the rate after each accepted path. errRATE.txt contains the error estimates.

* Pcross.txt
  The estimate of the total crossing probability as a function of the order parameter.

* ploc_unscaled.txt
  The local crossing probabilities in each of the plus ensembles ([0+], [1+], ...,).

* runav_flux.txt
  The running estimate of the flux through the first interfaces $\lambda_0$. errFLUX.txt contains the errors.

* runav_Ptot.txt
  The running estimate of the total crossing probability. This is the last number of contained in Pcross.txt, estimated after each accepted path. errPtot.txt contains the errors.

* runav_ploc.txt
  The running average of the local crossing probabilities in each ensemble. These are the values where the curves from ploc_unscaled.txt reach the neighboring interface. errploc.txt contains the errors.
- pathdistr.txt
  A distribution of the weighted path lengths in each ensemble.

- pathlengths.txt
  The average path length in each ensemble.
* runav_L0.txt
  The running averages of paths in either [0-] or [0+], which are used for the flux calculation. errL0.txt contains the errors.
- ploc_pointmatch.txt
  The same as ploc_unscaled.txt, except that the values are scaled by point matching, which is used to construct the total crossing probabilities in Pcross.txt.

- ploc_WHAM.txt
  Same as above, but scaled with the WHAM weights, which is more accurate than point-matching.

The first line of these files tells you what the different columns are.

#### The transition rate
The running *estimate* of the transition rate can be found in the file `wham/runav_rate.txt`. The last column in this file should be the most accurate. If you plot this, you get the estimate of the transition rate after each accepted path in your simulation. The last value of this running estimate is the best estimate, as that point uses all of the available data.

This rate is in infretis units, meaning the rate is the number of transitions per unit of infretis time.

$$\Delta t = \text{subcycles }\times \text{timestep}$$

So the rate in correct time units is

$$r = r' / \Delta t$$

where $r'$ is the rate in `runav_rate.txt`.

Below is a plot of the running estimate of this rate, and the corresponding error estimates.
![transition-rate](https://github.com/user-attachments/assets/24a55fca-4e04-4ecb-a9dc-0d5c8ef97152)


#### The crossing probability
Another interesting property we can inspect is the crossing probability (right figure below), which is the probability of reaching a given value of the order parameter, given that you start in state A. The crossing probability is constructed from the local crossing probabilities using either point matching or WHAM. These are the probabilities that you cross interface $i+1$ given that you are part of the [i+] ensemble. The local crossing probabilities are a useful check that the interfaces are placed adequately. In an ideal case, all interfaces have the same crossing probability. 

![probabilites](https://github.com/user-attachments/assets/7cfff00e-960d-433d-8c15-0732ba721d32)
