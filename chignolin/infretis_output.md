### WHAM output

All of the processed simulation results can be found in .txt files in the `wham/` folder. Each of the running estimated properties has corresponding errors that are calculated with a block-averaging procedure.

* runav_rate.txt
  The running estimate of the transition rate, meaning an estimate of the rate after each accepted path. errRATE.txt contains the error estimates.

* Pcross.txt
  The estimate of the total crossing probability as a function of the order parameter.

* ploc_unscaled.txt
  The local crossing probabilities in each of the plus ensembles ([0+], [1+], ...,).

* runav_flux.txt
  The running estimate of the flux trough the first interfaces $\lambda_0$. errFLUX.txt contains the errors.

* runav_Ptot.txt
  The running estimate of the total crossing probability. This is the last number of contained in Pcross.txt, estimated after each accepted path. errPtot contains the errors.

* runav_ploc.txt
  The running average of the local crossing probabilites in each ensemble. These are the values where the curves from ploc_unscaled.txt reach the neighboring interface. errploc.txt contains the errors.
- pathdistr.txt
  A distribution of the weighted pathlengths in each ensemble.

- pathlengths.txt
  The average path-length in each ensemble.
* runav_L0.txt
  The running averages of paths in either [0-] or [0+], which are used for the flux calculation. errL0.txt contains the errors.
- ploc_pointmatch.txt
  The same as ploc_unscaled.txt, except that the values are scaled by point matching, which is used to construct the total crossing proabilities in Pcross.txt.

- ploc_WHAM.txt
  Same as above, but scaled with the WHAM weights, which is more accurate than point-macthing.

The first line of these files tell you what the different columns are.

#### The transition rate

The running *estimate* of the transition rate can be found in the file `wham/runav_rate.txt`. The third column in this files should be the most accurate. If you plot this, you get the estimate of the transition rate after each accepted path in your simulation. The last value of this running estimate is the best estimate, as that point uses all of the available data.

These are in infretis units, meaning the rate is the number of transitions per unit of infretis time 

$$\Delta t = \text{subcycles }\times \text{timestep}$$

So the rate in correct time units is

$$r = r' / \Delta t$$

where $r'$ is the rate in `runav_rate.txt`.
