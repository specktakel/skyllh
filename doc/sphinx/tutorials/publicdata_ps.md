---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.14.1
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

Working with the public 10-year IceCube point-source data
==


This tutorial shows how to use the IceCube public 10-year point-source data with SkyLLH.

<!-- #raw raw_mimetype="text/restructuredtext" -->
**Disclaimer**

    The released 10-year IceCube point-source data can reproduce the published results only within
    a certain amount of uncertainty due to the limited instrument response function binning 
    provided in the data release. The IceCube collaboration is able to reproduce the published 
    results using detailed direct simulation data, as done for the publication.
<!-- #endraw -->

<!-- #raw -->

<!-- #endraw -->

```python
import numpy as np
from matplotlib import pyplot as plt
import scipy.stats
```

```python
from hierarchical_nu.events import Events
from hierarchical_nu.utils.roi import NorthernSkyROI
from hierarchical_nu.source.parameter import Parameter
from hierarchical_nu.detector.icecube import IC40, IC59, IC79, IC86_I, IC86_II
from astropy import units as u
from pathlib import Path
```

```python
skyllh_path = Path("/Users/David/Documents/phd/icecube/skyllh")
```

```python
Emin_det = Parameter(3e2*u.GeV, "Emin_det", fixed=True)
```

```python
NorthernSkyROI()
```

```python
events = Events.from_ev_file(IC40, IC59, IC79, IC86_I, IC86_II)
```

```python
events.export_to_csv(skyllh_path / "icecube_10year_ps/events")
```

Creating a configuration
---


First we have to create a configuration instance.

```python
from skyllh.core.config import Config
cfg = Config()
```

```python
cfg
```

Getting the datasets
---

```python
cfg
```

Now we import the dataset definition of the public 10-year point-source data set:

```python
grl = []
for p in ["IC86_II", "IC86_III", "IC86_IV", "IC86_V", "IC86_VI", "IC86_VII"]:
    grl.append(np.loadtxt(skyllh_path / f"icecube_10year_ps/uptime/{p}_exp.csv"))
```

```python
np.savetxt(skyllh_path / f"icecube_10year_ps/uptime/IC86_II_exp.csv", np.vstack(grl))
```

```python
from skyllh.datasets.i3.PublicData_10y_ps import create_dataset_collection
```

The collection of datasets can be created using the ``create_dataset_collection`` function. This function requires the base path to the data repository. It's the path where the public point-source data is stored. The public point-source data can be downloaded from the [IceCube website](http://icecube.wisc.edu/data-releases/20210126_PS-IC40-IC86_VII.zip).

```python
dsc = create_dataset_collection(
    cfg=cfg, 
    base_path=skyllh_path)
```

The ``dataset_names`` property provides a list of all the data sets defined in the data set collection of the public point-source data.

```python
dsc.dataset_names
```

The individual data sets ``IC86_II``, ``IC86_III``, ``IC86_IV``, ``IC86_V``, ``IC86_VI``, and ``IC86_VII`` are also available as a single combined data set ``IC86_II-VII``, because these data sets share the same detector simulation and event selection. Hence, we can get a list of data sets via the access operator ``[dataset1, dataset2, ...]`` of the ``dsc`` instance:

```python
cfg
```

```python
datasets = dsc['IC40', 'IC59', 'IC79', 'IC86_I', 'IC86_II']
```

```python
datasets[0].datafields
```

Getting the analysis
---


The analysis used for the published PRL results is referred in SkyLLH as "*traditional point-source analysis*" and is pre-defined:

```python
from skyllh.analyses.i3.publicdata_ps.time_integrated_ps import create_analysis
```

```python
help(create_analysis)
```

As source we use TXS 0506+056.

```python
from skyllh.core.source_model import PointLikeSource
```

```python
source = PointLikeSource(ra=np.deg2rad(77.35), dec=np.deg2rad(5.7))
```

```python
ana = create_analysis(cfg=cfg, datasets=datasets, source=source)
```

Initializing a trial
---


After the `Analysis` instance was created trials can be run. To do so the analysis needs to be initialized with some trial data. For instance we could initialize the analysis with the experimental data to "unblind" the analysis afterwards. Technically the `TrialDataManager` of each log-likelihood ratio function, i.e. dataset, is initialized with data.

The `Analysis` class provides the method `initialize_trial` to initialize a trial with data. It takes a list of `DataFieldRecordArray` instances holding the events. If we want to initialize a trial with the experimental data, we can get that list from the `Analysis` instance itself:

```python
events_list = [ data.exp for data in ana.data_list ]
ana.initialize_trial(events_list)
```

```python
np.deg2rad(-5)
```

Maximizing the log-likelihood ratio function
---


After initializing a trial, we can maximize the LLH ratio function using the `maximize_llhratio` method of the `Analysis` class. This method requires a ``RandomStateService`` instance in case the minimizer does not succeed and a new set of initial values for the fit parameters need to get generated. The method returns a 4-element tuple. The first element is the set of fit parameters used in the maximization. The second element is the value of the LLH ration function at its maximum. The third element is the array of the fit parameter values at the maximum, and the forth element is the status dictionary of the minimizer.

```python
from skyllh.core.random import RandomStateService
rss = RandomStateService(seed=1)
```

```python
(log_lambda_max, fitparam_values, status) = ana.llhratio.maximize(rss)
```

```python
print(f'log_lambda_max = {log_lambda_max}')
print(f'fitparam_values = {fitparam_values}')
print(f'status = {status}')
```

Calculating the test-statistic
---


Using the maximum of the LLH ratio function and the fit parameter values at the maximum we can calculate the test-statistic using the `calculate_test_statistic` method of the `Analysis` class:

```python
TS = ana.calculate_test_statistic(log_lambda_max, fitparam_values)
print(f'TS = {TS:.3f}')
```

## Unblinding the data


After creating the analysis instance we can unblind the data for the choosen source. Hence, we initialize the analysis with a trial of the experimental data, maximize the log-likelihood ratio function for all given experimental data events, and calculate the test-statistic value. The analysis instance has the method ``unblind`` that can be used for that. This method requires a ``RandomStateService`` instance in case the minimizer does not succeed and a new set of initial values for the fit parameters need to get generated.

```python
from skyllh.core.random import RandomStateService
rss = RandomStateService(seed=1)
```

```python
help(ana.unblind)
```

The ``unblind`` method returns the test-statistic value, the best-fit fit parameter values, and a status dictionary of the minimizer.

```python
(ts, x, status) = ana.unblind(minimizer_rss=rss)
```

```python
print(f'TS = {ts:.3f}')
print(f'ns = {x["ns"]:.2f}')
print(f'gamma = {x["gamma"]:.2f}')
```

## Calculating the corresponding flux normalization 


By default the analysis is created with a flux normalization of 1 GeV$^{-1}$s$^{-1}$cm$^{-2}$sr$^{-1}$ (see `refplflux_Phi0` argument of the `create_analysis` method). The analysis instance has the method `calculate_fluxmodel_scaling_factor` that calculates the scaling factor the reference flux normalization has to be multiplied with to represent a given analysis result, i.e. $n_{\text{s}}$ and $\gamma$ value. This function takes the detected mean $n_{\text{s}}$ value as first argument and the list of source parameter values as second argument:

```python
scaling_factor = ana.calculate_fluxmodel_scaling_factor(x['ns'], [x['ns'], x['gamma']])
print(f'Flux scaling factor = {scaling_factor:.3e}')
```

Hence, our result corresponds to a power-law flux of:

```python
print(f'{scaling_factor:.3e}'' (E/1000 GeV)^{-'f'{x["gamma"]:.2f}'+'} 1/(GeV s cm^2 sr)')
```

Evaluating the log-likelihood ratio function
---


Sometimes it is useful to be able to evaluate the log-likelihood ratio function, e.g. for creating a likelihood contour plot. Because SkyLLH's structure is based on the mathematical structure of the likelihood function, the `Analysis` instance has the property `llhratio` which is the class instance of the used log-likelihood ratio function. This instance has the method `evaluate`. The method takes an array of the fit parameter values as argument at which the LLH ratio function will be evaluated. It returns the value of the LLH ratio function at the given point and its gradients w.r.t. the fit parameters.

In our case this is the number of signal events, $n_{\mathrm{s}}$ and the spectral index $\gamma$. If we evaluate the LLH ratio function at the maximum, the gradients should be close to zero.

```python
help(ana.llhratio.evaluate)
```

```python
(llhratio_value, (grad_ns, grad_gamma)) = ana.llhratio.evaluate([x["ns"], x["gamma"]])
print(f'llhratio_value = {llhratio_value:.3f}')
print(f'grad_ns = {grad_ns:.3f}')
print(f'grad_gamma = {grad_gamma:.3f}')
```

Using the `evaluate` method of the `LLHRatio` class we can scan the log-likelihood ratio space and create a contour plot showing the best fit and the 68%, 90%, and 95% quantile assuming Wilks-theorem.

```python
(ns_min, ns_max, ns_step) = (0, 80, 0.5)
(gamma_min, gamma_max, gamma_step) = (1.5, 4.0, 0.1)

ns_edges = np.linspace(ns_min, ns_max, int((ns_max-ns_min)/ns_step)+1)
ns_vals = 0.5*(ns_edges[1:] + ns_edges[:-1])

gamma_edges = np.linspace(gamma_min, gamma_max, int((gamma_max-gamma_min)/gamma_step+1))
gamma_vals = 0.5*(gamma_edges[1:] + gamma_edges[:-1])

delta_ts = np.empty((len(ns_vals), len(gamma_vals)), dtype=np.double)
for (ns_i, ns) in enumerate(ns_vals):
    for (gamma_i, gamma) in enumerate(gamma_vals):

        delta_ts[ns_i, gamma_i] = (
            ana.calculate_test_statistic(llhratio_value, [x["ns"], x["gamma"]]) -
            ana.calculate_test_statistic(ana.llhratio.evaluate([ns, gamma])[0], [ns, gamma])
        )

# Determine the best fit ns and gamma values from the scan.
index_max = np.argmin(delta_ts)
ns_i_max = int(index_max / len(gamma_vals))
gamma_i_max = index_max % len(gamma_vals)
ns_best = ns_vals[ns_i_max]
gamma_best = gamma_vals[gamma_i_max]
```

```python
# Determine the delta lambda value for the 95% quantile assuming a chi-sqaure
# distribution with 2 degrees of freedom (i.e. assuming Wilks theorem).
chi2_68_quantile = scipy.stats.chi2.ppf(0.68, df=2)
chi2_90_quantile = scipy.stats.chi2.ppf(0.90, df=2)
chi2_95_quantile = scipy.stats.chi2.ppf(0.95, df=2)
```

```python
from matplotlib.colors import LogNorm
plt.figure(figsize=(8,6))
plt.pcolormesh(gamma_edges, ns_edges, delta_ts, cmap='nipy_spectral')
cbar = plt.colorbar()
cbar.set_label(r'$\Delta$TS')
plt.contour(gamma_vals, ns_vals, delta_ts, [chi2_68_quantile], colors='#FFFFFF')
plt.contour(gamma_vals, ns_vals, delta_ts, [chi2_90_quantile], colors='#AAAAAA')
plt.contour(gamma_vals, ns_vals, delta_ts, [chi2_95_quantile], colors='#444444')
plt.plot(gamma_best, ns_best, marker='x', color='white', ms=10)
plt.xlabel(r'$\gamma$')
plt.ylabel(r'$n_{\mathrm{s}}$')
plt.ylim(ns_min, ns_max)
plt.xlim(gamma_min, gamma_max)
```

Calculating the significance (local p-value)
---


The significance of the source, i.e. the local p-value, can be calculated by generating the test-statistic distribution of background-only data trials, i.e. for zero injected signal events. SkyLLH provides the helper function ``create_trial_data_file`` to do that:

```python
from skyllh.core.utils.analysis import create_trial_data_file
```

```python
help(create_trial_data_file)
```

At first we will generate 10k trials and look at the test-statistic distribution. We will time the trial generation using the ``TimeLord`` class.

```python
from skyllh.core.timing import TimeLord
tl = TimeLord()
```

```python
rss = RandomStateService(seed=1)
(_, _, _, trials) = create_trial_data_file(
    ana=ana,
    rss=rss,
    n_trials=1e4,
    mean_n_sig=0,
    pathfilename=skyllh_path / 'icecube_10year_ps/txs_bkg_trails.npy',
    ncpu=4,
    tl=tl)
print(tl)
```

After generating the background trials, we can histogram the test-statistic values and plot the TS distribution.

```python
(h, be) = np.histogram(trials['ts'], bins=np.arange(0, np.max(trials['ts'])+0.1, 0.1))
plt.plot(0.5*(be[:-1]+be[1:]), h, drawstyle='steps-mid', label='background')
plt.vlines(ts, 1, np.max(h), label=f'TS(TXS 0506+056)={ts:.3f}')
plt.yscale('log')
plt.xlabel('TS')
plt.ylabel('#trials per bin')
plt.legend()
pass
```

We can see that the TS value of the unblinded data for TXS is rather large and 10k trials are not enough to calculate a reliable estimate for the p-value. Hence, we will generate a few more trials. SkyLLH provides also a helper function to extend the trial data file we just created. It is called ``extend_trial_data_file``: 

```python
from skyllh.core.utils.analysis import extend_trial_data_file
```

```python
help(extend_trial_data_file)
```

```python
tl = TimeLord()
rss = RandomStateService(seed=2)
trials = extend_trial_data_file(
    ana=ana,
    rss=rss,
    n_trials=4e4,
    trial_data=trials,
    pathfilename='/home/mwolf/projects/publicdata_ps/txs_bkg_trails.npy',
    ncpu=8,
    tl=tl)
```

```python
print(tl)
```

The local p-value is defined as the fraction of background trials with TS value greater than the unblinded TS value of the source. 

```python
minus_log10_pval = -np.log10(len(trials[trials['ts'] > ts]) / len(trials))
print(f'-log10(p_local) = {minus_log10_pval:.2f}')
```

```python
(h, be) = np.histogram(trials['ts'], bins=np.arange(0, np.max(trials['ts'])+0.1, 0.1))
plt.plot(0.5*(be[:-1]+be[1:]), h, drawstyle='steps-mid', label='background')
plt.vlines(ts, 1, np.max(h), label=f'TS(TXS 0506+056)={ts:.3f}')
plt.yscale('log')
plt.xlabel('TS')
plt.ylabel('#trials per bin')
plt.legend()
pass
```
