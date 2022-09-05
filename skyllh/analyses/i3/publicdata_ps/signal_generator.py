# -*- coding: utf-8 -*-

import numpy as np
from scipy import interpolate

from skyllh.core.py import (
    issequenceof,
    float_cast,
    int_cast
)
from skyllh.core.llhratio import LLHRatio
from skyllh.core.dataset import Dataset
from skyllh.core.source_hypothesis import SourceHypoGroupManager
from skyllh.core.storage import DataFieldRecordArray

from skyllh.analyses.i3.publicdata_ps.utils import (
    psi_to_dec_and_ra,
    PublicDataSmearingMatrix,
)
from skyllh.analyses.i3.publicdata_ps.pd_aeff import PDAeff


class PublicDataDatasetSignalGenerator(object):

    def __init__(self, ds, src_dec, cache_effA=None, cache_sm=None, **kwargs):
        """Creates a new instance of the signal generator for generating
        signal events from a specific public data dataset.
        """
        super().__init__(**kwargs)

        if cache_sm is None:
            self.smearing_matrix = PublicDataSmearingMatrix(
                pathfilenames=ds.get_abs_pathfilename_list(
                    ds.get_aux_data_definition('smearing_datafile')))
        else:
            self.smearing_matrix = cache_sm

        if cache_effA is None:
            dec_idx = self.smearing_matrix.get_true_dec_idx(src_dec)
            (min_log_true_e,
             max_log_true_e) = \
                self.smearing_matrix.get_true_log_e_range_with_valid_log_e_pdfs(
                    dec_idx)
            kwargs = {
                'src_dec': src_dec,
                'min_log_e': min_log_true_e,
                'max_log_e': max_log_true_e
            }
            self.effA = PDAeff(
                pathfilenames=ds.get_abs_pathfilename_list(
                    ds.get_aux_data_definition('eff_area_datafile')), **kwargs)

        else:
            self.effA = cache_effA

    def _generate_inv_cdf_spline(self, flux_model, log_e_min,
                                 log_e_max):
        """Sample the true neutrino energy from the power-law
        re-weighted with the detection probability.
        """
        m = (self.effA.log10_enu_bincenters >= log_e_min) & (
            self.effA.log10_enu_bincenters < log_e_max)
        bin_centers = self.effA.log10_enu_bincenters[m]
        low_bin_edges = self.effA._log10_enu_binedges_lower[m]
        high_bin_edges = self.effA._log10_enu_binedges_upper[m]

        # Flux probability P(E_nu | gamma) per bin.
        flux_prob = flux_model.get_integral(
            10**low_bin_edges, 10**high_bin_edges
        ) / flux_model.get_integral(
            10**low_bin_edges[0], 10**high_bin_edges[-1]
        )

        # Do the product and normalize again to a probability per bin.
        product = flux_prob * self.effA.det_prob
        prob_per_bin = product / np.sum(product)

        # The probability per bin cannot be zero, otherwise the cumulative
        # sum would not be increasing monotonically. So we set zero bins to
        # 1000 times smaller than the smallest non-zero bin.
        m = prob_per_bin == 0
        prob_per_bin[m] = np.min(prob_per_bin[np.invert(m)]) / 1000
        prob_per_bin /= np.sum(prob_per_bin)

        # Compute the cumulative distribution CDF.
        cum_per_bin = np.cumsum(prob_per_bin)
        cum_per_bin = np.concatenate(([0], cum_per_bin))
        if np.any(np.diff(cum_per_bin) == 0):
            raise ValueError(
                'The cumulative sum of the true energy probability is not '
                'monotonically increasing! Values of the cumsum are '
                f'{cum_per_bin}.')

        bin_centers = np.concatenate(([low_bin_edges[0]], bin_centers))

        # Build a spline for the inverse CDF.
        return interpolate.splrep(cum_per_bin, bin_centers, k=1, s=0)

    @staticmethod
    def _eval_spline(x, spl):
        values = interpolate.splev(x, spl, ext=3)
        return values

    def _generate_events(
            self, rss, src_dec, src_ra, dec_idx,
            log_true_e_inv_cdf_spl, n_events):
        """Generates `n_events` signal events for the given source location
        and flux model.

        Note:
            Some values can be NaN in cases where a PDF was not available!

        Parameters
        ----------
        rss : instance of RandomStateService
            The instance of RandomStateService to use for drawing random
            numbers.
        src_dec : float
            The declination of the source in radians.
        src_ra : float
            The right-ascention of the source in radians.

        Returns
        -------
        events : numpy record array of size `n_events`
            The numpy record array holding the event data.
            It contains the following data fields:
                - 'isvalid'
                - 'log_true_energy'
                - 'log_energy'
                - 'sin_dec'
            Single values can be NaN in cases where a pdf was not available.
        """

        # Create the output event DataFieldRecordArray.
        out_dtype = [
            ('isvalid', np.bool_),
            ('log_true_energy', np.double),
            ('log_energy', np.double),
            ('dec', np.double),
            ('ra', np.double),
            ('sin_dec', np.double),
            ('ang_err', np.double),
            ('time', int),
            ('azi', np.double),
            ('zen', np.double),
            ('run', int)
        ]

        data = dict(
            [(out_dt[0], np.empty(
                (n_events,),
                dtype=out_dt[1])
              ) for out_dt in out_dtype]
        )

        events = DataFieldRecordArray(data, copy=False)

        sm = self.smearing_matrix

        log_true_e = self._eval_spline(
            rss.random.uniform(size=n_events), log_true_e_inv_cdf_spl)

        events['log_true_energy'] = log_true_e

        log_true_e_idxs = (
            np.digitize(log_true_e, bins=sm.true_e_bin_edges) - 1
        )

        # Sample reconstructed energies given true neutrino energies.
        (log_e_idxs, log_e) = sm.sample_log_e(rss, dec_idx, log_true_e_idxs)
        events['log_energy'] = log_e

        # Sample reconstructed psi values given true neutrino energy and
        # reconstructed energy.
        (psi_idxs, psi) = sm.sample_psi(
            rss, dec_idx, log_true_e_idxs, log_e_idxs)

        # Sample reconstructed ang_err values given true neutrino energy,
        # reconstructed energy, and psi.
        (ang_err_idxs, ang_err) = sm.sample_ang_err(
            rss, dec_idx, log_true_e_idxs, log_e_idxs, psi_idxs)

        isvalid = np.invert(
            np.isnan(log_e) | np.isnan(psi) | np.isnan(ang_err))
        events['isvalid'] = isvalid

        # Convert the psf into a set of (r.a. and dec.). Only use non-nan
        # values.
        (dec, ra) = psi_to_dec_and_ra(rss, src_dec, src_ra, psi[isvalid])
        events['ra'][isvalid] = ra
        events['dec'][isvalid] = dec
        events['sin_dec'][isvalid] = np.sin(dec)

        # Add an angular error. Only use non-nan values.
        events['ang_err'][isvalid] = ang_err[isvalid]

        # Add fields required by the framework
        events['time'] = np.ones(n_events)
        events['azi'] = np.ones(n_events)
        events['zen'] = np.ones(n_events)
        events['run'] = -1 * np.ones(n_events)

        return events

    def generate_signal_events(
            self, rss, src_dec, src_ra, flux_model, n_events):
        """Generates ``n_events`` signal events for the given source location
        and flux model.

        Returns
        -------
        events : numpy record array
            The numpy record array holding the event data.
            It contains the following data fields:
                - 'isvalid'
                - 'log_true_energy'
                - 'log_energy'
                - 'dec'
                - 'ra'
                - 'ang_err'
        """
        sm = self.smearing_matrix

        # Find the declination bin index.
        dec_idx = sm.get_true_dec_idx(src_dec)

        # Determine the true energy range for which log_e PDFs are available.
        (min_log_true_e,
         max_log_true_e) = sm.get_true_log_e_range_with_valid_log_e_pdfs(
             dec_idx)
        # Build the spline for the inverse CDF and draw a true neutrino
        # energy from the hypothesis spectrum.
        log_true_e_inv_cdf_spl = self._generate_inv_cdf_spline(
            flux_model, min_log_true_e, max_log_true_e)

        events = None
        n_evt_generated = 0
        while n_evt_generated != n_events:
            n_evt = n_events - n_evt_generated

            events_ = self._generate_events(
                rss, src_dec, src_ra, dec_idx, log_true_e_inv_cdf_spl, n_evt)

            # Cut events that failed to be generated due to missing PDFs.
            events_ = events_[events_['isvalid']]
            if not len(events_) == 0:
                n_evt_generated += len(events_)
                if events is None:
                    events = events_
                else:
                    events.append(events_)

        return events


class PublicDataSignalGenerator(object):
    """This class provides a signal generation method for a point-like source
    seen in the IceCube detector using the 10 years public data release.
    """

    def __init__(self, src_hypo_group_manager, dataset_list, data_list=None, llhratio=None):
        self.src_hypo_group_manager = src_hypo_group_manager
        self.dataset_list = dataset_list
        self.data_list = data_list
        self.llhratio = llhratio
        self.cache_effA = list()
        self.cache_sm = list()
        [self.cache_effA.append(None) for i in self._dataset_list]
        [self.cache_sm.append(None) for i in self._dataset_list]

    @property
    def src_hypo_group_manager(self):
        """The SourceHypoGroupManager instance defining the source groups with
        their spectra.
        """
        return self._src_hypo_group_manager

    @src_hypo_group_manager.setter
    def src_hypo_group_manager(self, manager):
        if(not isinstance(manager, SourceHypoGroupManager)):
            raise TypeError('The src_hypo_group_manager property must be an '
                            'instance of SourceHypoGroupManager!')
        self._src_hypo_group_manager = manager

    @property
    def dataset_list(self):
        """The list of Dataset instances for which signal events should get
        generated for.
        """
        return self._dataset_list

    @dataset_list.setter
    def dataset_list(self, datasets):
        if(not issequenceof(datasets, Dataset)):
            raise TypeError('The dataset_list property must be a sequence of '
                            'Dataset instances!')
        self._dataset_list = list(datasets)

    @property
    def llhratio(self):
        """The log-likelihood ratio function for the analysis.
        """
        return self._llhratio

    @llhratio.setter
    def llhratio(self, llhratio):
        if llhratio is not None:
            if(not isinstance(llhratio, LLHRatio)):
                raise TypeError('The llratio property must be an instance of '
                                'LLHRatio!')
        self._llhratio = llhratio

    def generate_signal_events(self, rss, mean, poisson=True):
        shg_list = self._src_hypo_group_manager.src_hypo_group_list

        tot_n_events = 0
        signal_events_dict = {}

        for shg in shg_list:
            # This only works with power-laws for now.
            # Each source hypo group can have a different power-law
            gamma = shg.fluxmodel.gamma
            weights, _ = self.llhratio.dataset_signal_weights([mean, gamma])
            for (ds_idx, w) in enumerate(weights):
                w_mean = mean * w
                if(poisson):
                    n_events = rss.random.poisson(
                        float_cast(
                            w_mean,
                            '`mean` must be castable to type of float!'
                        )
                    )
                else:
                    n_events = int_cast(
                        w_mean,
                        '`mean` must be castable to type of int!'
                    )
                tot_n_events += n_events

                events_ = None
                for (shg_src_idx, src) in enumerate(shg.source_list):
                    ds = self._dataset_list[ds_idx]
                    sig_gen = PublicDataDatasetSignalGenerator(
                        ds, src.dec, self.cache_effA[ds_idx],
                        self.cache_sm[ds_idx])
                    if self.cache_effA[ds_idx] is None:
                        self.cache_effA[ds_idx] = sig_gen.effA
                    if self.cache_sm[ds_idx] is None:
                        self.cache_sm[ds_idx] = sig_gen.smearing_matrix
                    # ToDo: here n_events should be split according to some
                    # source weight
                    events_ = sig_gen.generate_signal_events(
                        rss,
                        src.dec,
                        src.ra,
                        shg.fluxmodel,
                        n_events
                    )
                    if events_ is None:
                        continue

                    if shg_src_idx == 0:
                        signal_events_dict[ds_idx] = events_
                    else:
                        signal_events_dict[ds_idx].append(events_)

        return tot_n_events, signal_events_dict
