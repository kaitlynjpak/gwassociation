"""Radial distance-overlap terms for two three-dimensional skymaps.

The :class:`RadialOverlap` calculator evaluates the line-of-sight component of a
coincident-localization Bayes factor when both events provide HEALPix distance
metadata.  If either map lacks usable distance information, the calculator
returns a neutral factor so spatial-only workflows continue to run.
"""

import numpy as np

from ..stats import prior_dl2, DL_MIN_MPC, DL_MAX_MPC


class RadialOverlap:
    """
    Compute radial (distance) overlap between two skymaps with distance information.
    Implements the line-of-sight radial Bayes factor described in
    "Coincident Detection Significance in Multimessenger Astronomy".
    """

    def __init__(self, dmin: float = DL_MIN_MPC, dmax: float = DL_MAX_MPC):
        self.dmin = dmin
        self.dmax = dmax

    def compute_map_overlap(self, primary_event, secondary_event) -> float:
        """
        Calculate the radial overlap integral between two events that both have 3D skymaps.
        """
        if not getattr(primary_event, "is_3d", False) or not getattr(secondary_event, "is_3d", False):
            # Distance information is unavailable for one of the skymaps
            return 1.0

        if primary_event.nside != secondary_event.nside:
            raise ValueError("Skymap NSIDE mismatch for radial overlap calculation.")

        if primary_event.skymap is None:
            primary_event.load_skymap()
        if secondary_event.skymap is None:
            secondary_event.load_skymap()

        prob_primary = self._normalize(self._event_prob(primary_event))
        prob_secondary = self._normalize(self._event_prob(secondary_event))

        try:
            mu1, sigma1, _norm1 = self._extract_distance_parameters(primary_event)
            mu2, sigma2, _norm2 = self._extract_distance_parameters(secondary_event)
        except ValueError:
            return 1.0

        # Work with the *physical* conditional-distance moments (mean, std), not
        # the raw ansatz parameters distmu/distsigma. In low-probability pixels
        # the raw parameters are pathological (distmu can be large-negative and
        # distnorm ~ 1e80); left unmasked they detonate the weighted sum. The
        # moments are finite and well-behaved wherever the pixel carries usable
        # distance information, and non-finite/unphysical elsewhere.
        mean1, std1 = self._conditional_moments(mu1, sigma1)
        mean2, std2 = self._conditional_moments(mu2, sigma2)

        valid = (
            np.isfinite(mean1) & np.isfinite(mean2)
            & np.isfinite(std1) & np.isfinite(std2)
            & (std1 > 0) & (std2 > 0)
            & (mean1 > 0) & (mean2 > 0)
            & (mean1 < self.dmax) & (mean2 < self.dmax)
        )

        weights = np.where(valid, prob_primary * prob_secondary, 0.0)
        weight_sum = weights.sum()
        if weight_sum <= 0:
            # No sky pixel carries usable, overlapping distance information;
            # the distance term is uninformative rather than a rejection.
            return 1.0
        weights = weights / weight_sum

        # Per line of sight, the two conditional distance PDFs are approximated
        # as Gaussians in luminosity distance. Their overlap Bayes factor is
        #   I_dL(n) = integral p1(d) p2(d) / prior(d) dd
        #          ~= N(mean1 - mean2; 0, sqrt(std1^2 + std2^2)) / prior(d_p),
        # with prior(d) proportional to d^2 (uniform in volume) evaluated at the
        # inverse-variance-weighted distance d_p where the product concentrates.
        # Substitute safe placeholders in invalid pixels so the vectorised math
        # never evaluates inf/nan (those pixels are zeroed out by ``valid``).
        mean1_s = np.where(valid, mean1, 0.0)
        mean2_s = np.where(valid, mean2, 0.0)
        var_safe = np.where(valid, std1 ** 2 + std2 ** 2, 1.0)
        overlap = np.where(
            valid,
            1.0 / np.sqrt(2.0 * np.pi * var_safe)
            * np.exp(-0.5 * (mean1_s - mean2_s) ** 2 / var_safe),
            0.0,
        )

        inv1 = 1.0 / np.where(valid, std1 ** 2, 1.0)
        inv2 = 1.0 / np.where(valid, std2 ** 2, 1.0)
        mean_p = np.where(valid, (mean1_s * inv1 + mean2_s * inv2) / (inv1 + inv2), 1.0)
        prior = np.clip(prior_dl2(mean_p, self.dmin, self.dmax), 1e-30, None)

        radial_factor = np.where(valid, overlap / prior, 0.0)
        return float(np.sum(weights * radial_factor))

    @staticmethod
    def _conditional_moments(distmu, distsigma):
        """Return the physical (mean, std) of the conditional distance PDF.

        Uses ``ligo.skymap.distance.parameters_to_moments`` to convert the
        ansatz parameters into the mean and standard deviation of
        ``p(d) proportional to N(d; distmu, distsigma) d^2`` -- always finite and
        positive where the pixel has usable distance information. Falls back to
        treating the ansatz parameters as approximate moments if
        ``ligo.skymap`` is unavailable.
        """
        distmu = np.asarray(distmu, dtype=float)
        distsigma = np.asarray(distsigma, dtype=float)
        try:
            import ligo.skymap.distance as lsd

            mean, std, _ = lsd.parameters_to_moments(distmu, distsigma)
            return np.asarray(mean, dtype=float), np.asarray(std, dtype=float)
        except Exception:
            return distmu, distsigma

    @staticmethod
    def _event_prob(event):
        if hasattr(event, "prob") and event.prob is not None:
            return event.prob
        if isinstance(getattr(event, "skymap", None), dict):
            return event.skymap.get("prob", event.skymap.get("data"))
        return event.skymap

    @staticmethod
    def _normalize(probabilities):
        arr = np.asarray(probabilities, dtype=float)
        total = np.sum(arr)
        if total <= 0:
            raise ValueError("Skymap probabilities must sum to a positive value.")
        return arr / total

    @staticmethod
    def _extract_distance_parameters(event):
        distances = getattr(event, "distances", None)
        if not distances:
            raise ValueError("No distance information available.")

        mu = distances.get("distmu")
        if mu is None:
            mu = distances.get("distmean") or distances.get("mean")

        sigma = distances.get("distsigma")
        if sigma is None:
            sigma = distances.get("diststd") or distances.get("std")
        norm = distances.get("distnorm")

        if mu is None or sigma is None or norm is None:
            raise ValueError("Distance maps must contain distmu/distsigma/distnorm.")

        mu = np.asarray(mu, dtype=float)
        sigma = np.asarray(sigma, dtype=float)
        norm = np.asarray(norm, dtype=float)

        return mu, sigma, norm
