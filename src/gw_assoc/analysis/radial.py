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

        prob_primary = self._normalize(primary_event.skymap)
        prob_secondary = self._normalize(secondary_event.skymap)

        joint = prob_primary * prob_secondary
        joint_sum = np.sum(joint)
        if joint_sum <= 0:
            return 0.0
        weights = joint / joint_sum

        try:
            mu1, sigma1, norm1 = self._extract_distance_parameters(primary_event)
            mu2, sigma2, norm2 = self._extract_distance_parameters(secondary_event)
        except ValueError:
            return 1.0

        variance = sigma1 ** 2 + sigma2 ** 2 + 1e-12
        gaussian_overlap = (
            1.0 / np.sqrt(2.0 * np.pi * variance)
        ) * np.exp(-0.5 * ((mu1 - mu2) ** 2) / variance)

        los_overlap = gaussian_overlap * norm1 * norm2
        mu_bar = 0.5 * (mu1 + mu2)
        prior = prior_dl2(mu_bar, self.dmin, self.dmax)
        prior = np.clip(prior, 1e-30, None)

        radial_factor = los_overlap / prior
        return float(np.sum(weights * radial_factor))

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

