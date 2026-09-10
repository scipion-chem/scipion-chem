"""
idw.py
======
Compact-support Inverse Distance Weighting interpolation of 3D displacement
vectors using a KD-Tree.

Weight formula:

    w_k(x) = ( max(0, R - d(x, x_k)) / (R * d(x, x_k)) ) ^ p

Interpolated displacement at query point x:

    u(x) = sum_i( w_i(x) * u_i ) / sum_i( w_i(x) )   if d(x, x_i) != 0 for all i
    u(x) = u_i                                          if d(x, x_i) = 0 for some i

where u_i = dst_CA[i] - src_CA[i]  (displacement vector at Ca i)

Properties:
  - w_k = 0 exactly when d >= R  (compact support, no discontinuity at boundary)
  - w_k -> 0 smoothly as d -> R
  - d = 0 handled explicitly (exact-hit special case, no division by zero)

Reference:
  invdisttree.py - Denis / stackoverflow.com/a/3119544
  gist.github.com/pyRobShrk/11c45b8d6aad6d1d5ecd7869f502d8f0
"""

import numpy as np
from scipy.spatial import cKDTree


class InvDistTree3D:
    """
    Parameters
    ----------
    src_ca : ndarray (N, 3)
        Source Ca positions (control points).
    leafsize : int
        KD-Tree leaf size (default 10).

    Usage
    -----
    idw = InvDistTree3D(src_ca)
    interpolated = idw(all_atom_coords, displacements, R=15.0, k=8)
    """

    def __init__(self, src_ca, leafsize=10):
        src_ca = np.asarray(src_ca, dtype=np.float64)
        self.tree = cKDTree(src_ca, leafsize=leafsize)

    def __call__(self, q, displacements, R=15.0, k=8, p=2.0):
        """
        Interpolate displacement vectors at query points.

        Parameters
        q : ndarray (M, 3)
            Query positions (atom coordinates).
        R : float
            Search radius in Angstroms. Only Ca within this sphere contribute.
            Weight is exactly 0 at d >= R.
        k : int
            Maximum number of nearest Ca to use within R.
            If fewer than k Ca are found within R, all of them are used.

        Returns
        ndarray (M, 3) - interpolated displacement vectors.
                         Add to q to get the displaced positions:
                         new_coords = q + idw(q, R, k)
        """
        if displacements.shape != self.tree.data.shape:
            raise ValueError(
                'displacements shape %s does not match src_ca shape %s'
                % (displacements.shape, self.tree.data.shape)
            )

        displacements = np.asarray(displacements, dtype=np.float64)
        q = np.asarray(q, dtype=np.float64)
        single = q.ndim == 1

        if single:
            q = q[np.newaxis, :]  # create matrix, no vector (from (4, ) to (1, 4))

        result = np.zeros_like(q)
        neighbours = self.tree.query_ball_point(q, r=R, workers=-1)

        for j in range(len(q)):
            idx = np.array(neighbours[j], dtype=np.intp)

            if len(idx) == 0:
                # No Ca within R -> zero displacement, atom stays in place
                continue

            dists = np.linalg.norm(self.tree.data[idx] - q[j], axis=1)

            # Keep only the k nearest within R
            if len(idx) > k:
                order = np.argsort(dists)[:k]
                idx = idx[order]
                dists = dists[order]

            # if d(x, x_i) = 0 for some i -> u(x) = u_i  (no division by zero)
            zero_mask = dists < 1e-10
            if np.any(zero_mask):
                result[j] = displacements[idx[np.argmax(zero_mask)]]
                continue

            # All neighbours have d < R (guaranteed by query_ball_point)
            w = ((R - dists) / (R * dists)) ** p
            w_sum = w.sum()
            if w_sum == 0.0:
                continue

            w /= w_sum
            result[j] = np.dot(w, displacements[idx])  # sumatorio

        return result[0] if single else result