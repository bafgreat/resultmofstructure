#!/usr/bin/python
from __future__ import print_function
__author__ = "Dr. Dinga Wonanke"
__status__ = "production"
import math
import numpy as np


def align_to_ref(ref_coord, coord_to_align):
    """
    Function that aligns a coordinate to a reference in order to minimise the RMSD
    Parameters:
    -----------
    ref_coord: reference coordinate 
    coord_to_align : coordinate to be aligned to the 
                     reference coordinate

    Returns
    -------
    aligned_coords 
    centered_ref
    """

    # find the center of each coordinates
    center_0 = np.mean(ref_coord, axis=0)
    center_1 = np.mean(coord_to_align, axis=0)

    # Center each coordinate
    centered_ref = ref_coord - center_0
    coord_to_com = coord_to_align - center_1

    # Find the rotation that will align coord to Ref
    # Computation of the covariance matrix
    m_metrix = np.dot(np.transpose(coord_to_com), centered_ref)
    # M = coord_1.transpose().dot(coord_0)

    # computing the SVD of the covariance matrix
    u_, s_, v_h = np.linalg.svd(m_metrix)

    # Decide whether we need to correct our rotation matrix to ensure a right-handed coordinate system
    d = (np.linalg.det(u_) * np.linalg.det(v_h)) < 0.0
    if d:
        s_[-1] = -s_[-1]
        u_[:, -1] = -u_[:, -1]

    # Compute rotation matrix
    rotation_matrix = np.dot(u_, v_h)
    # Align the two sets of coordinates and calculate the RMSD
    aligned_coords = coord_to_com.dot(rotation_matrix)

    return aligned_coords, centered_ref


def compute_rmsd(ref_coord, coord_to_align):
    """
    Function that computes the RMSD between two coordinates
    Parameters:
    -----------
    ref_coord: reference coordinate 
    coord_to_align : coordinate to be aligned to the 
                     reference coordinate
    Returns
    -------
    rmsd
    """
    aligned_coords, centered_ref = align_to_ref(ref_coord, coord_to_align)

    distance = list(map(lambda i, j: np.linalg.norm(i-j)
                    ** 2, aligned_coords, centered_ref))

    rmsd = math.sqrt(sum(distance)/float(len(centered_ref)))
    return rmsd
