"""trains ML algorithms for simple LLP detection"""

import numpy as np
import tensorflow as tf
from tensorflow import keras
from typing import Tuple, Iterable, Union
from sklearn.preprocessing import StandardScaler


def read_data(fpath: str) -> Tuple[np.ndarray, np.ndarray]:
    """Reads the data produced by my 2 hit skim stored in the npy"""
    pass


def augment_data(x: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray, StandardScaler]:
    pass


def train_bdt(x_trn: np.ndarray, y_trn: np.ndarray, x_val: np.ndarray = None, y_val: np.ndarray = None, **kwargs):
    pass


def train_nn(x_trn: np.ndarray, y_trn: np.ndarray, x_val: np.ndarray = None, y_val: np.ndarray = None, **kwargs):
    pass


if __name__ == "__main__":
    print("Training Classifiers")
    # Cuts on data before training
    # No leptons
    # No jets (if doing cluster + cluster)
    # No jets matched to clusters (if doing jet + cluster)
    # Cluster in time and low spread
    # Nstation == 1?

    # Info cols (DONT TRAIN ON):
    # runNum
    # evtNum
    # cluster number

    # columns:
    # nCscRechitClusters
    # nDtRechitClusters
    # MET
    # MET Phi
    # double these for each cluster
    # CSC or DT
    # NStation10
    # AvgStation10
    # ClusterPhi
    # ClusterEta
    # ClusterMet_dPhi
    # ClusterSize
    # if doing jet + cluster include
    # jetE
    # jetPt
    # jetPhi
    # jetEta

    # target
    # LLP match
    # LLP match E or Pt
