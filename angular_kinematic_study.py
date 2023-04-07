"""Study for developing cuts on angular variables"""

import os
import pathlib
from collections import defaultdict

import ROOT as rt
import numpy as np

from src import CMS_lumi, tdrstyle
from src.muon_system import MuonSystem
from src.helper_functions import lnot, land, lor, lxor, asum, aabs

from src.muon_system import get_lat_leg, H1D, H2D, multi_plot, make_cluster_eff_1D, draw_csc_z_boxes
from src.helper_functions import canvas

#-----------------------#
# Sci-kit Learn Imports #

import sklearn as skl

from sklearn.preprocessing import (StandardScaler, MinMaxScaler, RobustScaler)
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline, Pipeline

from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis

from sklearn.metrics import (roc_curve, roc_auc_score)

#-----------------------#

from primary_analysis import bins # Eventually move bins to a file in src?

if __name__ == '__main__':

    cur_dir = os.getcwd()
    in_data_dir = cur_dir + "/data/processed/"
    dstat = "ca_0p6"

    fname_vars_mc = in_data_dir + "dRdEtadPhi_mc_" + stat + ".npy"
    fname_vars_r3 = in_data_dir + "dRdEtadPhi_r3_" + stat + ".npy"

    vars_mc, vars_r3 = np.load(fname_var_mc), np.load(fname_var_r3)

    # Load Data #
    X = np.r_[vars_mc, vars_r3]
    y = np.r_[np.ones((vars_mc.shape[0], 1)), np.zeros((vars_r3.shape[0], 1))]

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.4)

    # Define Classifiers #

    names = [
        "Logistic",
        "RBF SVM",
        "Gaussian Process",
        "Naive Bayes",
        "QDA",
    ]

    classifiers = [
        LogisticRegression(),
        SVC(kernel='RBF'),
        GaussianProcessClassifier(),
        GaussianNB(),
        QuadratricDiscriminantAnalysis(),
    ]

    for name, clf in zip(names, classifiers):
        clf = make_pipeline(StandardScaler(), clf)
        clf.fit(X_train, y_train)
        score = clf.score(X_test, y_test)

        # Plot Data

        # DecisionBoundaryDisplay.from_estimator(
        #     clf, X, ax=ax # , alpha=0.8, eps=0.5
        # )

        ax.scatter(X_train[:,2], X_train[:,0], c=y_train, edgecolors="k", alpha=1.0)
        ax.scatter(X_test[:,2], X_test[:,0], c=y_train, edgecolors="k", alpha=0.8)
        #ax.sel_xlim(min, max)
        #ax.sel_ylim(min, max)

        ax.scatter(X_train[:,2], X_train[:,1], c=y_train, edgecolors="k", alpha=1.0)
        ax.scatter(X_test[:,2], X_test[:,1], c=y_train, edgecolors="k", alpha=0.8)
        #ax.sel_xlim(min, max)
        #ax.sel_ylim(min, max)

        ax.scatter(X_train[:,1], X_train[:,0], c=y_train, edgecolors="k", alpha=1.0)
        ax.scatter(X_test[:,1], X_test[:,0], c=y_train, edgecolors="k", alpha=0.8)
        #ax.sel_xlim(min, max)
        #ax.sel_ylim(min, max)


        ax.text(xmax, ymin, f"{score:.2f}", size=15), ha="right")

        plt.tight_layout()
        plt.show()


