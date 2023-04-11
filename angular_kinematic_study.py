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

import matplotlib as mpl
import matplotlib.pyplot as plt

# -----------------------#
# Sci-kit Learn Imports #

import sklearn as skl

from sklearn.preprocessing import StandardScaler, MinMaxScaler, RobustScaler
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline, Pipeline

from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis, QuadraticDiscriminantAnalysis

from sklearn.metrics import roc_curve, roc_auc_score

# -----------------------#

from primary_analysis import bins  # Eventually move bins to a file in src?

if __name__ == "__main__":
    print("STARTING")
    cur_dir = os.getcwd()
    in_data_dir = cur_dir + "/data/processed/"
    out_data_dir = cur_dir + "/data/processed/"
    out_dir = cur_dir + "/reports/weekly/apr10/"

    stat = "ca_0p6"
    ending = f"_{stat}_.png"

    pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)  # make out directory if it doesn't exist
    pathlib.Path(out_data_dir).mkdir(parents=True, exist_ok=True)  # make out data directory if it doesn't exist
    print(f"Using output directory '{out_dir}'")

    ########
    fname_mc = in_data_dir + "metdRdEtadPhi_mc_" + stat + ".npy"
    fname_r3 = in_data_dir + "metdRdEtadPhi_r3_" + stat + ".npy"

    print("LOADING DATA")
    vars_mc, vars_r3 = np.load(fname_mc), np.load(fname_r3)
    print(f"MC LEN: {len(vars_mc):,} | R3 LEN: {len(vars_r3):,}")
    
    #ns = int((len(vars_mc)+len(vars_r3))/2)
    #print("Resampling to balance classes")
    #vars_mc = skl.utils.resample(vars_mc, n_samples=ns)
   # vars_r3 = skl.utils.resample(vars_r3, n_samples=ns)
    
    # print(f"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    # print(f"! TAKING ABSOLUTE VALUE OF DATA !")
    # print(f"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    # vars_mc, vars_r3 = np.abs(vars_mc), np.abs(vars_r3)

    # Load Data #
    X = np.r_[vars_mc, vars_r3]
    y = np.r_[np.ones((vars_mc.shape[0],)), np.zeros((vars_r3.shape[0],))]
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.4)

    data_bin_names = ["met", "dR", "dEta", "dPhi"]
    data_plot_labels = [r"$E_{T}^{miss}$ [GeV]",r"$\Delta R$", r"$\Delta\eta$", r"$\Delta\Phi$"]
    # Define Classifiers #

    names = [
        "Logistic",
        "RBF SVM",
        "Gaussian Process",
        "Naive Bayes",
        "LDA",
        "QDA",
    ]

    classifiers = [
        LogisticRegression(),
        SVC(kernel="rbf", probability=True),
        GaussianProcessClassifier(),
        GaussianNB(),
        LinearDiscriminantAnalysis(),
        QuadraticDiscriminantAnalysis(),
    ]

    fig, axs = plt.subplots(6, len(names), figsize=((len(names)) * 3, 6 * 3))
    cmap = "RdBu"

    for iclf, (name, clf) in enumerate(zip(names, classifiers)):
        print(f"({iclf+1:,}/{len(classifiers):,}) - {name=}")
        clf = make_pipeline(StandardScaler(), clf)
        clf.fit(X_train, y_train)
        try:
            y_pred = clf.decision_function(X)
            y_test_pred = clf.decision_function(X_test)
        except:
            print('here')
            y_pred = clf.predict_proba(X)
            y_test_pred = clf.predict_proba(X_test)

        if len(y_pred.shape) > 1:
            y_pred = y_pred[:, 1].flatten()
            y_test_pred = y_test_pred[:, 1].flatten()
        score = clf.score(X_test, y_test)
        auc = roc_auc_score(y_test, y_test_pred)
        print(f"\tScore: {score:>0.3f}, AUC: {auc:>0.3f}")

        # Plot Data
        for iplot, (ix, iy) in enumerate(zip((3, 3, 2, 1, 2, 3), (1, 2, 1, 0, 0, 0))):
            ax = axs[iplot, iclf]

            # DecisionBoundaryDisplay.from_estimator(
            #     clf, X, ax=ax # , alpha=0.8, eps=0.5
            # )

            xl, yl = data_plot_labels[ix], data_plot_labels[iy]
            _, xmin, xmax = bins[data_bin_names[ix]]
            _, ymin, ymax = bins[data_bin_names[iy]]

            # ax.scatter(X_train[:,ix], X_train[:,iy], c=y_train, edgecolors="k", alpha=1.0)
            # ax.scatter(X_test[:,ix], X_test[:,iy], c=y_train, edgecolors="k", alpha=0.8)

            xmc, xr3 = X[y == 1], X[y == 0]
            ypmc, ypr3 = y_pred[y == 1], y_pred[y == 0]
            vmin, vmax = np.min(y_pred), np.max(y_pred)
            print(vmin, vmax)

            if iplot == 0:
                ax.set_title(name)
            ax.set_xlabel(xl)
            ax.set_ylabel(yl)

            ax.scatter(xr3[:, ix], xr3[:, iy], c=ypr3, cmap=cmap, vmin=vmin, vmax=vmax, marker="X", alpha=0.1)  # , edgecolors="k")
            ax.scatter(xmc[:, ix], xmc[:, iy], c=ypmc, cmap=cmap, vmin=vmin, vmax=vmax, marker="o", alpha=0.1, edgecolors="k")

            ax.set_xlim(xmin, xmax)
            ax.set_ylim(ymin, ymax)

            ax.text(xmax, ymin, f"{auc:.3f}", size=15, ha="right")

    plt.tight_layout()
    plt.savefig(out_dir + "ang_kin_compare_cuts" + ending)
    plt.show()
