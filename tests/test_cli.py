import itertools
import joblib
import numpy as np
import pandas as pd
import os
import unittest
import tempfile
import pathlib

from numpy import testing
from sklearn import datasets, preprocessing
from matplotlib import pyplot

from functional_partitioning import _cluster as cluster
from functional_partitioning import _datasets as datasets
from functional_partitioning import _metrics as metrics
from functional_partitioning import _rwrtoolkit as rwrtoolkit
from functional_partitioning import _plot as plot

from scipy.spatial import distance
from scipy.cluster import hierarchy

# Fixtures.
def _initialize_test_fullranks():
    path_to_fullranks = pathlib.Path(tempfile.mkdtemp()) / 'fullranks.tsv'
    os.system(f'functional_partitioning --init-test-fullranks {path_to_fullranks} --no-partition')
    return path_to_fullranks

def _initialize_outdir():
    outdir = pathlib.Path(tempfile.mkdtemp())
    return outdir


# Test entry point
def test_entrypoint():
    exit_status = os.system('functional_partitioning --help')
    assert exit_status == 0

def test_version():
    exit_status = os.system('functional_partitioning --version')
    assert exit_status == 0

# Test that the fullranks file is created.
def test_init_fullranks_exit():
    path_to_fullranks = pathlib.Path(tempfile.mkdtemp()) / 'fullranks.tsv'
    exit_status = os.system(f'functional_partitioning --init-test-fullranks {path_to_fullranks} --no-partition')
    assert exit_status == 0

def test_init_fullranks_exists():
    path_to_fullranks = pathlib.Path(tempfile.mkdtemp()) / 'fullranks.tsv'
    os.system(f'functional_partitioning --init-test-fullranks {path_to_fullranks} --no-partition')
    assert path_to_fullranks.exists()

# class TestInitFullranks(unittest.TestCase):
#     '''Basic example of how to use setUpClass and tearDownClass'''
#     @classmethod
#     def setUpClass(cls):
#         cls._path_to_fullranks = _initialize_test_fullranks()
# 
#     @classmethod
#     def tearDownClass(cls):
#         os.remove(cls._path_to_fullranks)
# 
#     def test_foo(self):
#         assert os.path.exists(self._path_to_fullranks)

# Test the CLI.
class TestFunctionalPartitioning(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._path_to_fullranks = _initialize_test_fullranks()
        cls._outdir = _initialize_outdir()
        cls._exit_code = os.system(f'functional_partitioning --rwr-fullranks {cls._path_to_fullranks} --outdir {cls._outdir}')

    @classmethod
    def tearDownClass(cls):
        os.remove(cls._path_to_fullranks)

    def test_exit_code(self):
        assert self._exit_code == 0

    def test_fullranks_exists(self):
        assert self._path_to_fullranks.exists()

    def test_outdir_exists(self):
        assert self._outdir.exists()

    def test_clusters_exists(self):
        assert (self._outdir / 'clusters.tsv').exists()

    def test_dendrogram_exists(self):
        assert (self._outdir / 'dendrogram.png').exists()

    def test_dissimilarity_matrix_exists(self):
        assert (self._outdir / 'dissimilarity-matrix.tsv').exists()

    def test_dissimilarity_stats_exists(self):
        assert (self._outdir / 'dissimilarity-stats.tsv').exists()

    def test_plot_exists(self):
        assert (self._outdir / 'distribution-of-pairwise-dissimilarities.png').exists()


class TestOutClusters(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._path_to_fullranks = _initialize_test_fullranks()
        cls._outdir = _initialize_outdir()
        cls._out_clusters = cls._outdir / 'clusters-foo.tsv'
        cls._exit_code = os.system(f'functional_partitioning --rwr-fullranks {cls._path_to_fullranks} --outdir {cls._outdir} --out-clusters {cls._out_clusters}')

    @classmethod
    def tearDownClass(cls):
        os.remove(cls._path_to_fullranks)

    def test_exit_code(self):
        assert self._exit_code == 0

    def test_clusters_exists(self):
        assert self._out_clusters.exists()


# END.
