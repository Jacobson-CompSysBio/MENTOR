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


# Test entry point
def test_entrypoint():
    exit_status = os.system('functional_partitioning --help')
    assert exit_status == 0

def test_entrypoint():
    exit_status = os.system('functional_partitioning --version')
    assert exit_status == 0


class TestInitFullranks(unittest.TestCase):
    def setUp(self):
        self.path_to_fullranks = pathlib.Path(tempfile.mkdtemp()) / 'fullranks.tsv'
        self.exit_status = os.system(f'functional_partitioning --init-test-fullranks {self.path_to_fullranks} --no-partition')

    def test_exit_status(self):
        assert self.exit_status == 0

    def test_fullranks_exists(self):
        assert self.path_to_fullranks.exists()


# END.
