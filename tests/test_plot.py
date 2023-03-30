import itertools
import joblib
import numpy as np
import pandas as pd
import os
import unittest

from numpy import testing
from sklearn import datasets, preprocessing
from matplotlib import pyplot

from functional_partitioning import _cluster as cluster
from functional_partitioning import _datasets as datasets
from functional_partitioning import _metrics as metrics
from functional_partitioning import _rwrtoolkit as rwrtoolkit
from functional_partitioning import _plot as plot

def test_make_label_mapper_labels_default():
    nodetable = pd.DataFrame(
        data=[2, 2, 1, 2, 1, 1, 1, 2, 1],
        index=['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I'],
        columns=['cluster']
    )
    label_mapper = plot.make_label_mapper(
        nodetable=nodetable,
        use_locs=[0, 1], # List.
    )
    md5sum = joblib.hash(label_mapper)
    assert  md5sum == '77cbc5e1b7a83faa49a017350c891987'

class TestLabelMapper(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._nodetable = datasets.make_nodetable()
        cls._lm_use_names_index_label = plot.make_label_mapper(nodetable=cls._nodetable, use_names=['index', 'label'])
        cls._lm_use_names_index_genename = plot.make_label_mapper(nodetable=cls._nodetable, use_names=['index', 'genename'])
        cls._lm_use_locs_0_1 = plot.make_label_mapper(nodetable=cls._nodetable, use_locs=[0,1])
        cls._lm_use_locs_0_2 = plot.make_label_mapper(nodetable=cls._nodetable, use_locs=[0,2])
        cls._lm_use_locs_2_1 = plot.make_label_mapper(nodetable=cls._nodetable, use_locs=[2,1])

    def hash_lm_use_names_index_label(self):
        md5sum = joblib.hash(self._lm_use_names_index_label)
        return md5sum

    def hash_lm_use_names_index_genename(self):
        md5sum = joblib.hash(self._lm_use_names_index_genename)
        return md5sum

    def hash_lm_use_locs_0_1(self):
        md5sum = joblib.hash(self._lm_use_locs_0_1)
        return md5sum

    def hash_lm_use_locs_0_2(self):
        md5sum = joblib.hash(self._lm_use_locs_0_2)
        return md5sum

    def hash_lm_use_locs_2_1(self):
        md5sum = joblib.hash(self._lm_use_locs_2_1)
        return md5sum

    def test_lm_use_names_index_label(self):
        md5sum = self.hash_lm_use_names_index_label()
        assert md5sum == '77cbc5e1b7a83faa49a017350c891987'

    def test_lm_use_names_index_genename(self):
        md5sum = self.hash_lm_use_names_index_genename()
        assert md5sum == '67e69b9586aab5e041d8f0b4c1d6565c'

    def test_lm_use_locs_0_1(self):
        md5sum = self.hash_lm_use_locs_0_1()
        assert md5sum == '77cbc5e1b7a83faa49a017350c891987'

    def test_lm_use_locs_0_2(self):
        md5sum = self.hash_lm_use_locs_0_2()
        assert md5sum == '67e69b9586aab5e041d8f0b4c1d6565c'

    def test_lm_use_locs_2_1(self):
        md5sum = self.hash_lm_use_locs_2_1()
        assert md5sum == '00955817cf42f979f27cd2d29ef95e0d'

# END.
