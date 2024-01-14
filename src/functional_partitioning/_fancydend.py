'''
Wrapper functions for fancy dendrograms
'''

import gzip as gzip_module
import logging
import numpy as np
import os
import pandas as pd
import pathlib
import shlex
import shutil
import subprocess

def fancy_dendrogram(
        
        distances,
        clusters,
        map,
        outdir,
        subcluster,
        increment,
        maxsize,
        export,
        heatmaps,
        pcutoff,
        squish,
        relwidths,
        plotwidth
    
    ):

    command = 'Rscript ' + os.path.dirname(os.path.realpath(__file__)) + '/fp_create_dendogram.R'
    if distances is not None:
        command += f' --distances "{distances}"'
    if clusters is not None:
        command += f' --clusters "{clusters}"'
    if map is not None:
        command += f' --map "{map}"'
    if outdir is not None:
        command += f' --outdir "{outdir}"'
    if subcluster is not None:
        command += f' --subcluster "{subcluster}"'
    if increment is not None:
        command += f' --increment "{increment}"'
    if maxsize is not None:
        command += f' --maxsize "{maxsize}"'
    if export is not None:
        command += f' --export "{export}"'
    if heatmaps is not None:
        command += f' --heatmaps "{heatmaps}"'
    if squish is not None:
        command += f' --squish "{squish}"'
    if relwidths is not None:
        command += f' --relwidths "{relwidths}"'
    if plotwidth is not None:
        command += f' --plotwidth "{plotwidth}"'
    subprocess.run(command,shell = True,capture_output = False)

