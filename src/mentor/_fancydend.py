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
        heatmaps,
        reordercols,
        pcutoff,
        squish,
        relwidths,
        plotwidth
    
    ):

    command = 'Rscript ' + os.path.dirname(os.path.realpath(__file__)) + '/create_dendogram.R'
    if distances is not None:
        command += f' --distances={distances}'
    if clusters is not None:
        command += f' --clusters={clusters}'
    if map is not None:
        command += f' --map={map}'
    if outdir is not None:
        command += f' --outdir={outdir}'
        command = command + '/'
    if subcluster is True:
        command += f' --subcluster'
    if increment is not None:
        command += f' --increment={increment}'
    if maxsize is not None:
        command += f' --maxsize={maxsize}'
    if heatmaps is not None:
        command += f' --heatmaps={heatmaps}'
    if reordercols is not None:
        command += f' --reordercols'
# note: need to add pcutoff option as well
    if squish is not None:
        command += f' --squish={squish}'
    if relwidths is not None:
        command += f' --relwidths={relwidths}'
    if plotwidth is not None:
        command += f' --plotwidth={plotwidth}'
    print('\nrunning dendrogram command: ' + command)
    result = subprocess.run(command,shell = True,capture_output = True,text = True,check = False)
    print(result.stdout)
    print(result.stderr)
