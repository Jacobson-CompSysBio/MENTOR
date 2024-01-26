'''
Wrapper functions for RWRtoolkit
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
import sys
from mentor import _metrics as metrics
LOGGER = logging.getLogger(__name__)


def rwr_singletons(
    data = None,
    geneset = None,
    method = 'singletons',
    folds = None,
    restart = None,
    tau = None,
    numranked = None,
    outdir = None,
    modname = None,
    out_fullranks = None,
    out_medianranks = None,
    threads = None,
    verbose = None
):
    command = 'Rscript'
    exe = os.path.dirname(os.path.realpath(__file__)) + '/run_cv.R'
    command += f' {exe}'
    if data is not None:
        command += f' --data "{data}"'
    if geneset is not None:
        command += f' --geneset "{geneset}"'
    if method is not None:
        command += f' --method "{method}"'
    if folds is not None:
        command += f' --folds "{folds}"'
    if restart is not None:
        command += f' --restart "{restart}"'
    if tau is not None:
        if isinstance(tau, (int, float)):
            tau = str(tau)
        elif isinstance(tau, str):
            pass
        else:
            tau = ','.join(map(str, tau))
        command += f' --tau "{tau}"'
    if numranked is not None:
        command += f' --numranked "{numranked}"'
    if modname is not None:
        command += f' --modname "{modname}"'
    if out_fullranks is not None:
        command += f' --out-fullranks "{out_fullranks}"'
    if out_medianranks is not None:
        command += f' --out-medianranks "{out_medianranks}"'
    if outdir is not None:
        command += f' --outdir "{outdir}"'
    if threads is not None:
        command += f' --threads "{threads}"'
    if verbose is not None:
        command += f' --verbose'
    return command

def run(commands,sep = ' && ',dry_run = False,verbose = 0):
    if isinstance(commands, str):
        command = commands
    else:
        command = sep.join(commands)
    result = dict(
        command = command,
        dry_run = dry_run,
        returncode = None,
        stderr = None,
        stdout = None
    )
    if not dry_run:
        print(command)
        res = subprocess.run(command,shell = True,capture_output = True)
        print(res)
        if res.returncode == 0:
            if verbose > 0:
                print('Success!')
        else:
            print('Failed')
        result.update(
            stdout = res.stdout.decode(),
            stderr = res.stderr.decode(),
            returncode = res.returncode,
            args = res.args,
            dry_run = dry_run
        )
    else:
        print('[DRY RUN]', command)
    return result

def gzip(source,target = None,compresslevel = 9,encoding = None,errors = None,newline = None,in_place = True):
    if isinstance(source,str):
        source = pathlib.Path(source)
    if target is None:
        target = source.with_suffix(source.suffix + '.gz')
    with open(source, 'rb') as f_in:
        with gzip_module.open(target,'wb',compresslevel = compresslevel,encoding = encoding,errors = errors,newline = newline) as f_out:
            shutil.copyfileobj(f_in,f_out)
    if in_place and target.exists():
        os.remove(source)

def compress_results(dirname,**kwargs):
    if isinstance(dirname,str):
        dirname = pathlib.Path(dirname)
    for f in dirname.glob('RWR*.tsv'):
        gzip(f,**kwargs)

def fullranks_to_matrix(path_or_dataframe,to = 'scores',drop_missing = True):
    if isinstance(path_or_dataframe,pd.DataFrame):
        fullranks = path_or_dataframe
    else:
        fullranks = pd.read_table(path_or_dataframe)
    if drop_missing:
        fullranks = fullranks.query('seed!="missing"')
    if to.lower().startswith('score'):
        mat = fullranks.pivot(index = 'seed',columns = 'NodeNames',values = 'Score')
    elif to.lower().startswith('rank'):
        mat = fullranks.pivot(index = 'seed',columns = 'NodeNames',values = 'rank')
    else:
        raise ValueError(f'Invalid value for `to`: {to}')
    return mat

def transform_fullranks(path_or_dataframe,drop_missing = True,max_rank = 'elbow'):
    if isinstance(path_or_dataframe,pd.DataFrame):
        fullranks = path_or_dataframe
    else:
        fullranks = pd.read_table(path_or_dataframe)
    if drop_missing:
        fullranks = fullranks.query('seed!="missing"')
    if max_rank == 'elbow':
        y = fullranks.groupby('rank')['Score'].mean()
        max_rank = metrics.get_elbow(y)
        LOGGER.info(f'Set max_rank to {max_rank}.')
    ranks = fullranks_to_matrix(fullranks,to = 'rank',drop_missing = drop_missing)
    scores = fullranks_to_matrix(fullranks,to = 'scores',drop_missing = drop_missing)
    labels = scores.index.to_list()
    assert (ranks.index == scores.index).all()
    assert (ranks.columns == scores.columns).all()
    mask = (ranks <= max_rank).fillna(False)
    col_mask = mask.any()
    features = scores.loc[:, col_mask]
    features = features.rank(axis = 1,method = 'first',ascending = False)
    for i in labels:
        if (i in features.columns.to_list()) and np.isnan(features.loc[i, i]):
            features.loc[i, i] = 0
    return features, labels
