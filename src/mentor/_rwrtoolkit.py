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

if 'PATH_TO_RWRTOOLKIT' in os.environ:
    PATH_TO_RWRTOOLKIT = os.getenv('PATH_TO_RWRTOOLKIT')
else:
    PATH_TO_RWRTOOLKIT = None

if 'PATH_TO_CONDA_ENV' in os.environ:
    PATH_TO_CONDA_ENV = os.getenv('PATH_TO_CONDA_ENV')
else:
    PATH_TO_CONDA_ENV = None

def _activate_env(
    path_to_conda_env=PATH_TO_CONDA_ENV
):
    '''
    Return the command to activate conda environment.
    Parameters
    ----------
    path_to_conda_env : str, None
        If provided, activate the conda environment.  This should be the *path*
        to the env, NOT the name of the env. Default None.
    path_to_rwrtoolkit : str, None
        If provided, set the environment variable PATH_TO_RWRTOOLS. Default None.

    Returns
    -------
    str
    '''
    path_to_conda_env = os.path.expanduser(path_to_conda_env)
    path_to_conda_env = os.path.realpath(path_to_conda_env)
    command = (
        'set +eu && '
        'PS1=dummy && '
        'source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh && '
        f'conda activate "{path_to_conda_env}" && '
        'echo "CONDA_PREFIX is $CONDA_PREFIX" '
    )
    return command

def rwr_singletons(
    # path_to_conda_env=None,
    # path_to_rwrtoolkit=PATH_TO_RWRTOOLKIT,
    data=None,
    geneset=None,
    method='singletons',
    folds=None,
    restart=None,
    tau=None,
    numranked=None,
    outdir=None,
    modname=None,
    #plot=None,
    out_fullranks=None,
    out_medianranks=None,
    threads=None,
    verbose=None
):
    '''
    Return the command to run RWR-singletons.
    Usage: /Users/m8z/src/RWRtoolkit/inst/scripts/run_cv.R [options]
    Options:
        -d DATA, --data=DATA
            The path to the .Rdata file for your combo of underlying functional
            networks. This file is produced by RWR_make_multiplex.
        -g GENESET, --geneset=GENESET
            The path to the gene set file. It must have the following first two
            columns with no headers tab-delimited: <setid> <gene> <weight>.
        --method=METHOD
            Cross-validation method. `kfold`, `loo`, or `singletons`. [default
            kfold]
        -f FOLDS, --folds=FOLDS
            Number (k) of folds to use in k-fold CV. [default 5]
        -r RESTART, --restart=RESTART
            Set the restart parameter [0,1). Higher value means the walker will
            jump back to a seed node more often. [default 0.7]
        --tau=TAU
            comma-separated list of values between that MUST add up to the
            number of network layers in the .Rdata file. One value per network
            layer that determines the probability that the random walker will
            restart in that layer. e.g. if there are three layers (A,B,C) in
            your multiplex network, then --tau '0.2,1.3,1.5' will mean that
            layer A is less likely to be walked on after a restart than layers
            B or C. [default 1.0]
        -n NUMRANKED, --numranked=NUMRANKED
            proportion of ranked genes to return [0,1]. e.g. 0.1 will return
            the top 10%. [default 1]
        -o OUTDIR, --outdir=OUTDIR
            Path to the output directory. Both 'fullranks' and 'medianranks'
            will be saved here with auto-generated filenames. (--out-fullranks
            and --out-medianranks override this.)
        -m MODNAME, --modname=MODNAME
            String to include in output filename. (--out-fullranks and
            --out-medianranks override this.)
        -p, --plot
            Output plots of ROC, PRC, NDCG etc. [default FALSE]
        --out-fullranks=OUT-FULLRANKS
            Specify the full path for full results. Ignore --outdir and
            --modname and use this instead.
        --out-medianranks=OUT-MEDIANRANKS
            Specify the full path for median results. Ignore --outdir and
            --modname and use this instead.
        -t THREADS, --threads=THREADS
            Number of threads to use. Default for your system is all cores - 1.
            [default 11]
        -v, --verbose
            Verbose mode. [default FALSE]
        -h, --help
            Show this help message and exit
    Paramters
    ---------
    tau : list
        List of values will be converted to string for RWR_KFOLD.R
    '''
    # if path_to_conda_env is not None:
    #     command = _activate_env(path_to_conda_env) + ' && Rscript'
    # else:
    command = 'Rscript'
    exe = os.path.dirname(os.path.realpath(__file__)) + '/run_cv.R'
    command += f' {exe}'
    # if path_to_rwrtoolkit is not None:
    #     exe = os.path.join(path_to_rwrtoolkit, 'inst', 'scripts', 'run_cv.R')
    #     command += f' "{exe}"'
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
    #if plot is not None:
    #    command += f' --plot'
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

def run(commands, sep=' && ', dry_run=False, verbose=0):
    '''
    Run a shell command in a subprocess.
    Parameters
    ----------
    commands : str, list
        List of command elements to pass to subprocess.
    dry_run : bool
        If True, do nothing.
    Returns
    -------
    dict
        A copy of `pset` with the following additional keys from the subprocess result:
            {command_key}_stdout
            {command_key}_stderr
            {command_key}_returncode
    '''
    if isinstance(commands, str):
        command = commands
    else:
        command = sep.join(commands)
    result = dict(
        command=command,
        dry_run=dry_run,
        returncode=None,
        stderr=None,
        stdout=None
    )
    if not dry_run:
        print(command)
        res = subprocess.run(command, shell=True, capture_output=True)
        if res.returncode == 0:
            if verbose > 0:
                print('Success!')
        else:
            print('here')
            #print('[WARNING] Command exited with non-zero status: {}'.format(res.returncode), file=sys.stderr)
            #print(command, file=sys.stderr)
            #print(res.stderr.decode(), file=sys.stderr)
        result.update(
            stdout=res.stdout.decode(),
            stderr=res.stderr.decode(),
            returncode=res.returncode,
            args=res.args,
            dry_run=dry_run
        )
    else:
        print('[DRY RUN]', command)
    return result

def gzip(source, target=None, compresslevel=9, encoding=None, errors=None, newline=None, in_place=True):
    if isinstance(source, str):
        source = pathlib.Path(source)
    if target is None:
        target = source.with_suffix(source.suffix + '.gz')
    with open(source, 'rb') as f_in:
        with gzip_module.open(target, 'wb', compresslevel=compresslevel, encoding=encoding, errors=errors, newline=newline) as f_out:
            shutil.copyfileobj(f_in, f_out)
    if in_place and target.exists():
        os.remove(source)

def compress_results(dirname, **kwargs):
    '''
    Gzip the RWR results.
    '''
    if isinstance(dirname, str):
        dirname = pathlib.Path(dirname)
    for f in dirname.glob('RWR*.tsv'):
        gzip(f, **kwargs)

def fullranks_to_matrix(path_or_dataframe, to='scores', drop_missing=True):
    '''
    Convert RWR "fullranks" file to matrix.
    Parameters
    ----------
    path_or_dataframe : str, pd.DataFrame
        Path to 'fullranks' file from `RWR-CV --method=singletons` or
        pandas.DataFrame.
    max_rank : int, str
        Maximum rank to use for clustering. If 'elbow', use elbow method to
        determine max_rank.
    drop_missing : bool
        Drop genes that are labeled "missing" in the fullranks file.
    Returns
    -------
    X : pd.DataFrame
    '''
    if isinstance(path_or_dataframe, pd.DataFrame):
        fullranks = path_or_dataframe
    else:
        fullranks = pd.read_table(path_or_dataframe)
    if drop_missing:
        fullranks = fullranks.query('seed!="missing"')
    if to.lower().startswith('score'):
        mat = fullranks.pivot(index='seed', columns='NodeNames', values='Score')
    elif to.lower().startswith('rank'):
        mat = fullranks.pivot(index='seed', columns='NodeNames', values='rank')
    else:
        raise ValueError(f'Invalid value for `to`: {to}')
    return mat

def transform_fullranks(path_or_dataframe, drop_missing=True, max_rank='elbow'):
    if isinstance(path_or_dataframe, pd.DataFrame):
        fullranks = path_or_dataframe
    else:
        fullranks = pd.read_table(path_or_dataframe)
    if drop_missing:
        fullranks = fullranks.query('seed!="missing"')
    if max_rank == 'elbow':
        y = fullranks.groupby('rank')['Score'].mean()
        max_rank = metrics.get_elbow(y)
        LOGGER.info(f'Set max_rank to {max_rank}.')
    ranks = fullranks_to_matrix(fullranks, to='rank', drop_missing=drop_missing)
    scores = fullranks_to_matrix(fullranks, to='scores', drop_missing=drop_missing)
    labels = scores.index.to_list()
    assert (ranks.index == scores.index).all()
    assert (ranks.columns == scores.columns).all()
    mask = (ranks <= max_rank).fillna(False)
    col_mask = mask.any()
    features = scores.loc[:, col_mask]
    features = features.rank(axis=1, method='first', ascending=False)
    for i in labels:
        if ( i in features.columns.to_list() ) and np.isnan(features.loc[i, i]):
            features.loc[i, i] = 0
    return features, labels
