import pandas as pd
import numpy as np
import logging

from functional_partitioning import _metrics as metrics

LOGGER= logging.getLogger(__name__)

def get_top_ranked(scores, ranks=None, max_rank='elbow'):
    if not isinstance(scores, pd.DataFrame):
        scores = pd.DataFrame(scores)

    if ranks is None:
        ranks = scores.rank(axis=1, ascending=False)
    elif not isinstance(ranks, pd.DataFrame):
        ranks = pd.DataFrame(ranks)

    if max_rank == 'elbow':
        # Find elbow and set max_rank.
        mean_scores = scores.mean()
        max_rank = metrics.get_elbow(mean_scores)
        LOGGER.info(f'Set max_rank to {max_rank}.')

    # Filter the rank vectors.
    mask = (ranks <= max_rank).fillna(False)
    col_mask = mask.any()

    features = scores.loc[:, col_mask]

    # Re-rank them.
    features = features.rank(axis=1, method='first', ascending=False)
    
    for i in features.index:
        # Fill the seed in each vector as 0; ie, give the seed the best rank.
        if np.isnan(features.loc[i, i]):
            features.loc[i, i] = 0
    
    return features


# END.
