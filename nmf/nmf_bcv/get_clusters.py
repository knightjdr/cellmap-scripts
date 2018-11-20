"""
Obtains cluster assignments from the given NMF scores.
"""

import argparse
import logging
import numpy as np
import pandas as pd
import sys


def get_arguments():
    """Obtains command-line arguments."""

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--input', type=argparse.FileType('rU'), required=True,
        metavar='INPUT',
        help='read coefficients from %(metavar)s')

    parser.add_argument(
        '--output', type=argparse.FileType('w'), default=sys.stdout,
        metavar='OUTPUT',
        help='write cluster assignments to %(metavar)s')

    parser.add_argument(
        '--allow-unassigned', default=False, action='store_true',
        help='allow patients to be unassigned')

    parser.add_argument(
        '--variance-coefficient', type=float, default=1.,
        metavar='COEFFICIENT',
        help='require that scores are %(metavar)s standard deviations above 0 '
        'to be counted in calculating cluster assignment')

    parser.add_argument(
        '--log', metavar='LOG',
        help='write logging information to %(metavar)s')

    return parser.parse_args()


def configure_logging(log=None):
    """
    Configures logging.

    :param str log
    """

    if log:

        logging.basicConfig(level=logging.DEBUG, filename=log,
                            filemode='w',
                            format='%(asctime)s %(levelname)-8s %(message)s')

        console = logging.StreamHandler()

        console.setLevel(logging.INFO)

        console.setFormatter(logging.Formatter('%(message)s'))

        logging.getLogger().addHandler(console)

    else:

        logging.basicConfig(level=logging.INFO, format='%(message)s')


def load_scores(handle):
    """
    Loads scores from the given handle.

    :param io.file handle

    :rtype pd.DataFrame
    """

    logging.info('Loading scores')

    result = pd.read_csv(handle, index_col=0)

    logging.info('Loaded a table with shape {}'.format(result.shape))

    return result


def get_mle_variance(series, mean=None):
    """
    Obtains the MLE variance given a series of values and a calculated mean.

    :param pd.Series[float] series

    :param float mean

    :rtype float
    """

    if mean is None:

        mean = series.mean()

    return 1 / series.size * ((series - mean) ** 2).sum()


def get_clusters(df, allow_unassigned, variance_coefficient):
    """
    Obtains cluster assignments from the given scores.

    :param pd.DataFrame df

    :param bool allow_unassigned

    :param float variance_coefficient

    :rtype pd.Series[int]
    """

    logging.info('Calculating cluster assignments')

    # Calculate minimum thresholds to call cluster assignments.

    min_thresholds = pd.Series(np.tile(1e-6, df.shape[1]), index=df.columns)

    if allow_unassigned:

        # Estimate the variance for each factor by fixing the MLE estimate of
        # the mean to 0.

        min_thresholds = df.apply(get_mle_variance, mean=0).apply(
            np.sqrt) * variance_coefficient

    unassigned_mask = (df >= min_thresholds).sum(axis=1) < 1

    result = df.apply(pd.Series.argmax, axis=1)

    result.loc[unassigned_mask] = 0

    return result


def write_output(series, handle):
    """
    Writes the given cluster assignments to the given handle.

    :param pd.Series[int] series

    :param str handle
    """

    logging.info('Writing output')

    df = series.reset_index()

    df.columns = ['id', 'classification']

    df.to_csv(handle, index=False)


if __name__ == '__main__':

    # Get arguments.

    args = get_arguments()

    # Configure logging.

    configure_logging(args.log)

    # Load the scores.

    scores = load_scores(args.input)

    # Obtain cluster assignments.

    clusters = get_clusters(
        scores, allow_unassigned=args.allow_unassigned,
        variance_coefficient=args.variance_coefficient)

    # Write the output.

    write_output(clusters, args.output)

    logging.info('Done')
