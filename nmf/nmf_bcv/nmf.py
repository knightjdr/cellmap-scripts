"""
Conducts NMF on input data.
"""

import argparse
import logging
import numpy as np
import pandas as pd

from sklearn.decomposition import NMF
from sklearn.externals import joblib


def get_arguments():
    """
    Obtains command-line arguments.

    :rtype argparse.Namespace
    """

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--input', type=argparse.FileType('rU'), required=True,
        metavar='INPUT',
        help='read input data from CSV file %(metavar)s')

    parser.add_argument(
        '--k', type=int, required=True, metavar='K',
        help='reduce the input to rank %(metavar)s')

    parser.add_argument(
        '--model-output', required=True, metavar='MODEL-OUTPUT',
        help='output the resulting model to Pickle file %(metavar)s')

    parser.add_argument(
        '--basis-output', type=argparse.FileType('w'), required=True,
        metavar='BASIS-OUTPUT',
        help='output the basis to %(metavar)s')

    parser.add_argument(
        '--score-output', type=argparse.FileType('w'), required=True,
        metavar='SCORE-OUTPUT',
        help='output the scores/coefficients to %(metavar)s')

    parser.add_argument(
        '--init', choices=('random', 'nndsvd', 'nndsvda', 'nndsvdar'),
        default='nndsvd', metavar='INIT',
        help='use method %(metavar)s to initialize values (default: '
        '%(default)s)')

    parser.add_argument(
        '--l1-ratio', type=float, default=0., metavar='L1-RATIO',
        help='use %(metavar)s as the regularization mixing parameter (use 0 '
        'to specify an L2 penalty, 1 to specify an L1 penalty, or a value in '
        '(0, 1) to specify a combination; default: %(default)s)')

    parser.add_argument(
        '--seed', type=int, default=301290193, metavar='SEED',
        help='set the seed to %(metavar)s')

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


def load_data(handle):
    """
    Loads data from the given handle.

    :param io.file handle

    :rtype pd.DataFrame
    """

    logging.info('Loading data')

    result = pd.read_csv(handle, index_col=0)

    logging.info('Loaded a table with shape {}'.format(result.shape))

    return result


def run_nmf(data, k, init, l1_ratio):
    """
    Runs NMF on the given data for a given number of factors.

    :param pd.DataFrame data

    :param int k

    :param str init

    :param float l1_ratio

    :rtype Tuple[NMF, pd.DataFrame, pd.DataFrame]
    """

    logging.info('Conducting NMF')

    nmf = NMF(n_components=k, init=init, l1_ratio=l1_ratio)

    scores = pd.DataFrame(
        nmf.fit_transform(data), index=data.index, columns=np.arange(k) + 1)

    basis = pd.DataFrame(
        nmf.components_.T, index=data.columns, columns=np.arange(k) + 1)

    basis.index.name = 'variable'

    return nmf, basis, scores


def write_output(model, basis, scores, model_output, basis_handle,
                 score_handle):
    """
    Writes the given model, basis, and scores to their respective given files.

    :param NMF model

    :param pd.DataFrame basis

    :param pd.DataFrame scores

    :param str model_output

    :param str basis_handle

    :param str score_handle
    """

    logging.info('Writing model to {}'.format(model_output))

    joblib.dump(model, model_output)

    logging.info('Writing basis to {}'.format(basis_handle.name))

    basis.to_csv(basis_handle)

    logging.info('Writing scores to {}'.format(score_handle.name))

    scores.to_csv(score_handle)


if __name__ == '__main__':

    # Get arguments.

    args = get_arguments()

    # Configure logging.

    configure_logging(args.log)

    # Conduct the analysis.

    data = load_data(args.input)

    np.random.seed(args.seed)

    nmf, basis, scores = run_nmf(data, args.k, args.init, args.l1_ratio)

    write_output(nmf, basis, scores, args.model_output,
                 args.basis_output, args.score_output)

    logging.info('Done')
