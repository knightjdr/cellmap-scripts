"""
Conducts cross-validated NMF to determine the number of factors to use,
using the bi-cross-validation procedure outlined in Owen and Perry.
"""

import argparse
import feather
import functools
import logging
import multiprocessing
import numpy as np
import pandas as pd
import tqdm

from sklearn.decomposition import NMF


def get_arguments():
    """Obtains command-line arguments."""

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--input', type=argparse.FileType('rU'), required=True,
        metavar='INPUT',
        help='read input data from CSV file %(metavar)s')

    parser.add_argument(
        '--output', required=True, metavar='OUTPUT',
        help='write Q2 values to Feather file %(metavar)s')

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
        '--k', type=int, nargs='+', metavar='K',
        help='calculate Q2 for given ranks %(metavar)s')

    parser.add_argument(
        '--folds', type=int, default=3, metavar='FOLDS',
        help='run bi-cross-validation with %(metavar)s folds (default: '
        '%(default)s)')

    parser.add_argument(
        '--iterations', type=int, default=2000, metavar='ITERATIONS',
        help='run bi-cross-validation with %(metavar)s iterations (default: '
        '%(default)s)')

    parser.add_argument(
        '--seed', type=int, default=290348203, metavar='SEED',
        help='set the seed to %(metavar)s')

    parser.add_argument(
        '--cores', type=int, metavar='CORES',
        help='use %(metavar)s cores for the analysis')

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


def get_k(data, k, folds):
    """
    Determines the values of k to use.

    :param pd.DataFrame data

    :param List[int] k

    :param int folds

    :rtype List[int]
    """

    result = k if k is not None else np.arange(
        np.floor(min(data.shape) * (1 - 1 / folds)) - 1, dtype=int) + 1

    logging.info('Setting k = {}'.format(', '.join(str(i) for i in result)))

    return list(result)


def get_seeds(seed, iterations):
    """
    Obtains seeds for the given number of iterations.

    :param int seed

    :param int iterations

    :rtype np.array<int>
    """

    logging.info('Generating seeds')

    np.random.seed(seed)

    int_info = np.iinfo(np.uint32)

    return np.random.randint(0, int_info.max, iterations)


def run_k(train_data, data_left, data_right, test_data, k, init, l1_ratio):
    """Runs cross-validation with the given training data and given test data
    for a given number of factors, returning Q2.

    :param pd.DataFrame train_data

    :param pd.DataFrame data_left

    :param pd.DataFrame data_right

    :param int k

    :parma str init

    :param float l1_ratio

    :rtype float
    """

    nmf = NMF(n_components=int(k), tol=1e-6,
              max_iter=200, init=init, l1_ratio=l1_ratio)

    coefficients = nmf.fit_transform(train_data)

    test_reconstructions = data_left @ np.linalg.pinv(nmf.components_) @ np.linalg.pinv(coefficients) @ data_right

    return 1 - ((test_data - test_reconstructions) ** 2).sum().sum() / (test_data ** 2).sum().sum()


def run_fold(data, k, fold_number, sample_folds, feature_folds, init,
             l1_ratio):
    """
    Runs cross-validation with the given training data and the given test
    data, returning the squared reconstruction error for that fold.

    :param pd.DataFrame data

    :param List[int] k

    :param int fold_number

    :param np.array[int] sample_folds

    :param np.array[int] feature_folds

    :param str init

    :param float l1_ratio

    :rtype pd.DataFrame[int]
    """

    fold_number += 1

    train_data = data.loc[
        sample_folds != fold_number, feature_folds != fold_number]

    data_left = data.loc[
        sample_folds == fold_number, feature_folds != fold_number]

    data_right = data.loc[
        sample_folds != fold_number, feature_folds == fold_number]

    test_data = data.loc[
        sample_folds == fold_number, feature_folds == fold_number]

    result = [run_k(train_data, data_left, data_right, test_data, k_, init,
                    l1_ratio) for k_ in k]

    return pd.DataFrame.from_items([
        ('fold', fold_number),
        ('k', k),
        ('q2', result)
    ])


def get_folds(num_items, folds):
    """Obtains fold numbers for the given number of items and folds.

    :param int num_items

    :param int folds

    :rtype np.array<int>
    """

    n_per_fold = num_items // folds

    fold_numbers = np.zeros(num_items, dtype=int)

    remaining_indices = np.arange(num_items)

    for fold in range(folds):

        if fold + 1 >= folds:

            fold_numbers[remaining_indices] = fold + 1

        else:

            fold_indices = np.random.choice(
                remaining_indices, n_per_fold, replace=False)

            remaining_indices = np.setdiff1d(remaining_indices, fold_indices)

            fold_numbers[fold_indices] = fold + 1

    return fold_numbers


def run_iteration(data, folds, k, init, l1_ratio, iteration_seed):
    """
    Runs cross-validation for an iteration on the given data, given a number
    of folds and a random seed.

    :param pd.DataFrame data

    :param int folds

    :param Int[list] k

    :param Tuple[int, int] iteration_seed

    :param str init

    :param float l1_ratio

    :rtype pd.DataFrame
    """

    iteration, seed = iteration_seed

    np.random.seed(seed)

    sample_folds = get_folds(data.shape[0], folds)

    feature_folds = get_folds(data.shape[1], folds)

    q2 = []

    for fold in range(folds):

        train_data = data.loc[sample_folds != (fold + 1)]

        test_data = data.loc[sample_folds == (fold + 1)]

        q2.append(
            run_fold(data, k, fold, sample_folds, feature_folds, init,
                     l1_ratio))

    q2 = pd.concat(q2)

    q2.insert(0, 'iteration', iteration + 1)

    return q2


def write_output(q2, filename):
    """
    Writes the given data to the given Feather output.

    :param pd.DataFrame q2

    :param str filename
    """

    logging.info('Writing output to {}'.format(filename))

    feather.write_dataframe(q2, filename)


if __name__ == '__main__':

    # Get arguments.

    args = get_arguments()

    # Configure logging.

    configure_logging(args.log)

    # Load the data.

    data = load_data(args.input)

    # Determine the number of ranks to use.

    k = get_k(data, args.k, args.folds)

    # Generates a bunch of seeds.

    seeds = get_seeds(args.seed, args.iterations)

    # Conducts cross-validated NMF.

    logging.info('Conducting cross-validation')

    f_partial = functools.partial(
        run_iteration, data, args.folds, k, args.init, args.l1_ratio)

    pool = multiprocessing.Pool(args.cores)

    # q2 = pd.concat(
    # list(tqdm.tqdm(map(f_partial, enumerate(seeds)),
    # total=args.iterations)))

    q2 = pd.concat(list(tqdm.tqdm(
        pool.imap_unordered(f_partial, enumerate(seeds)),
        total=args.iterations)))

    pool.close()

    pool.join()

    # Write the output.

    write_output(q2, args.output)

    logging.info('Done')
