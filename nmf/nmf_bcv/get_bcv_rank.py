"""
Obtains the rank with the maximum mean Q2 from bi-cross-validation.

Some runs may have extremely poor values of Q2. By default, only Q2 values
greater than -1 will be considered. Specify --q2-threshold to change this
threshold.
"""

import argparse
import feather
import logging


def get_arguments():
    """
    Obtains command-line arguments.

    :rtype argparse.Namespace
    """

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--input', type=argparse.FileType('r'), required=True, metavar='INPUT',
        help='read input Q2 values from Feather file %(metavar)s')

    parser.add_argument(
        '--output', type=argparse.FileType('w'), required=True,
        metavar='OUTPUT',
        help='output the optimal rank to text file %(metavar)s')

    parser.add_argument(
        '--q2-threshold', type=float, default=-1., metavar='THRESHOLD',
        help='exclude values of Q2 less than %(metavar)s from consideration')

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
    Loads values of Q2 from the given handle.

    :param io.file handle

    :rtype pandas.DataFrame
    """

    logging.info('Loading Q2 values from {}'.format(handle.name))

    result = feather.read_dataframe(handle.name)

    logging.info('Loaded a table with shape {}'.format(result.shape))

    return result


def filter_runs(df, threshold):
    """
    Filters runs in the given data frame whose Q2 is less than the given
    threshold.

    :param pandas.DataFrame df

    :param float threshold

    :rtype pandas.DataFrame
    """

    logging.info('Filtering runs with Q2 < {}'.format(threshold))

    result = df.loc[df['q2'] > threshold]

    logging.info('Result is a table with shape {}'.format(result.shape))

    return result


def get_optimal_rank(df):
    """
    Obtains the rank with the maximum mean Q2 from the given data frame.

    :param pandas.DataFrame df

    :rtype int
    """

    logging.info('Calculating optimal rank')

    means = df.groupby('k')['q2'].mean()

    return means.index[means == means.max()][0]


def write_output(k, handle):
    """
    Writes the optimal rank to the given handle.

    :param int k

    :param str handle
    """

    logging.info('Writing output')

    handle.write(str(k))


if __name__ == '__main__':

    # Get arguments.

    args = get_arguments()

    # Configure logging.

    configure_logging(args.log)

    # Conduct the analysis.

    q2 = load_data(args.input)

    q2 = filter_runs(q2, args.q2_threshold)

    k = get_optimal_rank(q2)

    write_output(k, args.output)

    logging.info('Done')
