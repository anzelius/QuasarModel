from datetime import datetime
import time
import gzip
import os
import logging
import logging.handlers

from gauss import Gaussian, GaussList
from numpy import radians
from astropy.io import fits
from constants import one_source1, one_source2, one_source3, one_source4, one_source5, one_source6, \
    two_sources

import csv
from logging import getLogger
from source_model import SourceModel
from numpy import sqrt


"""
##################################################
Tools for logging 
##################################################
"""


# Custom filter use to format records
class ContextFilter(logging.Filter):
    def filter(self, record):
        setattr(record, 'utc', datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S.%f')[:-3])
        return True


# Set default logger
def set_logger(log_path='', prefix='', console=False, size=1000000, debug=False):
    # Functions needed to provide name of new compress file
    def namer(name):
        folder = os.path.dirname(name)
        return os.path.join(folder, datetime.utcnow().strftime(f'{prefix}%Y-%m-%d.%H%M%S.gz'))

    # Functions needed to created file rotator with gzip compression
    def rotator(source, destination):
        with open(source, "rb") as sf, open(destination, "wb") as df:
            df.write(gzip.compress(sf.read(), 9))
        os.remove(source)

    logger = logging.getLogger('QuasarModelLog')
    logger.setLevel(logging.DEBUG)
    logger.addFilter(ContextFilter())
    formatter = logging.Formatter('%(message)s')
    # Add File handler
    if log_path:
        fh = logging.FileHandler(log_path, 'w')
        fh.setLevel(logging.INFO)
        fh.setFormatter(formatter)
        fh.rotator = rotator
        fh.namer = namer
        logger.addHandler(fh)
    # Add console filter
    if console:
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG if debug else logging.INFO)
        ch.setFormatter(formatter)
        logger.addHandler(ch)
    return logger

# decorator to calculate duration taken by any function.
def time_it(func):
    # added arguments inside the inner1,
    # if function takes any arguments,
    # can be added like this.
    def inner_fnc(*args, **kwargs):
        # storing time before function execution
        begin = time.time()

        results = func(*args, **kwargs)

        # storing time after function execution
        end = time.time()
        logger = logging.getLogger('QuasarModelLog')
        logger.info(f'{func.__name__} took {end - begin:.2f} seconds')
        return results

    return inner_fnc


"""
##################################################
Tools for image data acquisition 
##################################################
"""


def get_analytical_image():
    gauss_vec = GaussList()
    while True:
        try:
            num_comp = int(input('Number of components: '))
            break
        except ValueError:
            print('invalid number, try again')

    for i in range(num_comp):
        gauss = Gaussian()
        while True:
            try:
                gauss.amp = float(input('Enter in amp: '))
                x0, y0 = input('Enter in offsets: ').split(' ')
                gauss.x0, gauss.y0 = int(x0), int(y0)
                sig_x, sig_y = input('Enter in sig_x and sig_y: ').split(' ')
                gauss.a = 1 / float(sig_x) ** 2
                gauss.b = 1 / float(sig_y) ** 2
                gauss.theta = radians(float(input('Enter in rotation (deg): ')))
                break
            except ValueError:
                print('invalid input, try again.')
        gauss_vec.append(gauss)
    image = gauss_vec.build_image()
    gauss_vec.sort()

    return image, gauss_vec


def check_if_test_source(source_path):
    test_sources = {'1': one_source1, '2': one_source2, '3': one_source3, '4': one_source4, '5': one_source5,
                    '6': one_source6, '7': two_sources}
    if source_path in test_sources:
        source_path = test_sources[source_path]
    return source_path


def get_image_from_path(source_path):
    image = fits.getdata(source_path, ext=0)
    image.shape = image.shape[2:]
    return image


"""
##################################################
Tools for running and logging several files
##################################################
"""


def model_several_sources(base_path):
    """
    Model several sources by providing the directory path to where they are located. Will create two .csv files.
    precision.csv stores a summary of the modelled precision of all sources.
    parameters.csv stores the parameters of the modelled gaussians
    :param base_path: path to directory with source images.
    directory structure should be: base_path / source_dir / source.fits
    :return: None
    """
    logger = getLogger('QuasarModelLog')

    header1 = ['Source', 'Mean error', 'Total error', 'mean error as % of max val', 'total error as % of max val']
    header2 = ['Source', 'Amplitude', 'X0', 'Y0', 'Sigma_x', 'Sigma_y', 'theta']

    source_model = SourceModel()
    with open('precisions.csv', 'w', encoding='UTF8') as f1, open('parameters.csv', 'w', encoding='UTF8') as f2:
        writer1 = csv.writer(f1)
        writer1.writerow(header1)

        writer2 = csv.writer(f2)
        writer2.writerow(header2)

        for index, source_path in enumerate([os.path.join(path, name)
                                             for path, _, files in os.walk(base_path) for name in files if name.endswith('map.fits')], 1):
            try:
                image = get_image_from_path(source_path)
                org, mdl, _, gauss_fnd = source_model.process(image)
                precision = source_model.precision(org, mdl)

                data1 = [source_path, precision['mean error'],
                         precision['total error'],
                         precision['mean error as % of max val'],
                         precision['total error as % of max val']]
                writer1.writerow(data1)

                writer2.writerow([source_path])
                for index, gauss in enumerate(gauss_fnd):
                    data2 = [f"Gauss #{index+1}", gauss.amp, gauss.x0, gauss.y0, 1/sqrt(gauss.a), 1/sqrt(gauss.b), gauss.theta]
                    writer2.writerow(data2)

                logger.info(f'{source_path} ok')
            except Exception as exc:
                logger.error(f'{source_path} {str(exc)}')
