#import packages
import sys
import PolyLibScan
import logging
import time

# Set logging config
logging.basicConfig(format='[%(asctime)s] %(message)s', level=logging.INFO, datefmt='%Y/%m/%d-%H:%M:%S')

if __name__ == '__main__':
    config = sys.argv[1]
    pls.job.run_with_config(config)
