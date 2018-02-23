#import packages
import sys
import PolyLibScan
import logging
import time

# Set logging config
logging.basicConfig(format='[%(asctime)s] %(message)s', level=logging.INFO, datefmt='%Y/%m/%d-%H:%M:%S')


def main(config):
	logging.info('Starting Job..')

	j = PolyLibScan.Job(config)

	j.setup_env()
	logging.info('Finished setting up environment.')

	j.run()
	time.sleep(5)
	logging.info('Finished simulations. Starting to clean up.')
	j.setup_job_save()
	j.save()
	logging.info('Finished saving to DB.')
	j.clean_up()
	logging.info('Finished cleaning up.')
	logging.info('Finished LAMMPS Job.')
	

if __name__ == '__main__':
    config = sys.argv[1]
    main(config)
