#!/usr/bin/env python


import logging
LOG_FILENAME = '/tmp/logging_example.out'
logging.basicConfig(filename=LOG_FILENAME,level=logging.DEBUG,filemode="w")

logging.debug('This message should go to the log file')
logging.info('tester')
