import logging

# Standard output handler
stdout_handler = logging.StreamHandler()
stdout_handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(filename)s %(levelname)s %(message)s')
stdout_handler.setFormatter(formatter)

# File handler
file_handler = logging.FileHandler('dear_diary.log')
file_handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s %(filename)s %(levelname)s %(message)s')
file_handler.setFormatter(formatter)

# Custom logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
logger.addHandler(stdout_handler)
logger.addHandler(file_handler)