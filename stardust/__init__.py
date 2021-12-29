import sys

if (sys.version_info >= (3, 7)) | (sys.version_info <= (3, 5)):
    print(f'Detected Python {sys.version[:6]}, some functionality might be broken. Please use Python 3.6.10')

from . import main
from . import utils