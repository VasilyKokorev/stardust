from setuptools import setup
from setuptools.extension import Extension


setup(
    name = "stardust",
    author = "Vasily Kokorev",
    author_email = "vasily.kokorev.astro@gmail.com",
    description = "Composite Template SED Fitting Tool",
    license = "MIT",
    url = "https://github.com/vasilykokorev/stardust",  
    packages=['stardust'],
    scripts=['stardust/utils.py','stardust/filter_importer.py'],
    install_requires = ['numpy',
                        'scipy',
                        'matplotlib',
                        'astropy',
                        'multiprocess',
                        'tqdm'],

    package_data={'stardust': ['filters/*', 'templates/*',
                           'templates/AGN/*','templates/QSO/*',
                           'templates/DL07/*']},
)