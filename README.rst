``Stardust``: Composite Template Fitting Software
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**Stardust** is a Python package which is designed to extract galaxy properties by fitting their multiwavelength data to a set of linearly combined templates. 

This algorithm brings together three different families of templates:

1) UV+Optical emission from dust unobscured stellar light

2) AGN heated dust in the MIR

3) IR dust reprocessed stellar light in the NIR-FIR

One of the features of **Stardust** is that the template fitting does not rely on energy balance. As a result the total luminosity of dust obscured and dust unobscured stellar light do not rely on each other, and it is possible to fit object (e.g. SMGs) where the energy balance approach might not be necessarily applicable. Detailed description of **Stardust** and its first application are presented in  `Kokorev et al. 2021 <https://ui.adsabs.harvard.edu/abs/2021ApJ...921...40K/abstract>`__.



Requirements: 
~~~~~~~~~~~~~
    .. code:: 
    
       python>=3.9.0
       numpy
       scipy
       matplotlib
       astropy
       multiprocess
       tqdm
       
Installation:
~~~~~~~~~~~~~
For MacOS users it is recommended to create a separate **Python 3.6.10**
environment with miniconda. See instructions `here <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`__.

You can then pull stardust and install it with pip.

    .. code:: bash
    
        $ git clone https://github.com/VasilyKokorev/stardust.git
        $ cd stardust
        $ pip install .
  
Usage:
~~~~~~
The primary use of this code is to extract the parameters from the infrared photometry. 
The program also has an ability to fit AGN and Stellar emission templates, if the user desires to do so.
It is recommended to use this code within a jupyter notebook.

See `here <https://github.com/VasilyKokorev/stardust/blob/master/example/Stardust_Example.ipynb>`__ for a quick example of how to fit a catalogue. 


An example data-set, a subset from the `COSMOS Super-Deblended Catalogue 2 (Jin+18), <https://ui.adsabs.harvard.edu/abs/2018ApJ...864...56J/abstract>`__ along with the example configuration files are provided in the `example folder <https://github.com/VasilyKokorev/ctf/tree/master/example>`__.

See the `quickstart guide <https://github.com/VasilyKokorev/ctf/blob/master/docs/README.md>`__ and `docs folder <https://github.com/VasilyKokorev/ctf/tree/master/docs>`__ for more detailed instructions.

Acknowledgements:
~~~~~~~~~~~~~
Filters curves are adopted from `eazy-py <https://github.com/gbrammer/eazy-py>`__ by Gabe Brammer.

Troubleshooting:
~~~~~~
**1. The module does not install:** If you are on linux, try using Python 3.9.1 instead.

Contacts:
~~~~~~

Vasily Kokorev: vasily.kokorev.astro@gmail.com
