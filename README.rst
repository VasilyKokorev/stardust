``Stardust``: Composite Template Fitting Software
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This repository provides an SED fitting tool, primarily focused on the IR part of the spectrum.


Requirements: 
~~~~~~~~~~~~~
    .. code:: 
    
       python==3.6.10
       numpy
       scipy
       matplotlib
       astropy
       multiprocess
       tqdm
       
Installation:
~~~~~~~~~~~~~
It is recommended to create a separate (**Python==3.6.10**) environment with miniconda. See instructions 'here, <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>'
    .. code:: bash
    
        $ git git@github.com:VasilyKokorev/stardust.git
        $ cd stardust
        $ pip install .
  
Usage:
~~~~~~
The primary use of this code is to extract the parameters from the infrared photometry. 
The program also has an ability to fit AGN and Stellar emission templates, if the user desires to do so.


An example data-set, a subset from the `COSMOS Super-Deblended Catalogue 2 (Jin+18), <https://ui.adsabs.harvard.edu/abs/2018ApJ...864...56J/abstract>`__ along with the example configuration files are provided in the `example folder <https://github.com/VasilyKokorev/ctf/tree/master/example>`__.

See the `quickstart guide <https://github.com/VasilyKokorev/ctf/blob/master/docs/README.md>`__ and `docs folder <https://github.com/VasilyKokorev/ctf/tree/master/docs>`__ for more detailed instructions.

Contacts:
~~~~~~

Vasily Kokorev: vasily.kokorev.astro@gmail.com
