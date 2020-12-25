``CTF``: Composite Template Fitting Software
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This repository provides an SED fitting tool, primarily focused on the IR part of the spectrum.


Requirements: 
~~~~~~~~~~~~~
    .. code:: 
    

       numpy
       scipy
       astropy
       multiprocessing
       
Installation:
~~~~~~~~~~~~~
Not required.

    .. code:: bash
    
        $ git git@github.com:VasilyKokorev/ctf.git
        $ cd ctf
  
Usage:
~~~~~~
The primary use of this code is to extract the parameters from the infrared photometry. 
The program also has an ability to fit AGN and Stellar emission templates, if the user desires to do so.


The example data-set, a subset from the COSMOS Super-Debldended Catalogue 2  `Jin+18 <https://ui.adsabs.harvard.edu/abs/2018ApJ...864...56J/abstract>`__ along with the example configuration files are provided in the _`example <https://github.com/VasilyKokorev/ctf/tree/master/example>`__ folder.

See the /docs folder for detailed instructions.
