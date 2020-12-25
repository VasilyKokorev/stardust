# Defining band names and parameters: 

Bands:
-----

Create a text file that contains the following information:

**Column I:** Filter *ID* as given in the [filters/filters.info](https://github.com/VasilyKokorev/ctf/blob/master/filters/filters.info) file

**Column II:** Name of the flux array corresponding to that filter in your input table

**Column III:** Error array

**<ins>The input flux and flux errors must be in mJy!</ins>**


Use the template below as an example.
Params:

Currently the code needs the galaxy id and the redshift in order to fit a source.
Stellar mass Mstar is not required, but if available it can be used to compute the gas-to-dust ratio.

Stellar mass in your catalogue should be in units of Msol.

Define the params file as shown below:

Column I: Internal naming convention. Do not change.
Column II: Names for these params in your input table.

If stellar mass is not available in your catalogue, write None, in the second column, instead of Mstar. 


  
