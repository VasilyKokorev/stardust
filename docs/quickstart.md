# Defining band names and parameters: 

Bands:
-----

Create a text file that contains the following information:

**Column I:** Filter *ID* as given in the [filters/filters.info](https://github.com/VasilyKokorev/ctf/blob/master/filters/filters.info) file.

**Column II:** Name of the flux array corresponding to that filter in your input table.

**Column III:** Error array corresponding to the flux array.

**<ins>The input flux and flux errors must be in mJy!</ins>**


Use the template below as an example for the input bands:


18 FCH1 DFCH1

19 FCH2 DFCH2

20 FCH3 DFCH3

21 FCH4 DFCH4

325 F24 DF24

329 F100 DF100

330 F160 DF160

331 F250 DF250

332 F350 DF350

333 F500 DF500

324 F850 DF850

348 F1100 DF1100

350 F1200 DF1200


Params:
-----

Currently the code needs the galaxy id and the redshift in order to fit a source.
Stellar mass Mstar is not required, but if available it can be used to compute the gas-to-dust ratio.

Stellar mass in your catalogue should be in units of Msol.

Define the params file as shown below:

Column I: Internal naming convention. Do not change.
Column II: Names for these params in your input table.

If stellar mass is not available in your catalogue, write None, in the second column, instead of Mstar. 


  
