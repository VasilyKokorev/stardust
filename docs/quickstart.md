# Defining band names and parameters: 

Bands:
-----

Create a text file that contains the following information:

**Column I:** Filter *ID* as given in the [filters/filters.info](https://github.com/VasilyKokorev/ctf/blob/master/filters/filters.info) file.

**Column II:** Name of the flux array corresponding to that filter in your input table.

**Column III:** Error array corresponding to the flux array.

**<ins>The input flux and flux errors must be in mJy!</ins>**

(this would be changed in future versions)


Use this [template](https://github.com/VasilyKokorev/ctf/blob/master/example/example.bands) as an example for the input bands.

Any text format that can be read by numpy.loadtxt should work.


Params:
-----

The software requires for a given object to have an **ID** and a **redshift**. 



Stellar mass in your catalogue should be in units of Msol.

Define the params file as shown below:

Column I: Internal naming convention. Do not change.
Column II: Names for these params in your input table.

If stellar mass is not available in your catalogue, write None, in the second column, instead of Mstar. 


  
