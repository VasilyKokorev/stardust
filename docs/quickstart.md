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

Additionally it is possible to give the stellar mass (in units of Msol), if one wishes to compute the gas-to-dust ratio and gas mass.

If stellar mass is not available in your catalogue the mstar_colname should be given as *None*

Create the params file as shown [here](https://github.com/VasilyKokorev/ctf/blob/master/example/example.param)
 
