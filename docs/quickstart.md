# Defining band names and parameters: 

Bands:
-----

Create a text file that contains the following information:

**Column I:** Filter **ID** as given in the [filters/filters.info](https://github.com/VasilyKokorev/ctf/blob/master/filters/filters.info) file.

**Column II:** Flux column name.

**Column III:** Error column name.

**<ins>The input flux and flux errors must be in mJy!</ins>**

(this would be changed in future versions)


Use this [template](https://github.com/VasilyKokorev/ctf/blob/master/example/example.bands) as an example for the input bands.

Any text format that can be read by numpy.loadtxt should work.


Params:
-----

The software requires for a given object to have an **ID** and a **redshift**. 

Additionally it is possible to give the stellar mass (in units of Msol), if one wishes to compute the gas-to-dust ratio and gas mass.

If stellar mass is not available in your catalogue the mstar_colname should be given as *None*

Create the params file as shown [here](https://github.com/VasilyKokorev/ctf/blob/master/example/example.param).


Non-standard photometric bands:
-----

The code also allows for a manual input of a wavelength for a given measurement. This can be useful when e.g. some of the objects in your catalogue have ALMA photometry available for them.

To parse those a separate text file is required, see [example](https://github.com/VasilyKokorev/ctf/blob/master/example/example.bands_extra):

**Column I:** Wavelength column name.

**Column II:** Flux column name.

**Column III:** Error column name.

**<ins>Wavelength should be given in microns.</ins>**

This creates a square-wave, centred at the defined wavelegnth with a width of 10 microns.

# Preparing the input catalogue: 

Currently the only accepted format is a standard **astropy .fits Table**


 
