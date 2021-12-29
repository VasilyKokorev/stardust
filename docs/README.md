# Preparing a run
-----
## Translate File
The translate **.bands** file should contain the filter id and the correspinding names for flux and flux uncertainty entries in the catalogue. Please follow [this example](https://github.com/VasilyKokorev/stardust/blob/master/example/stellar%2Bir/example.bands).


Also see here [filters/filters.info](https://github.com/VasilyKokorev/stardust/blob/master/stardust/filters/filters.info) for filter ids.

## Parameters

Stardust requires for each object to have an **ID** and a **redshift**

Additionally it is possible to give the stellar mass (in units of Msol), if one wishes to compute the gas-to-dust ratio and gas mass in the abscence of UV-optical photometry.

If stellar mass is not available in your catalogue the mstar_colname should be given as *None*

Create the params file as shown [here](https://github.com/VasilyKokorev/stardust/blob/master/example/stellar%2Bir/example.param).


## Non-standard photometric bands

The code also allows for a manual input of a wavelength for a given measurement. This can be useful when e.g. some fluxes in your catalogue come from collapsing the spectrum, and do not have a defined filter (e.g. ALMA).

To parse those a separate text file is required, see [example](https://github.com/VasilyKokorev/stardust/blob/master/example/extra_bands/fv.bands_extra):

**Column I:** Wavelength column name.

**Column II:** Flux column name.

**Column III:** Error column name.

**<ins>Wavelength should be given in microns.</ins>**

This creates a square-wave, centred at the defined wavelegnth with a width of 10 microns.

Do not forget to set **EXTRA_BANDS=1** in your config file.

# Preparing the input catalogue: 
-----

The catalogue should contain all the columns that you listed in the ***.param** and **.bands** files.
Ensure that the all flux columns have consistent units.

In principle any astropy.Table readable table format, with column names is parsable. It is however strongly encouraged to use .fits files.

Missing photometry/redshift etc. should have a value of **-99**.

# Creating a config file: 
-----
A config file contains paths to your data, the translate and param files, as well as the output folders. This is also a place where you can specify the units of the catalogue. See the example [here](https://github.com/VasilyKokorev/stardust/blob/master/example/example.config). You can also change the config parameters manually when you load in the module.

# Starting a run: 
-----

See a [jupyter notebook example](https://github.com/VasilyKokorev/stardust/blob/master/example/Stardust_Example.ipynb) on how to fit single objects/catalogues with stardust, and display the outputs.
