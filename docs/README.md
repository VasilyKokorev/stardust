# Preparing a run
-----
## Translate File
Create a translate file between the photometry and stardust filter file. Should contain the filter id and the correspinding names for flux and flux uncertainty entries in the catalogue. Follow [this example](https://github.com/VasilyKokorev/stardust/blob/master/example/stellar%2Bir/example.bands).


Also see here [filters/filters.info](https://github.com/VasilyKokorev/ctf/blob/master/filters/filters.info) for filter ids.

## Parameters

Stardust requires for each object to have an **ID** and a **redshift**

Additionally it is possible to give the stellar mass (in units of Msol), if one wishes to compute the gas-to-dust ratio and gas mass in the abscence of UV-optical photometry.

If stellar mass is not available in your catalogue the mstar_colname should be given as *None*

Create the params file as shown [here](https://github.com/VasilyKokorev/stardust/blob/master/example/stellar%2Bir/example.param).


## Non-standard photometric bands

The code also allows for a manual input of a wavelength for a given measurement. This can be useful when e.g. some fluxes in your catalogue come from collapsing the spectrum, and do not have a defined filter (e.g. ALMA).

To parse those a separate text file is required, see [example](https://github.com/VasilyKokorev/stardust/blob/master/example/ir_only/example.bands_extra):

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


# Starting a run: 
-----

(Very preliminary)

Navigate to the [init.py](https://github.com/VasilyKokorev/ctf/blob/master/config.py) and provide the paths for the input catalogue, bands and param files, as well as a path to the output directory.

Run ``main.py`` with **Python 3.6** or above.
