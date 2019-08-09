# PHX_files

A set of `python` functions for reading `.TSn` and `.TBL` files produced by a Phoenix Geophysics brand `MTU-5A`.

This library may work with data files produced by other Phoenix brand receivers, but is only tested with `MTU-5A` data.

Note that the `MTU-5A` data comes in 24 bit 2's complement, whereas the newer receivers (e.g. the `MTU-5C`) are 32 bit (I think...), hence this package won't work with those data files (I don't think the changes to make it work should be too difficult though, but I don't have any data to try it out on).
If anybody is interested in seeing this library updated to work with the 32 bit data (and is willing to furnish data files for testing), do let me know.

Given a path to the Phoenix utility `SysCal.exe` (along with the directories containing the requisite `.clb` and `.clc` box and sensor calibration files), the function `get_syscal` can generate system response functions that can be used for calibration.
Note that prior to performing the deconvolution, the data first should be transformed to unit scale (i.e. scaled by 1.0/2**23 for 24 bit data).
`get_syscal` is also currently only written to work on a Mac or Linux machine that has Wine installed, but could be easily modified to work on a windows machine as well (lmk if that is something that you’d be interested in).
Also important for calibration: the Phoenix `MTC-150` coils (and possibly their other models as well) actually record the NEGATIVE of the magnetic field!
So, you need to either take the negative of the H-field channels prior to processing, or instruct your processing routine to rotate the H-field data by 180 degrees.

I am also trying to include utilities that I have developed during my time working with Phoenix data, that I think might be useful for other people.
These can be found in the `examples` dir.

## Latitude and Longtitude formats

The format for the `LATG` and `LNGG` fields in the `.TBL` files is a little confusing (well... it confused me at first).
They are given as `degrees` and `decimal minutes`, with no seperator between the two i.e.

`LATG = (DD MM.MM)`

and

`LNGG = (DDD MM.MM)`

## Datetime formats

The format for the `STIM` and `ETIM` fields is as follows

`val = (sec, min, hour, day, month, year, day of week, century)`
