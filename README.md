# PHX_files

A set of `python` functions for reading `.TSn` and `.TBL` files produced by a Phoenix Geophysics brand `MTU-5A`.

This library may work with data files produced by other Phoenix brand receivers, but is only tested with `MTU-5A` data.

## Latitude and Longtitude formats

The format for the `LATG` and `LNGG` fields in the `.TBL` files is a little confusing (well... it confused me at first).
They are given as `degrees` and `decimal minutes`, with no seperator between the two i.e.

`LATG = (DD MM.MM)`

and

`LNGG = (DDD MM.MM)`

## Datetime formats

The format for the `STIM` and `ETIM` fields is as follows

`val = (sec, min, hour, day, month, year, day of week, century)`
