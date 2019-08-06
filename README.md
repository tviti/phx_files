# PHX_files

A set of `python` functions for reading `.TSn` and `.TBL` files produced by a Phoenix Geophysics brand `MTU-5A`.

This library may work with data files produced by other Phoenix brand receivers, but is only tested with `MTU-5A` data.
Note that the `MTU-5A` data comes in 24 bit 2's complement, whereas the newer receivers (e.g. the `MTU-5C`) are 32 bit (I think...), hence this package won't work with those data files (I don't think the changes to make it work should be too difficult though, but I don't have any data to try it out on).

I am also trying to include utilities that I have developed during my time working with Phoenix data, that I think might be useful for other people.

## Latitude and Longtitude formats

The format for the `LATG` and `LNGG` fields in the `.TBL` files is a little confusing (well... it confused me at first).
They are given as `degrees` and `decimal minutes`, with no seperator between the two i.e.

`LATG = (DD MM.MM)`

and

`LNGG = (DDD MM.MM)`

## Datetime formats

The format for the `STIM` and `ETIM` fields is as follows

`val = (sec, min, hour, day, month, year, day of week, century)`
