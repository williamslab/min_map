# Minimal viable genetic map

Given an input genetic map and a maximum error rate, greedily drops positions
to produce a map with fewer entries. Following processing, the only remaining
map positions are those that, if dropped, would produce an error greater than
the threshold (assuming linear interpolation for finding genetic positions).

## Example

To run on chromosome 17 of the [HapMap II genetic map](ftp://ftp.ncbi.nlm.nih.gov/hapmap/recombination/) ([http map link](http://bit.ly/33s6XQG)), with 0.02 cM of error tolerated, run:

    ./min_map.py -error 0.02 -chr 17 -mapfile /path/to/genetic_map_GRCh37_chr17.txt

This produces `min_viable_map17.txt` with the subsetted genetic map. (Use
`-out` to change the output prefix; default is `min_viable_map`.)

## Options

Use `./min_map.py -h` for a list of options. If using a genetic map that
doesn't match the format of the HapMap map, be sure to set `-physcol` for the
physical position column and `-genetcol` for the genetic position column.

By default, the script assumes the input map has a header, which it deletes.
Use `-noheader` to treat all the input lines (except those that begin with
`#`) as map entries.

To run on a sex-specific genetic map, combine the male and female maps into
separate columns in a single input file (for each chromosome). You can then
call `./min_map.py` with `-sexspecific` and set `-genetcol` to the column with
one of the two genetic positions (male or female) and `-genet2col` to the
other.
