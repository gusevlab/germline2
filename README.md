# :herb: GREMLIN

Efficiently identifying shared genetic segments in large-scale data.

## Usage

The `boost/1.57.0` library is required.

```
make
g2 [options] <haps file> <sample file> <genetic map file> <output file>
```

Required inputs:

| field | description |
| --- | --- |
| haps/sample | SHAPEIT/IMPUTE format input phased haplotypes, alleles can be anything as long as haplotype entries are 0/1 |
| genetic map file | Each row has three fields: [`physical position`] [`cm/Mb`] [`cM`], and the 2nd field is ignored |
| output file | Pointer to where the outputs will go, will generate an $OUT.match file |

Optional switches:

| switch | description |
| --- | --- |
| -b | Binary output for large files, see parse_bmatch [default off] |
| -d | Dynamic hash seed cutoff (for big N) [default = 0/off] |
| -f | Minimum minor allele frequency [default = 0.0] |
| -g | Allowed gaps between seeds [default = 1] |
| -h | Haploid mode, do not allow switches between haplotypes [default off] |
| -m | Minimum match length [default = 1.0] |
| -s | Skip words with (seeds/samples) less than than this value (for big N) [default = 0.0] |

## Output

Output goes into a $OUT.match file with each row containing the following entries:

| ID1 | ID2 | P0 | P1 | cM | # words |
| --- | --- | --- | --- | --- | --- |

If haploid mode is on (`-h`) then ".0" or ".1" is appended to the IDs to indicate a match along the first or second haplotype.

For large data, you can enable binary outputs by adding the `-b` switch, which will generate three files (`$OUT.bmatch/bmid/bsid`) that can be parsed using the provided `parse_bmatch` program (~3x reduction in file size).

## Example

`make test` runs sample data in the `example/` directory using the following command:
```
./g2 -m 0.9 \
example/SIM.NE_20000.MATCH_FREQ.SHAPEIT.haps \
example/SIM.NE_20000.MATCH_FREQ.SHAPEIT.sample \
example/genMap.1KG.b37.chr1.map \
example/SIM.NE_20000.MATCH_FREQ.INFERRED.match
```

The output segments are then evaluated for accuracy using the `example/accuracy.sh` script.

This data was simulated using the ARGON software as shown in `example/sim.sh`, down-sampled to a HapMap3 allele frequency distribution, and phased with SHAPEIT2.
