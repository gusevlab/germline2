# GREMLIN

## Usage

```
make
g2 [options] <haps file> <sample file> <genetic map file> <output file>
```

Required:

| haps/sample | SHAPEIT/IMPUTE format input phased haplotypes |
| genetic map file | rows = <physical position> <cm/Mb> <cM> |

Options:

| -b | Binary output for large files, see parse_bmatch [default = 0/off] |
| -d | Dynamic hash seed cutoff (for big N) [default = 0/off] |
| -f | Minimum minor allele frequency [default = 0.0] |
| -g | Allowed gaps between seeds [default = 1] |
| -m | Minimum match length [default = 1.0] |
| -s | Skip words with (seeds/samples) less than than this value (for big N) [default = 0.0] |

## Output

| ID1 | ID2 | P0 | P1 | cM | # words |

## Example

`make test` runs sample data in the `example/` directory using the following command:
```
-m 0.9 example/SIM.NE_20000.MATCH_FREQ.SHAPEIT.haps example/SIM.NE_20000.MATCH_FREQ.SHAPEIT.sample example/genMap.1KG.b37.chr1.map example/SIM.NE_20000.MATCH_FREQ.INFERRED.match
```

The output segments are then evaluated for accuracy using the `example/accuracy.sh` script.

This data was simulated using the ARGON software as shown in `example/sim.sh`, down-sampled to a HapMap3 allele frequency distribution, and phased with SHAPEIT2.
