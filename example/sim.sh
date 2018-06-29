#!/bin/sh

POP=$1
OUT="SIM.NE_$POP"

java -Xmx18G -jar ARGON.0.1.jar -N $POP -pop $p -map genMap.1KG.b37.chr1.map -size 10 -out $OUT -seq 0.01 -haps-out -IBD 1 -seed $SEED

mv $OUT.ibd $OUT.TRUE.ibd

cat $OUT.TRUE.ibd \
| awk '{ FAM1=int(($1-1)/2+1); SUFF1=($1-1)%2; FAM2=int(($2-1)/2+1); SUFF2=($2-1)%2; print FAM1,FAM2,$3,$4,$7 }' \
| sort -g -k1 -k2 -k3 | awk -f stitch.awk | awk 'NF != 0' | tr ' ' '\t' > $OUT.TRUE.match
