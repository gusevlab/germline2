IBD_TRU=$1
IBD_EST=$2
OUT=$3

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# increase inferred segment length and measure false %
cat \
<(cat $IBD_EST | awk -v l=0 '($5) > l && $1 != $2 { if($2 > $1) id=$1":"$2; else id=$2":"$1; print 1,id,$3,$4,$5 }') \
<(cat $IBD_TRU | awk -v l=0 '($5) > l && $1 != $2 { if($2 > $1) id=$1":"$2; else id=$2":"$1; print 0,id,$3,$4,$5 }') \
| sort -k2,2 -k1,1 | $DIR/accuracy | awk -v p=$OUT '{ print $0 > p".accuracy."$1".out" }'

for l in `seq 1 10`; do
fp=`cat $OUT.accuracy.1.out | awk -v l=$l '$NF > l { print 1 - $2 }' | awk -f $DIR/avg.awk`
echo "FP cM>$l : $fp" | tr ' ' '\t'
done
echo "---"

for l in `seq 1 10`; do
fp=`cat $OUT.accuracy.0.out | awk -v l=$l '$NF > l { print $2 }' | awk -f $DIR/avg.awk`
echo "TP cM>$l : $fp" | tr ' ' '\t'
done
echo "---"

rm $OUT.accuracy.*
