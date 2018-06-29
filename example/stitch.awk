{
if($1==pre_a && $2 == pre_b && $3 - pre_end < 1e3) { if($3>pre_end && $4>pre_end) {pre_end=$4; len += $5} else if ($4 > pre_end) { len += $5/($4-$3) * ($4-pre_end); pre_end=$4; } } else { print pre_a, pre_b, pre_start, pre_end, len; pre_a = $1; pre_b=$2; pre_start=$3; pre_end=$4; len=$5 }
}
