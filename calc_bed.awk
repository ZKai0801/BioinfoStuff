#!/usr/bin/awk -f

{sum+=$3-$2}END{print sum}
