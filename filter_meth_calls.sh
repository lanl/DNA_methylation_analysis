#!/usr/bin/bash
#Script to cull Pico methylation calls with read depth less than 20
#creates two files:XYZ_agt20.txt and XYZ_not_agt20.txt
awk 'BEGIN{FS="\t"};NR==1||($5>19&&$7>19&&$9>19&&$11>19&&$13>19&&$15>19){print $0}' $1 > $1"_agt20.txt"
awk 'BEGIN{FS="\t"};NR==1||!($5>19&&$7>19&&$9>19&&$11>19&&$13>19&&$15>19){print $0}' $1 >$1"_not_agt20.txt"
