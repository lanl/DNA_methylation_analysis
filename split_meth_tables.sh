#!/bin/bash
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$4,$5}' ../Calls/7-26/CpG_result_table.txt >../Calls/7-26-Aza1_CpG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$4,$5}' ../Calls/7-26/CHG_result_table.txt >../Calls/7-26-Aza1_CHG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$4,$5}' ../Calls/7-26/CHH_result_table.txt >../Calls/7-26-Aza1_CHH.txt
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$4,$5}' ../Calls/7-27/CpG_result_table.txt >../Calls/7-27-Aza1_CpG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$4,$5}' ../Calls/7-27/CHG_result_table.txt >../Calls/7-27-Aza1_CHG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$4,$5}' ../Calls/7-27/CHH_result_table.txt >../Calls/7-27-Aza1_CHH.txt
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$4,$5}' ../Calls/7-28/CpG_result_table.txt >../Calls/7-28-Aza1_CpG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$4,$5}' ../Calls/7-28/CHG_result_table.txt >../Calls/7-28-Aza1_CHG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$4,$5}' ../Calls/7-28/CHH_result_table.txt >../Calls/7-28-Aza1_CHH.txt
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$4,$5}' ../Calls/7-31/CpG_result_table.txt >../Calls/7-31-Aza2_CpG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$4,$5}' ../Calls/7-31/CHG_result_table.txt >../Calls/7-31-Aza2_CHG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$4,$5}' ../Calls/7-31/CHH_result_table.txt >../Calls/7-31-Aza2_CHH.txt
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$6,$7}' ../Calls/7-25/CpG_result_table.txt >../Calls/7-25-Aza1_CpG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$6,$7}' ../Calls/7-25/CHG_result_table.txt >../Calls/7-25-Aza1_CHG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$6,$7}' ../Calls/7-25/CHH_result_table.txt >../Calls/7-25-Aza1_CHH.txt
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$6,$7}' ../Calls/7-26/CpG_result_table.txt >../Calls/7-26-Aza2_CpG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$6,$7}' ../Calls/7-26/CHG_result_table.txt >../Calls/7-26-Aza2_CHG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$6,$7}' ../Calls/7-26/CHH_result_table.txt >../Calls/7-26-Aza2_CHH.txt
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$6,$7}' ../Calls/7-27/CpG_result_table.txt >../Calls/7-27-Aza2_CpG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$6,$7}' ../Calls/7-27/CHG_result_table.txt >../Calls/7-27-Aza2_CHG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$6,$7}' ../Calls/7-27/CHH_result_table.txt >../Calls/7-27-Aza2_CHH.txt
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$6,$7}' ../Calls/7-28/CpG_result_table.txt >../Calls/7-28-Aza2_CpG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$6,$7}' ../Calls/7-28/CHG_result_table.txt >../Calls/7-28-Aza2_CHG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$6,$7}' ../Calls/7-28/CHH_result_table.txt >../Calls/7-28-Aza2_CHH.txt
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$6,$7}' ../Calls/7-31/CpG_result_table.txt >../Calls/7-31-Aza3_CpG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$6,$7}' ../Calls/7-31/CHG_result_table.txt >../Calls/7-31-Aza3_CHG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$6,$7}' ../Calls/7-31/CHH_result_table.txt >../Calls/7-31-Aza3_CHH.txt
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$8,$9}' ../Calls/7-25/CpG_result_table.txt >../Calls/7-25-Aza2_CpG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$8,$9}' ../Calls/7-25/CHG_result_table.txt >../Calls/7-25-Aza2_CHG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$8,$9}' ../Calls/7-25/CHH_result_table.txt >../Calls/7-25-Aza2_CHH.txt
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$8,$9}' ../Calls/7-26/CpG_result_table.txt >../Calls/7-26-Aza3_CpG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$8,$9}' ../Calls/7-26/CHG_result_table.txt >../Calls/7-26-Aza3_CHG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$8,$9}' ../Calls/7-26/CHH_result_table.txt >../Calls/7-26-Aza3_CHH.txt
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$8,$9}' ../Calls/7-27/CpG_result_table.txt >../Calls/7-27-Aza3_CpG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$8,$9}' ../Calls/7-27/CHG_result_table.txt >../Calls/7-27-Aza3_CHG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$8,$9}' ../Calls/7-27/CHH_result_table.txt >../Calls/7-27-Aza3_CHH.txt
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$8,$9}' ../Calls/7-28/CpG_result_table.txt >../Calls/7-28-Aza3_CpG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$8,$9}' ../Calls/7-28/CHG_result_table.txt >../Calls/7-28-Aza3_CHG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$8,$9}' ../Calls/7-28/CHH_result_table.txt >../Calls/7-28-Aza3_CHH.txt
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$8,$9}' ../Calls/7-31/CpG_result_table.txt >../Calls/7-31-Aza1_CpG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$8,$9}' ../Calls/7-31/CHG_result_table.txt >../Calls/7-31-Aza1_CHG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($5>50 && $7>50 && $9>50){print $1,$2,$3,$8,$9}' ../Calls/7-31/CHH_result_table.txt >../Calls/7-31-Aza1_CHH.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$10,$11}' ../Calls/7-25/CpG_result_table.txt >../Calls/7-25-Con1_CpG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$10,$11}' ../Calls/7-25/CHG_result_table.txt >../Calls/7-25-Con1_CHG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$10,$11}' ../Calls/7-25/CHH_result_table.txt >../Calls/7-25-Con1_CHH.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$10,$11}' ../Calls/7-26/CpG_result_table.txt >../Calls/7-26-Con2_CpG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$10,$11}' ../Calls/7-26/CHG_result_table.txt >../Calls/7-26-Con2_CHG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$10,$11}' ../Calls/7-26/CHH_result_table.txt >../Calls/7-26-Con2_CHH.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$10,$11}' ../Calls/7-27/CpG_result_table.txt >../Calls/7-27-Con3_CpG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$10,$11}' ../Calls/7-27/CHG_result_table.txt >../Calls/7-27-Con3_CHG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$10,$11}' ../Calls/7-27/CHH_result_table.txt >../Calls/7-27-Con3_CHH.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$10,$11}' ../Calls/7-28/CpG_result_table.txt >../Calls/7-28-Con2_CpG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$10,$11}' ../Calls/7-28/CHG_result_table.txt >../Calls/7-28-Con2_CHG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$10,$11}' ../Calls/7-28/CHH_result_table.txt >../Calls/7-28-Con2_CHH.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$10,$11}' ../Calls/7-31/CpG_result_table.txt >../Calls/7-31-Con2_CpG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$10,$11}' ../Calls/7-31/CHG_result_table.txt >../Calls/7-31-Con2_CHG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$10,$11}' ../Calls/7-31/CHH_result_table.txt >../Calls/7-31-Con2_CHH.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$12,$13}' ../Calls/7-25/CpG_result_table.txt >../Calls/7-25-Con2_CpG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$12,$13}' ../Calls/7-25/CHG_result_table.txt >../Calls/7-25-Con2_CHG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$12,$13}' ../Calls/7-25/CHH_result_table.txt >../Calls/7-25-Con2_CHH.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$12,$13}' ../Calls/7-26/CpG_result_table.txt >../Calls/7-26-Con3_CpG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$12,$13}' ../Calls/7-26/CHG_result_table.txt >../Calls/7-26-Con3_CHG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$12,$13}' ../Calls/7-26/CHH_result_table.txt >../Calls/7-26-Con3_CHH.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$12,$13}' ../Calls/7-27/CpG_result_table.txt >../Calls/7-27-Con1_CpG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$12,$13}' ../Calls/7-27/CHG_result_table.txt >../Calls/7-27-Con1_CHG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$12,$13}' ../Calls/7-27/CHH_result_table.txt >../Calls/7-27-Con1_CHH.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$12,$13}' ../Calls/7-28/CpG_result_table.txt >../Calls/7-28-Con1_CpG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$12,$13}' ../Calls/7-28/CHG_result_table.txt >../Calls/7-28-Con1_CHG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$12,$13}' ../Calls/7-28/CHH_result_table.txt >../Calls/7-28-Con1_CHH.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$12,$13}' ../Calls/7-31/CpG_result_table.txt >../Calls/7-31-Con3_CpG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$12,$13}' ../Calls/7-31/CHG_result_table.txt >../Calls/7-31-Con3_CHG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$12,$13}' ../Calls/7-31/CHH_result_table.txt >../Calls/7-31-Con3_CHH.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$14,$15}' ../Calls/7-25/CpG_result_table.txt >../Calls/7-25-Con3_CpG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$14,$15}' ../Calls/7-25/CHG_result_table.txt >../Calls/7-25-Con3_CHG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$14,$15}' ../Calls/7-25/CHH_result_table.txt >../Calls/7-25-Con3_CHH.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$14,$15}' ../Calls/7-26/CpG_result_table.txt >../Calls/7-26-Con1_CpG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$14,$15}' ../Calls/7-26/CHG_result_table.txt >../Calls/7-26-Con1_CHG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$14,$15}' ../Calls/7-26/CHH_result_table.txt >../Calls/7-26-Con1_CHH.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$14,$15}' ../Calls/7-27/CpG_result_table.txt >../Calls/7-27-Con2_CpG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$14,$15}' ../Calls/7-27/CHG_result_table.txt >../Calls/7-27-Con2_CHG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$14,$15}' ../Calls/7-27/CHH_result_table.txt >../Calls/7-27-Con2_CHH.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$14,$15}' ../Calls/7-28/CpG_result_table.txt >../Calls/7-28-Con3_CpG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$14,$15}' ../Calls/7-28/CHG_result_table.txt >../Calls/7-28-Con3_CHG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$14,$15}' ../Calls/7-28/CHH_result_table.txt >../Calls/7-28-Con3_CHH.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$14,$15}' ../Calls/7-31/CpG_result_table.txt >../Calls/7-31-Con1_CpG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$14,$15}' ../Calls/7-31/CHG_result_table.txt >../Calls/7-31-Con1_CHG.txt
awk 'BEGIN{FS="\t";OFS="\t"};($11>50 && $13>50 && $15>50){print $1,$2,$3,$14,$15}' ../Calls/7-31/CHH_result_table.txt >../Calls/7-31-Con1_CHH.txt