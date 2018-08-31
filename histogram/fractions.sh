#!/bin/bash

grep -v '#' hist_600 | \
  awk '{total += $3; metal += $4; parts += $7; m1 += $8; m2 += $10; m3 += $12; m4 += $14; m5 += $16; m6 += $18; m7 += $20} END \
       {avg = (m1+m2+m3+m4+m5+m6+m7)/7.0; var = ((m1-avg)^2 + (m2-avg)^2 + (m3-avg)^2 + (m4-avg)^2 + (m5-avg)^2 + (m6-avg)^2 + (m7-avg)^2)/7.0; \
       print "Total mass: ", total/2e33, "    DM/Baryon ratio: ", 4.55e7/(total/2e33); \
       print  "Total metals: ", metal/2e33, "    Fraction within virial radius: ", metal/2e33/16.7; \
       print  "Total particles: ", parts/2e33, "    Fraction within virial radius: ", parts/2e33/16.7; \
       print  "SN mass within virial radius: ", m1/2e33, m2/2e33, m3/2e33, m4/2e33, m5/2e33, m6/2e33, m7/2e33; \
       print  "Average mass +/- error: ", avg/2e33, sqrt(var)/2e33
       print  "SN fractions within virial radius: ", m1/2e33/2.5, m2/2e33/2.5, m3/2e33/2.5, m4/2e33/2.5, m5/2e33/2.5, m6/2e33/2.35, m7/2e33/1.85 }'
