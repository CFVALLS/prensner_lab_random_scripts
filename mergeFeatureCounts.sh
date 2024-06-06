#!/usr/bin/env bash
# get the count
ls -1  *geneCounts.txt | parallel 'cat {} | sed '1d' | cut -f7 {} > {/.}_clean.txt' 
ls -1  *geneCounts.txt | head -1 | xargs cut -f1 > genes.txt
paste genes.txt *geneCounts_clean.txt > sampleCount.txt
