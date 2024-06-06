merged_gtf=$1
gffcomapre_tmap=$2

## How many genes have at least one transcript found
echo "Count of Annotated genes with at least one transcript assembled:"
#TODO: Reeplace stringtie with variable 
cat stringtie_merged.gtf | perl -ne 'if ($_ =~ /gene_id\s+\"(\S+)\"\;/){print "$1\n"}' | sort | uniq | wc -l


## How many novel transcripts per gene were detected:
echo "Count of genes with novel transcripts assembled:"
#TODO: replace gffcompare file with tmap 
grep "j" gffcompare.stringtie_merged.gtf.tmap | cut -f 1 | sort | uniq | wc -l


## Display the transcripts that correspond to intergenic regions with the highest read support (candidate novel regions of transcription)
echo "Novel transcripts candidates wit high support"
grep -w "u" gffcompare.stringtie_merged.gtf.tmap | sort -n -k 10 | column -t



