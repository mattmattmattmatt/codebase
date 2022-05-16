for f in `cat SRA_list`; do echo $f; ~/software/suites/sratoolkit.2.9.6-1-ubuntu64/bin/fastq-dump --gzip --split-3 $f; done
