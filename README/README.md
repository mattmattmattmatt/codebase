
# GATK Calls

1. Prepare unmapped bams according to GATK best practice workflow
2. Map reads using bwa according to GATK best practices
3. Call variants to produce g.vcf files
```bash
java -Xmx16G -Xms16g -jar $GATK_HOME/GenomeAnalysisTK.jar \
    	-T HaplotypeCaller \
    	-I $bam \
    	-R $fasta \
        -nct 4 \
    	-gt_mode DISCOVERY \
    	-stand_emit_conf 10 \
    	-stand_call_conf 30 \
    	-ERC GVCF \
    	-o $gvcf \
        -log ${bam%.bam}_haplotypecaller.log
```
4. Mark duplicates
5. Call genotypes
```bash
java -Xmx128G -Xms128g -jar $GATK_HOME/GenomeAnalysisTK.jar \
    -T GenotypeGVCFs \
    -R $GENOME \
    -nt 50 \
    --variant combined.g.vcf \
    -o combined_genotypes.vcf.gz
```

```
ln -s ../gatk3/combined_genotypes.vcf.gz gatk3.vcf.gz
```

# Freebayes Calls

1. Use mapped bams from GATK pipeline
2. Run freebayes
```bash
freebayes-parallel aten_regions.txt 46 -f aten_final_0.1.fasta -L lowcov_bams.txt \
	--genotype-qualities -V -X -E 0 --strict-vcf --populations lowcov_populations.txt  > lowcov_wb_pops.vcf
```

Link to the calls

```
ln -s ../freebayes_qc/lowcov_wb_down_pops.vcf freebayes.vcf
```

# ANGSD Calls

1. Use mapped bams from GATK pipeline
2. Run ANGSD in SAMtools calling mode
```bash
angsd -bam lowcov_down.bamlist -GL 1 -out angsd_samtools -doMaf 2 -SNP_pval 1e-6 -doMajorMinor 1
```
3. Run ANGSD in GATK calling mode
```bash
angsd -bam lowcov_down.bamlist -GL 2 -out angsd_gatk -doMaf 2 -SNP_pval 1e-6 -doMajorMinor 1
```

