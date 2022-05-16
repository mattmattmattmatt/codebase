Matt Cook collaborator Chris Kiely looking at samples for IBD/PSC.

Details in PSC_cohort_info.xlsx

Downloaded single VASP reports.
Ran VASP on all 25 samples as well as for IBD gene list and PSC gene list.

Quick look myself yielding the following possible causal variants:

cohort144 -> EPOR missense mutation
cohort142 -> CDK12 missense mutation
cohort141 -> TACC2 missense mutation (if mutation validates via Sanger)
cohort149 -> UGT2B7 missense mutation (homozygous in child and ref in parent)
cohort143 -> GALC two adjacent rare mutations in same codon

Make de novo file
>head -1 APOSLE_cohort141_sg1_3_samples_snv_23-10-2017.tsv > de_novo_PSC.tsv
>grep 'ref_parents het_child' APOSLE_cohort*snv_23-10-2017.tsv | grep PASS >> de_novo_PSC.tsv 

Make other non-mendelian cases
>head -1 APOSLE_cohort141_sg1_3_samples_snv_23-10-2017.tsv > hom_child_ref_parent.tsv 
>grep -v mendelian_rules APOSLE_cohort14*snv*.tsv | grep PASS  | grep -v 'ref_parents het_child' | cut -d':' -f2- | grep -v '	X	' | grep -v '	Y	' >> hom_child_ref_parent.tsv 