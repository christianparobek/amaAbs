#projectdir=/proj/julianog/users/ChristianP/amaAbs
#vcfdir=/proj/julianog/users/ChristianP/amaAbs/variants/split_indivs_HC



ref=/proj/julianog/refs/Pf3D7_v13.0/PlasmoDB-13.0_Pfalciparum3D7_Genome.fasta
gatk=/nas02/apps/biojars-1.0/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar
vcfdir=/proj/julianog/users/ChristianP/amaAbs/variants/split_indivs_HC

###########################################
###### GET THE SELECTED INTERSECTION ######
###########################################

java -jar $gatk -T CombineVariants \
	-R $ref \
	-V:nb3 $vcfdir/NB3_HC.sans0.vcf \
	-V:nb4 $vcfdir/NB4_HC.sans0.vcf \
	-V:oa1 $vcfdir/OA1_HC.sans0.vcf \
	-V:ob1 $vcfdir/OB1_HC.sans0.vcf \
	-V:qvn $vcfdir/QVN_HC.sans0.vcf \
	-V:qvo $vcfdir/QVO_HC.sans0.vcf \
	-o sel_union.vcf

java -jar $gatk -T SelectVariants \
	-R $ref \
	-V:variant sel_union.vcf \
	-select 'set == "Intersection";' \
	-o sel_intersect.vcf

###########################################
##### GET THE UNSELECTED INTERSECTION #####
###########################################

java -jar $gatk -T CombineVariants \
	-R $ref \
	-V:pb4 $vcfdir/PB4_HC.sans0.vcf \
	-V:pc1 $vcfdir/PC1_HC.sans0.vcf \
	-V:pbsn $vcfdir/PBSN_HC.sans0.vcf \
	-V:pbso $vcfdir/PBSO_HC.sans0.vcf \
	-o unsel_union.vcf

java -jar $gatk -T SelectVariants \
	-R $ref \
	-V:variant unsel_union.vcf \
	-select 'set == "Intersection";' \
	-o unsel_intersect.vcf


###########################################
##### GET VARIANTS THAT AREN'T SHRAED #####
###########################################

bedtools intersect -v -a sel_intersect.vcf \
	-b unsel_intersect.vcf \
	> in_all_SEL_not_in_all_UNSEL.vcf
		# -v reports only entries in A w no overlaps in B

bedtools intersect -header -v -a sel_intersect.vcf \
	-b unsel_union.vcf \
	> in_all_SEL_in_no_UNSEL.vcf
		# -v reports only entries in A w no overlaps in B

bedtools intersect -v -a unsel_intersect.vcf \
	-b sel_union.vcf \
	> in_all_UNSEL_in_no_SEL.vcf
		# -v reports only entries in A w no overlaps in B

###########################################
######### WHICH VARIANTS IN EXONS #########
###########################################

## Make the necessary exon and non-exon GFF definition files
grep "exon\|##" PlasmoDB-13.0_Pfalciparum3D7.gff | grep -v "gene" | head -n -1 > PlasmoDB-13.0_Pfalciparum3D7_exons.gff
    # Make an exon-only GFF file (includes header)
    # Grep out any gene entries that slipped in because they had "exon" in their description
    # Remove the last line of the file, since it starts with '##' 

## Use bedtools to subset exonic and non-exonic regions 
bedtools intersect -a in_all_SEL_in_no_UNSEL.vcf -b PlasmoDB-13.0_Pfalciparum3D7_exons.gff > in_all_SEL_in_no_UNSEL_exons.vcf


###########################################
############ ANNOTATE VARIANTS ############
###########################################

## RUN snpEff
java -Xmx4g -jar ~/snpEff/snpEff.jar \
	-v \
	-o gatk \
	Pf3D7v90 \
	in_all_SEL_in_no_UNSEL.vcf \
	> in_all_SEL_in_no_UNSEL.anno.vcf



