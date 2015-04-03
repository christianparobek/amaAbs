## A bash script to do the PLINKSEQ analysis
## Hopefully PLINKSEQ will actually do what we need it to
## Started 2 April 2015
## Christian P

projectdir=/proj/julianog/users/ChristianP/amaAbs
vcfdir=/proj/julianog/users/ChristianP/amaAbs/variants/split_indivs_HC

# Declare a new project
pseq amaAbs-assoc new-project

# Load the VCFs into the new project
pseq amaAbs-assoc load-vcf --vcf $vcfdir/NB3_HC.sans0.vcf $vcfdir/NB4_HC.sans0.vcf \
	$vcfdir/OA1_HC.sans0.vcf $vcfdir/OB1_HC.sans0.vcf \
	$vcfdir/PB4_HC.sans0.vcf $vcfdir/PBSN_HC.sans0.vcf \
	$vcfdir/PBSO_HC.sans0.vcf $vcfdir/PC1_HC.sans0.vcf \
	$vcfdir/QVN_HC.sans0.vcf $vcfdir/QVO_HC.sans0.vcf


pseq amaAbs-assoc load-vcf --vcf $vcfdir/../all10_HC.pass.vcf


# Add the pheno info
pseq amaAbs-assoc load-pheno --file selected.phe 
