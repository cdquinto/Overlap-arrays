########Intersect of MEGA arrays and other commercial arrays
#####  /data/users/coquinto/intersect_arrays

###I downloaded the list of SNPs in different arrays from the UCSC

1) affy6.0.txt,  909297 sites (lines in files)
2) illumina650.txt, 660388 sites
3) illumina1M-Duo.txt, 1217520 sites
4) illuminaOmni1-Q.txt, 1169872 sites 

###I obtained the list of SNPs from the Lazaridis paper
## Human Origins 1 Array (Affymetrix)
##in /data/users/coquinto/EuropeFullyPublic: 600841 release.snp

###And I got a bed file for the OmniExpress Illumina (from Ariella)
### annotOmniEx.txt, 726482 sites

####I found a bed file with the SNPs in the Infinium OmniExpress-24 Kit.
InfiniumOmniExpress-24v1-2_A1.bed

###Cut the columns in each file, so I have the chromosome info, physical position and rsID

for file in illumina650.txt illumina1M-Duo.txt illuminaOmni1-Q.txt
do
sed 1d $file > out
cut -f2,4,5 out > ${file}_2
done 

sed 1d affy6.0.txt > out
cut -f2,4,9 out > out2
mv out2 affy6.0.txt_2
rm out 

sed 's/ \+/\t/g' annotOmniEx.txt > out 
sed 1d out > out2
cut -f2,3,6 out2 > annotOmniEx.txt_2

cut -f2 annotOmniEx.txt_2 > out1
cut -f3 annotOmniEx.txt_2 > out2
cut -f1 annotOmniEx.txt_2 > out3
paste out1 out2 out3 > out 
mv out annotOmniEx.txt_2

sed -i 's/ \+/\t/g' release.snp
sed -i 's/^\t//' release.snp 

cut -f2 release.snp > out1
cut -f4 release.snp > out2
cut -f1 release.snp > out3
paste out1 out2 out3 > out
mv out release.snp_2 

sed -i 's/^chr23/chrX/' release.snp_2
sed -i 's/^chr24/chrY/' release.snp_2

sed 1d InfiniumOmniExpress-24v1-2_A1.bed > out 
cut -f1,3,4 out > InfiniumOmniExpress-24v1-2_A1.bed_2

####Sort by physical position

for file in affy6.0.txt_2 annotOmniEx.txt_2 illumina1M-Duo.txt_2 illumina650.txt_2 illuminaOmni1-Q.txt_2 release.snp_2
do
sort -k1,1 -k2,2n ${file} > ${file}_temp
mv ${file}_temp ${file}
done 

###Add chr at the beginning
sed -i 's/^/chr/g' annotOmniEx.txt_2
sed -i 's/^/chr/g' release.snp_2 

sed -i 1d HumanOmni5-4v1_C.bed
cut -f1,3,4 HumanOmni5-4v1_C.bed > HumanOmni5-4v1_C.bed_2

###Final files: affy6.0.txt_2, annotOmniEx.txt_2, illumina1M-Duo.txt_2, illumina650.txt_2, illuminaOmni1-Q.txt_2, release.snp_2, InfiniumOmniExpress-24v1-2_A1.bed_2, entrenamiento_16sample.map_2, HumanOmni5-4v1_C.bed_2

###Files from the MEGA arrays (Consortium, and AMR-AFR)
##Consortium: entrenamiento_clusterGLOBAL.95.snps.rs.map 
##AMR-AFR: maavp1v1.95.snps.rs.map

for file in entrenamiento_clusterGLOBAL.95.snps.rs.map 
do
cut -f1 ${file} > out1
cut -f4 ${file} > out2
cut -f2 ${file} > out3
paste out1 out2 > out
mv out ${file}
done

sed -i 's/^/chr/g' maavp1v1.95.snps.rs.map
sed -i 's/^/chr/g' entrenamiento_clusterGLOBAL.95.snps.rs.map ##(b37)

mv maavp1v1.95.snps.rs.map maavp1v1.95.snps.rs.map_2
mv entrenamiento_clusterGLOBAL.95.snps.rs.map entrenamiento_clusterGLOBAL.95.snps.rs.map_2

cut -f1,4 entrenamiento_16sample.map > entrenamiento_16sample.map_2
sed -i 's/^/chr/' entrenamiento_16sample.map_2

##AMR-AFR: Multi-EthnicAMR-AFR-8v1-0_A1.bed
##Global: Multi-EthnicGlobal_B1.bed

sed -i 1d Multi-EthnicAMR-AFR-8v1-0_A1.bed
cut -f1,3,4 Multi-EthnicAMR-AFR-8v1-0_A1.bed > Multi-EthnicAMR-AFR-8v1-0_A1.bed_2

sed -i 1d Multi-EthnicGlobal_B1.bed
cut -f1,3,4 Multi-EthnicGlobal_B1.bed > Multi-EthnicGlobal_B1.bed_2

#####Results
#####Find SNPs in common (autosomes, X and Y)

perl get_common_snps_doron.pl affy6.0.txt_2 annotOmniEx.txt_2 illumina1M-Duo.txt_2 illumina650.txt_2 illuminaOmni1-Q.txt_2 release.snp_2 maavp1v1.95.snps.rs.map_2 entrenamiento_clusterGLOBAL.95.snps.rs.map_2
##1977 SNPs in common between all arrays

perl get_common_snps_doron.pl affy6.0.txt_2 InfiniumOmniExpress-24v1-2_A1.bed_2 illumina1M-Duo.txt_2 illumina650.txt_2 illuminaOmni1-Q.txt_2 release.snp_2 maavp1v1.95.snps.rs.map_2 entrenamiento_clusterGLOBAL.95.snps.rs.map_2
##1984 SNPs

perl get_common_snps_doron.pl maavp1v1.95.snps.rs.map_2 entrenamiento_clusterGLOBAL.95.snps.rs.map_2
###1330634 SNPs

perl get_common_snps_doron.pl Multi-EthnicAMR-AFR-8v1-0_A1.bed_2 Multi-EthnicGlobal_B1.bed_2
###1412629 SNPs

for file in affy6.0.txt_2 annotOmniEx.txt_2 illumina1M-Duo.txt_2 illumina650.txt_2 illuminaOmni1-Q.txt_2 release.snp_2 InfiniumOmniExpress-24v1-2_A1.bed_2 Multi-EthnicAMR-AFR-8v1-0_A1.bed_2 Multi-EthnicGlobal_B1.bed_2 HumanOmni5-4v1_C.bed_2
do
echo $file 
perl get_common_snps_doron.pl ${file} entrenamiento_clusterGLOBAL.95.snps.rs.map.bed_2
done 

1) affy6.0.txt_2 (b37): 161709
2) annotOmniEx.txt_2 (b37): 301211
3) illumina1M-Duo.txt_2 (b37): 287696
4) illumina650.txt_2 (b37): 206736
5) illuminaOmni1-Q.txt_2 (b37): 328366
6) release.snp_2 (b37): 110563
7) InfiniumOmniExpress-24v1-2_A1.bed_2 (b37): 300725
8) Multi-EthnicAMR-AFR-8v1-0_A1.bed_2 (b37): 1298920
9) Multi-EthnicGlobal_B1.bed_2 (b37): 1333406
10) HumanOmni5-4v1_C.bed_2 (b37): 610522

########

for file in affy6.0.txt_2 annotOmniEx.txt_2 illumina1M-Duo.txt_2 illumina650.txt_2 illuminaOmni1-Q.txt_2 release.snp_2 InfiniumOmniExpress-24v1-2_A1.bed_2 Multi-EthnicAMR-AFR-8v1-0_A1.bed_2 HumanOmni5-4v1_C.bed_2
do 
perl get_common_snps_doron.pl ${file} Multi-EthnicGlobal_B1.bed_2
done 

1) affy6.0.txt_2 (b37): 165795
2) annotOmniEx.txt_2 (b37): 313327
3) illumina1M-Duo.txt_2 (b37): 308002
4) illumina650.txt_2 (b37): 211845
5) illuminaOmni1-Q.txt_2 (b37): 347757
6) release.snp_2 (b37): 105771
7) InfiniumOmniExpress-24v1-2_A1.bed_2 (b37): 311759
8) Multi-EthnicAMR-AFR-8v1-0_A1.bed_2 (b37): 1412629
9) HumanOmni5-4v1_C.bed_2 (b37): 720300

