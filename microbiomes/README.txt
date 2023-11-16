Mining Fish Microbiome Studies for Ichthycolid Sequences

For each study run this:
```
blastn -task megablast -num_threads 4 -outfmt "6 qseqid qlen length bitscore qcovs evalue pident sseqid" -db Kim_ASVs_DB -query ../ichty_apico_ssu.txt -out blast.txt
```
keep anything above 85% similarity

Extract all matches in fasta form and compile into master fasta:
```
grep -w -A 1 -f  list.txt rep_set.fa | grep -v -- "^--$" > extracted.fa

cat extracted.fa ../recovered_ichthycolids.fa >> ../recovered_ichthycolids.fa 
```

1. Kim, Pil Soo, et al. "Host habitat is the major determinant of the gut microbiome of fish." Microbiome 9.1 (2021): 166.

2. Sylvain, François-Étienne, et al. "Fish skin and gut microbiomes show contrasting signatures of host species and habitat." Applied and environmental microbiology 86.16 (2020): e00789-20.
Found 


Map to Apico SSU Tree
```
muscle -in ../recovered_ichthycolids.fa  -out queries.afa
muscle -profile -in1 Jan_seqs_plusmore.trim.fa -in2 queries.afa -out epa_all.aligned.fa
raxmlHPC -T 2 -m GTRCAT -f v -G 0.2 -n EPARUN -s epa_all.aligned.fa -t RAxML_bipartitions.EPA.tre 
```

More studies
3. Chiarello, M., Auguet, JC., Bettarel, Y. et al. Skin microbiome of coral reef fish is highly variable and driven by host phylogeny and diet. Microbiome 6, 147 (2018). https://doi.org/10.1186/s40168-018-0530-4
```
makeblastdb -in RefSeq_chiarello_et_al_2018.fasta -dbtype nucl -out chiarelloDB
blastn -task megablast -num_threads 4 -outfmt "6 qseqid qlen length bitscore qcovs evalue pident sseqid" -db chiarelloDB -query ../ichty_apico_ssu.txt -out blast.txt
```
keep anything above 85% similarity

Extract all matches in fasta form and compile into master fasta:
```
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' RefSeq_chiarello_et_al_2018.fasta > RefSeq_chiarello_et_al_2018_oneline.fasta
grep -w -A 1 -f  keep.txt RefSeq_chiarello_et_al_2018_oneline.fasta | grep -v -- "^--$" > extracted.fa
cat extracted.fa ../recovered_ichthycolids.fa >> ../recovered_ichthycolids2.fa 
```

OTU01441 found in 15 of 44 screened fish species SKIN. Not water samples. 

4. Christopher R J Kavazos and others, Intestinal Microbiome Richness of Coral Reef Damselfishes (Actinopterygii: Pomacentridae), Integrative Organismal Biology, Volume 4, Issue 1, 2022, obac026, https://doi.org/10.1093/iob/obac026
- Exported DNAseqs from Phyloseq Object
```
makeblastdb -in output.fasta -dbtype nucl -out KavanosDB
blastn -task megablast -num_threads 4 -outfmt "6 qseqid qlen length bitscore qcovs evalue pident sseqid" -db KavanosDB -query ../ichty_apico_ssu.txt -out blast.txt
```
keep anything above 85% similarity

Extract all matches in fasta form and compile into master fasta:
```
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' output.fasta > output_oneline.fasta
grep -w -A 1 -f  keep.txt output_oneline.fasta | grep -v -- "^--$" > extracted.fa
cat extracted.fa ../recovered_ichthycolids2.fa >> ../recovered_ichthycolids3.fa 
```
Lots of blast macthes!


5. Chiarello, Marlène et al. (2020), Data from: Exceptional but vulnerable microbial diversity in coral reef animal surface microbiomes, Dryad, Dataset, https://doi.org/10.5061/dryad.wh70rxwjw
```
makeblastdb -in Chiarello.fasta -dbtype nucl -out ChiarelloDB
blastn -task megablast -num_threads 4 -outfmt "6 qseqid qlen length bitscore qcovs evalue pident sseqid" -db ChiarelloDB -query ../ichty_apico_ssu.txt -out blast.txt
```
keep anything above 85% similarity

Extract all matches in fasta form and compile into master fasta:
```
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' Chiarello.fasta > Chiarello_oneline.fasta
grep -w -A 1 -f  keep.txt Chiarello_oneline.fasta | grep -v -- "^--$" > extracted.fa
cat extracted.fa ../recovered_ichthycolids3.fa >> ../recovered_ichthycolids4.fa 
```


6. Galand, P.E., Ruscheweyh, HJ., Salazar, G. et al. Diversity of the Pacific Ocean coral reef microbiome. Nat Commun 14, 3039 (2023). https://doi.org/10.1038/s41467-023-38500-x

```
makeblastdb -in galand.fa -dbtype nucl -out galandDB
blastn -task megablast -num_threads 4 -outfmt "6 qseqid qlen length bitscore qcovs evalue pident sseqid" -db galandDB -query ../ichty_apico_ssu.txt -out blast.txt
```
keep anything above 85% similarity

Extract all matches in fasta form and compile into master fasta:
```
grep -w -A 1 -f  keep.txt galand.fa | grep -v -- "^--$" > extracted.fa
cat extracted.fa ../recovered_ichthycolids4.fa >> ../recovered_ichthycolids5.fa 
```

7. León-Zayas R, McCargar M, Drew JA, Biddle JF. 2020. Microbiomes of fish, sediment and seagrass suggest connectivity of coral reef microbial populations. PeerJ 8:e10026 https://doi.org/10.7717/peerj.10026
- Raw SRA data is being searched here
```
makeblastdb -in ERR1019614.fasta -dbtype nucl -out ZayasDB
blastn -task megablast -num_threads 4 -outfmt "6 qseqid qlen length bitscore qcovs evalue pident sseqid" -db ZayasDB -query ../ichty_apico_ssu.txt -out blast.txt
```
keep anything above 85% similarity

Extract all matches in fasta form and compile into master fasta:
```
grep -w -A 1 -f  keep.txt ERR1019614.fasta | grep -v -- "^--$" > extracted.fa
cat extracted.fa ../recovered_ichthycolids5.fa >> ../recovered_ichthycolids6.fa 
```

Going to map these 1903 seqs to get preliminary look.

Map to Apico SSU Tree
```
muscle -in ../recovered_ichthycolids6.fa  -out queries.afa
muscle -profile -in1 Jan_seqs_plusmore.trim.fa -in2 queries.afa -out epa_all.aligned.fa
seqkit rmdup -s < epa_all.aligned.fa > epa_all.aligned_nodups.fa
awk '/^>/ { f = !a[$0]++ } f' epa_all.aligned.fa > epa_all.aligned_nodups.fa
raxmlHPC -T 2 -m GTRCAT -f v -G 0.2 -n EPARUN -s epa_all.aligned_nodups.fa -t RAxML_bipartitions.9E_api_plastid_ssu_tree11
sed 's/QUERY___//g' RAxML_labelledTree.EPARUN | sed 's/\[I[0-9]*\]//g' > RAxML_placement_tree_round2.tre
```

Tree is bad... need to re-run blast on final set then try again...


Toadfish
```
makeblastdb -in ASVs.fa -dbtype nucl -out toadfishDB
blastn -task megablast -num_threads 4 -outfmt "6 qseqid qlen length bitscore qcovs evalue pident sseqid" -db toadfishDB -query ../ichty_apico_ssu.txt -out blast.txt
```
keep anything above 85% similarity

Extract all matches in fasta form and compile into master fasta:
```
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' ASVs.fa > ASVs_oneline.fasta
grep -w -A 1 -f  keep.txt ASVs_oneline.fasta | grep -v -- "^--$" > extracted.fa
```
No matches in blast, won't add. 

8. Minich, J.J., Härer, A., Vechinski, J. et al. Host biology, ecology and the environment influence microbial biomass and diversity in 101 marine fish species. Nat Commun 13, 6978 (2022). https://doi.org/10.1038/s41467-022-34557-2
```
makeblastdb -in Minich_ASVs.fa -dbtype nucl -out MinichDB
blastn -task megablast -num_threads 4 -outfmt "6 qseqid qlen length bitscore qcovs evalue pident sseqid" -db MinichDB -query ../ichty_apico_ssu.txt -out blast.txt
```
keep anything above 85% similarity

Extract all matches in fasta form and compile into master fasta:
```
grep -w -A 1 -f  keep.txt Minich_ASVs.fa| grep -v -- "^--$" > extracted.fa
cat extracted.fa ../recovered_ichthycolids6.fa >> ../recovered_ichthycolids7.fa 
```



## Now Re-run Blast on this new recovered database and take anything above 95%?
```
makeblastdb -in recovered_ichthycolids7.fa  -dbtype nucl -out recoveredDB
blastn -task megablast -num_threads 4 -outfmt "6 qseqid qlen length bitscore qcovs evalue pident sseqid" -db recoveredDB -query ichty_apico_ssu.txt -out blast.txt -max_target_seqs 10000
grep -w -A 1 -f  95.txt recovered_ichthycolids7.fa | grep -v -- "^--$" > 95.fa
```
Now map to tree again:
```
muscle -in ../../95.fa -out queries.afa
muscle -profile -in1 Jan_seqs_plusmore.trim.fa -in2 queries.afa -out epa_all.aligned.fa
seqkit rmdup -s < epa_all.aligned.fa > epa_all.aligned_nodups.fa
awk '/^>/ { f = !a[$0]++ } f' epa_all.aligned.fa > epa_all.aligned_nodups.fa
raxmlHPC -T 2 -m GTRCAT -f v -G 0.2 -n EPARUN -s epa_all.aligned_nodups.fa -t RAxML_bipartitions.9E_api_plastid_ssu_tree11
sed 's/QUERY___//g' RAxML_labelledTree.EPARUN | sed 's/\[I[0-9]*\]//g' > RAxML_placement_tree_round2.tre
```

Grabbing ones that branched near Ichthy for full blast search...
```
grep -w -A 1 -f  Probable_Ichthys.txt 95.fa | grep -v -- "^--$" > probable_ichthys.fa
```



# Maybe incorporate some metagenomics samples
1. Fergus W. J. Collins, Calum J. Walsh, Beatriz Gomez-Sala, Elena Guijarro-García, David Stokes, Klara B. Jakobsdóttir, Kristján Kristjánsson, Finlay Burns, Paul D. Cotter, Mary C. Rea, Colin Hill & R. Paul Ross (2021) The microbiome of deep-sea fish reveals new microbial species and a sparsity of antibiotic resistance genes, Gut Microbes, 13:1, DOI: 10.1080/19490976.2021.1921924
testing on one metagenome first
- Pachystomias microdon (ERR6050740)
```
makeblastdb -in ERR6050740.fasta -dbtype nucl -out ERR6050740DB
blastn -task megablast -num_threads 4 -outfmt "6 qseqid qlen length bitscore qcovs evalue pident sseqid" -db ERR6050740DB -query ../ichty_apico_ssu.txt -out blast.txt
```
keep anything above 85% similarity
Extract all matches in fasta form and compile into master fasta:
```
grep -w -A 1 -f  keep.txt ERR6050740.fasta | grep -v -- "^--$" > extracted.fa
cat extracted.fa ../recovered_ichthycolids7.fa >> ../recovered_ichthycolids8.fa 
```

Actually no, lets do MGnify fish metagenomes: https://www.ebi.ac.uk/metagenomics/browse/studies?biome=root%3AHost-associated%3AFish&page_size=50 
1. Iwatsuki, Toshihide, et al. "16S rRNA gene amplicon sequencing of gut microbiota in three species of deep-sea fish in Suruga Bay, Japan." Microbiology resource announcements 10.1 (2021): e01260-20. (https://www.ebi.ac.uk/metagenomics/studies/MGYS00005732#overview)
```
makeblastdb -in iwatsuki.fa -dbtype nucl -out iwatsukiDB
blastn -task megablast -num_threads 4 -outfmt "6 qseqid qlen length bitscore qcovs evalue pident sseqid" -db iwatsukiDB -query ../ichty_apico_ssu.txt -out blast.txt
```
keep anything above 85% similarity
Extract all matches in fasta form and compile into master fasta:
```
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' iwatsuki.fa > iwatsuki_ol.fa
grep -w -A 1 -f  keep.txt iwatsuki_ol.fa | grep -v -- "^--$" > extracted.fa
cat extracted.fa  >> ../recovered_ichthycolids_magnify.fa 
```
Looks great!!!!

2. Thépot V, Slinger J, Rimmer MA, Paul NA, Campbell AH. Is the Intestinal Bacterial Community in the Australian Rabbitfish Siganus fuscescens Influenced by Seaweed Supplementation or Geography? Microorganisms. 2022 Feb 23;10(3):497. doi: 10.3390/microorganisms10030497. PMID: 35336073; PMCID: PMC8954549. https://www.ebi.ac.uk/metagenomics/studies/MGYS00003166#overview
Grabbed small sampling of these...
```
makeblastdb -in thepot.fa -dbtype nucl -out thepotDB
blastn -task megablast -num_threads 4 -outfmt "6 qseqid qlen length bitscore qcovs evalue pident sseqid" -db thepotDB -query ../ichty_apico_ssu.txt -out blast.txt
```
Nothing there above 85%

3. Salmon
```
makeblastdb -in salmon.fa -dbtype nucl -out salmonDB
blastn -task megablast -num_threads 4 -outfmt "6 qseqid qlen length bitscore qcovs evalue pident sseqid" -db salmonDB -query ../ichty_apico_ssu.txt -out blast.txt
```
Nothing there above 85%

4. Metagenome of grass carp intestinal contents and mucosa (https://www.ebi.ac.uk/metagenomics/studies/MGYS00000380#overview)

```
makeblastdb -in carp.fa -dbtype nucl -out carpDB
blastn -task megablast -num_threads 4 -outfmt "6 qseqid qlen length bitscore qcovs evalue pident sseqid" -db carpDB -query ../ichty_apico_ssu.txt -out blast.txt
```
Not good either, Proceed with tree with deep sea fish bc that's very interesting


```
muscle -in recovered_ichthycolids_magnify.fa -out queries.afa
muscle -profile -in1 Jan_seqs_plusmore.trim.fa -in2 queries.afa -out epa_all.aligned.fa
awk '/^>/ { f = !a[$0]++ } f' epa_all.aligned.fa > epa_all.aligned_nodups.fa
raxmlHPC -T 2 -m GTRCAT -f v -G 0.2 -n EPARUN -s epa_all.aligned_nodups.fa -t RAxML_bipartitions.9E_api_plastid_ssu_tree11
sed 's/QUERY___//g' RAxML_labelledTree.EPARUN | sed 's/\[I[0-9]*\]//g' > RAxML_placement_tree.tre
```

Grabbing ones that branched near Ichthy for full blast search...
```
grep -w -A 1 -f  probable_ichtys.txt recovered_ichthycolids_magnify.fa | grep -v -- "^--$" > probable_ichthys_magnify.fa
```


# Another batch from Jake Minich
```
makeblastdb -in ASVs.fa -dbtype nucl -out ASVsDB
blastn -task megablast -num_threads 4 -outfmt "6 qseqid qlen length bitscore qcovs evalue pident sseqid" -db ASVsDB -query ../../ichty_apico_ssu.txt -out blast.txt

grep -w -A 1 -f  keep.txt ASVs.fa | grep -v -- "^--$" > extracted.fa
```
Gonna re-run the tree fresh from here. 
```
muscle -in ../extracted.fa -out queries.afa
muscle -profile -in1 Jan_seqs_plusmore.trim.fa -in2 queries.afa -out epa_all.aligned.fa
awk '/^>/ { f = !a[$0]++ } f' epa_all.aligned.fa > epa_all.aligned_nodups.fa
raxmlHPC -T 2 -m GTRCAT -f v -G 0.2 -n EPARUN -s epa_all.aligned_nodups.fa -t RAxML_bipartitions.9E_api_plastid_ssu_tree12_GTRCAT
sed 's/QUERY___//g' RAxML_labelledTree.EPARUN | sed 's/\[I[0-9]*\]//g' > RAxML_placement_tree.tre

grep -w -A 1 -f  ichthys.txt extracted.fa | grep -v -- "^--$" > Minich_ichthys.fa
```


HMM WAY
```
perl ~/Desktop/scripts/fasta2stockholm.pl Jan_seqs_plusmore.aligned.fa > Jan_seqs_plusmore.aligned.sto
hmmbuild Jan_seqs_plusmore.hmm Jan_seqs_plusmore.aligned.sto
hmmalign --dna --mapali Jan_seqs_plusmore.aligned.fa Jan_seqs_plusmore.hmm all_recovered_ichtys.fa > epa_all.aligned.sto
perl ~/Desktop/scripts/stockholm2fasta.pl epa_all.aligned.sto > epa_all.aligned.fa
awk '/^>/ {print; next} {gsub("\\.", "-"); print}' epa_all.aligned.fa > cleaned_epa_all.aligned.fa
raxmlHPC -T 4 -m GTRCAT -f v -G 0.2 -n EPARUN -s cleaned_epa_all.aligned.fa -t RAxML_bipartitions.9E_api_plastid_ssu_tree13.tre 
sed 's/QUERY___//g' RAxML_labelledTree.EPARUN | sed 's/\[I[0-9]*\]//g' > RAxML_placement_tree.tre
```


Lungfish skin explant sample from https://www.science.org/doi/full/10.1126/sciadv.abj0829 
```
vsearch -sortbylength SRR12533720.fasta -output SRR12533720.sorted.fasta -minseqlength 150
makeblastdb -in SRR12533720.sorted.fasta -dbtype nucl -out LungfishDB
blastn -task megablast -num_threads 4 -outfmt "6 qseqid qlen length bitscore qcovs evalue pident sseqid" -db LungfishDB -query ../ichty_apico_ssu.txt -out blast.txt
```
keep anything above 90% similarity
Extract all matches in fasta form and compile into master fasta:
```
grep -w -A 1 -f  keep.txt SRR12533720.fasta | grep -v -- "^--$" > extracted.fa
```

Make Tree
```
hmmalign --dna --mapali Jan_seqs_plusmore.aligned.fa Jan_seqs_plusmore.hmm extracted.fa > epa_all.aligned.sto
perl ~/Desktop/scripts/stockholm2fasta.pl epa_all.aligned.sto > epa_all.aligned.fa
awk '/^>/ {print; next} {gsub("\\.", "-"); print}' epa_all.aligned.fa > cleaned_epa_all.aligned.fa
raxmlHPC -T 4 -m GTRCAT -f v -G 0.2 -n EPARUN -s cleaned_epa_all.aligned.fa -t RAxML_bipartitions.9E_api_plastid_ssu_tree13.tre 
sed 's/QUERY___//g' RAxML_labelledTree.EPARUN | sed 's/\[I[0-9]*\]//g' > RAxML_placement_tree.tre
```



Now shark study
https://animalmicrobiome.biomedcentral.com/articles/10.1186/s42523-022-00168-x#availability-of-data-and-materials
```
vsearch -sortbylength shark.fa -output shark.sorted.fasta -minseqlength 150
makeblastdb -in shark.sorted.fasta -dbtype nucl -out sharkDB
blastn -task megablast -num_threads 4 -outfmt "6 qseqid qlen length bitscore qcovs evalue pident sseqid" -db sharkDB -query ../ichty_apico_ssu.txt -out blast.txt
```
Nothing above 82%

Seminal and blood plasma microbiome in the in Scyliorhinus canicula
```
vsearch -sortbylength Scan.fa -output Scan.sorted.fasta -minseqlength 150
makeblastdb -in Scan.sorted.fasta -dbtype nucl -out ScanDB
blastn -task megablast -num_threads 4 -outfmt "6 qseqid qlen length bitscore qcovs evalue pident sseqid" -db ScanDB -query ../ichty_apico_ssu.txt -out blast.txt
```








