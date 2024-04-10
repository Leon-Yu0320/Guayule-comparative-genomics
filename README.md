# **Guayule genome project**


**Author: Li'ang Yu, Andrew Nelson** \
**Date: Oct 31st, 2023**

<img width="1011" alt="Screenshot 2024-04-10 at 3 48 22 PM" src="https://github.com/Leon-Yu0320/Guayule-comparative-genomics/assets/69836931/36e7a99d-27dc-4799-852d-38e1316023ef">

Guayule (Parthenium argentatum), a woody desert shrub indigenous to the Chihuahuan desert of North America, is under development in the semi-arid southwestern U.S., driven by a need for high quality natural rubber, resin, and bioenergy feedstock produced by the plant. 

Here we performed the comparative genomics analysis of two genomes of guayule,  a diploid guayule **(2X)**, CAL-3, which is a disease-resistant selection from W6-429 9, as well as AZ-2, a allo-tetraploid guayule hybrid **(4X)**, and AZ-6, a pure guayule tetraploid **(2X)**. Both AZ-2 and AZ-6 reproduce apomictically. 

While a draft genome assembly is available for W6-429/CAL-3 1, most breeding efforts have been placed in the agronomically important accessions AZ-2 and AZ-6. AZ-2 has been the focus of most functional genomics as it is transformable, whereas AZ-6 contains significantly more rubber at harvest than AZ-2. In addition to producing more rubber, AZ-6 is also morphologically distinct from AZ-2, in that it is a smaller, more compact plant at maturity, allowing for denser planting and possibly more cold tolerance

<img width="369" alt="Screenshot 2024-04-10 at 3 47 13 PM" src="https://github.com/Leon-Yu0320/Guayule-comparative-genomics/assets/69836931/da737278-c281-467a-be3f-ece6f74a576e">


- [**Guayule genome project**](#guayule-genome-project)
  - [Genomic analysis of the genome synteny](#genomic-analysis-of-the-genome-synteny)
    - [Conduct the synteny analyiss by nucmer using genomic sequences](#conduct-the-synteny-analyiss-by-nucmer-using-genomic-sequences)
    - [Using nucmer and mummer for plotting](#using-nucmer-and-mummer-for-plotting)
  - [Test the mapping depth of resequencing data derived from 10 groups of crosses (Diploid \* Tetraploid)](#test-the-mapping-depth-of-resequencing-data-derived-from-10-groups-of-crosses-diploid--tetraploid)
    - [Trim reads](#trim-reads)
    - [Mapping reads to reference genome](#mapping-reads-to-reference-genome)
    - [Remove multiple-mapped reads](#remove-multiple-mapped-reads)
    - [Calculte the coverage in sliding windows](#calculte-the-coverage-in-sliding-windows)
      - [Manually prepare the bed file for each chromosomes (see example)](#manually-prepare-the-bed-file-for-each-chromosomes-see-example)

## Genomic analysis of the genome synteny

**AZ2 genome basci information**
```
Chr01A	98981724	8	80	81
Chr01B	92884050	100219012	80	81
Chr01C	92689182	194264121	80	81
Chr02A	100908417	288111926	80	81
Chr02B	90500992	390281707	80	81
Chr02C	89934939	481913970	80	81
Chr02D	88028799	572973104	80	81
Chr03A	92058258	662102271	80	81
Chr03B	88360942	755311266	80	81
Chr03C	88235095	844776728	80	81
Chr03D	61413590	934114770	80	81
Chr03E	60419288	996296038	80	81
Chr04A	94338614	1057470576	80	81
Chr04B	85535101	1152988431	80	81
Chr04C	84514908	1239592729	80	81
Chr04D	82754871	1325164082	80	81
Chr05A	84685864	1408953397	80	81
Chr05B	83425079	1494697843	80	81
Chr05C	81788810	1579165744	80	81
Chr05D	77648727	1661976923	80	81
Chr06A	85280400	1740596268	80	81
Chr06B	75467198	1826942681	80	81
Chr06C	74680594	1903353227	80	81
Chr06D	71923996	1978967337	80	81
Chr07A	77088775	2051790391	80	81
Chr07B	75821094	2129842784	80	81
Chr07C	73542210	2206611650	80	81
Chr07D	73314558	2281073146	80	81
Chr08A	87508957	2355304144	80	81
Chr08B	75465782	2443906971	80	81
Chr08C	74973946	2520316084	80	81
Chr08D	71924250	2596227213	80	81
Chr09A	72402239	2669050525	80	81
Chr09B	71919192	2742357800	80	81
Chr09C	60905134	2815175990	80	81
Chr10A	84808527	2876842447	80	81
Chr10B	75012979	2962711089	80	81
Chr10C	74595959	3038661739	80	81
Chr11A	76421451	3114190156	80	81
Chr11B	72489106	3191566884	80	81
Chr11C	71502228	3264962112	80	81
Chr12A	85739862	3337358126	80	81
Chr12B	70918791	3424169745	80	81
Chr12C	70266948	3495975029	80	81
Chr12D	69963381	3567120322	80	81
Chr13A	69995893	3637958254	80	81
Chr13B	69325164	3708829104	80	81
Chr13C	66260975	3779020841	80	81
Chr13D	65544011	3846110087	80	81
Chr14A	68925583	3912473407	80	81
Chr14B	67378108	3982260568	80	81
Chr14C	66985537	4050480911	80	81
Chr15A	73824657	4118303776	80	81
Chr15B	65796143	4193051250	80	81
Chr15C	65411185	4259669853	80	81
Chr15D	64528613	4325898686	80	81
Chr16A	68507878	4391233915	80	81
Chr16B	67424057	4460598150	80	81
Chr16C	65161587	4528865016	80	81
Chr16D	53647601	4594841131	80	81
Chr17A	67462663	4649159336	80	81
Chr17B	66488367	4717465291	80	81
Chr17C	65711542	4784784771	80	81
Chr17D	63374426	4851317716	80	81
Chr18A	65377935	4915484331	80	81
Chr18B	55739520	4981679499	80	81
Chr18C	53604530	5038115771	80	81
Chr18D	53290559	5092390366	80	81
```

**Cal3 Hap1-Basic information**
```
Chr01   95457745        7       80      81
Chr02   91772066        96650981        80      81
Chr03   88676301        189570205       80      81
Chr04   87773215        279354967       80      81
Chr05   83151666        368225355       80      81
Chr06   79118053        452416424       80      81
Chr07   77941552        532523460       80      81
Chr08   73639970        611439289       80      81
Chr09   73447079        685999766       80      81
Chr10   75296575        760364941       80      81
Chr11   72769860        836602731       80      81
Chr12   71501908        910282222       80      81
Chr13   67988348        982677911       80      81
Chr14   67308586        1051516121      80      81
Chr15   65991018        1119666072      80      81
Chr16   65851447        1186481985      80      81
Chr17   62634404        1253156583      80      81
Chr18   57812506        1316573925      80      81
scaffold_146    51856   1375109102      80      81
scaffold_154    50231   1375161621      80      81
```

**Cal3 Hap2-Basic information**
```
Chr01   94250353        7       80      81
Chr02   90979286        95428497        80      81
Chr03   87804457        187545032       80      81
Chr04   86616899        276447052       80      81
Chr05   81626467        364146670       80      81
Chr06   78374608        446793475       80      81
Chr07   75387401        526147773       80      81
Chr08   73406134        602477524       80      81
Chr09   73293307        676801242       80      81
Chr10   73788245        751010723       80      81
Chr11   72462808        825721329       80      81
Chr12   70934005        899089930       80      81
Chr13   65977075        970910618       80      81
Chr14   65619891        1037712414      80      81
Chr15   65224688        1104152561      80      81
Chr16   65093094        1170192565      80      81
Chr17   60262870        1236099330      80      81
Chr18   57574277        1297115493      80      81
scaffold_61     109269  1355409462      80      81
```




### Conduct the synteny analyiss by nucmer using genomic sequences
**Extract haplotype sequences for each chromosome group for comparisons**

**For chromosomes with 4 haplotypes**
```bash

genome_dir="/home/lyu/01_project/04_Guayule/01_Genome"
anno_dir="/home/lyu/01_project/04_Guayule/02_Annotation"
work_dir="/home/lyu/01_project/04_Guayule/03_Nucmer/02_Group"

for i in 02 04 05 06 07 08 12 13 15 16 17 18;
do
  ### For AZ-2 genome
  samtools faidx $genome_dir/Parthenium_argentatum_var_AZ2_D.mainGenome.fasta Chr${i}A > $work_dir/AZ2-group${i}.fasta
  samtools faidx $genome_dir/Parthenium_argentatum_var_AZ2_D.mainGenome.fasta Chr${i}B >> $work_dir/AZ2-group${i}.fasta
  samtools faidx $genome_dir/Parthenium_argentatum_var_AZ2_D.mainGenome.fasta Chr${i}C >> $work_dir/AZ2-group${i}.fasta
  samtools faidx $genome_dir/Parthenium_argentatum_var_AZ2_D.mainGenome.fasta Chr${i}D >> $work_dir/AZ2-group${i}.fasta

  ### For Cal-3 genome
  samtools faidx $genome_dir/Parthenium_argentatum_var_Cal3_B2.HAP1.mainGenome.fasta Chr${i} > $work_dir/Cal3-group${i}-Hap1.fasta
  sed -i 's/Chr/Chr_Hap1_/g' $work_dir/Cal3-group${i}-Hap1.fasta

  samtools faidx $genome_dir/Parthenium_argentatum_var_Cal3_B2.HAP2.mainGenome.fasta Chr${i} > $work_dir/Cal3-group${i}-Hap2.fasta
  sed -i 's/Chr/Chr_Hap2_/g' $work_dir/Cal3-group${i}-Hap2.fasta

  cat $work_dir/Cal3-group${i}-Hap1.fasta $work_dir/Cal3-group${i}-Hap2.fasta > $work_dir/Cal3-group${i}-combine.fasta

  ### Perform the number analysis for each chromosome groups
  nucmer -p AZ-Cal3-Group${i} $work_dir/Cal3-group${i}-combine.fasta $work_dir/AZ2-group${i}.fasta
  show-coords AZ-Cal3-Group${i}.delta AZ-Cal3-Group${i}.coords

 ### Perform the alignment analysis for each chromosome groups using mummer
  mummer -mum -b -c -l 5000 $work_dir/Cal3-group${i}-combine.fasta $work_dir/AZ2-group${i}.fasta > $work_dir/Guayule-AZ2-Cal3-Group${i}.mums
  
  ### Plot all Hap-1 relative to other homologs 
  mummerplot -postscript -q Chr${i}A -r Chr_Hap1_${i} -p Guayule-AZ2-Cal3-Hap1-Group${i}-1 $work_dir/01_Guayule-AZ2-Cal3-Group${i}.mums
  mummerplot -postscript -q Chr${i}B -r Chr_Hap1_${i} -p Guayule-AZ2-Cal3-Hap1-Group${i}-2  $work_dir/01_Guayule-AZ2-Cal3-Group${i}.mums
  mummerplot -postscript -q Chr${i}C -r Chr_Hap1_${i} -p Guayule-AZ2-Cal3-Hap1-Group${i}-3  $work_dir/01_Guayule-AZ2-Cal3-Group${i}.mums
  mummerplot -postscript -q Chr${i}D -r Chr_Hap1_${i} -p Guayule-AZ2-Cal3-Hap1-Group${i}-4  $work_dir/01_Guayule-AZ2-Cal3-Group${i}.mums

  ### Plot all Hap-2 relative to other homologs
  mummerplot -postscript -q Chr${i}A -r Chr_Hap2_${i} -p Guayule-AZ2-Cal3-Hap2-Group${i}-1 $work_dir/01_Guayule-AZ2-Cal3-Group${i}.mums
  mummerplot -postscript -q Chr${i}B -r Chr_Hap2_${i} -p Guayule-AZ2-Cal3-Hap2-Group${i}-2  $work_dir/01_Guayule-AZ2-Cal3-Group${i}.mums
  mummerplot -postscript -q Chr${i}C -r Chr_Hap2_${i} -p Guayule-AZ2-Cal3-Hap2-Group${i}-3  $work_dir/01_Guayule-AZ2-Cal3-Group${i}.mums
  mummerplot -postscript -q Chr${i}D -r Chr_Hap2_${i} -p Guayule-AZ2-Cal3-Hap2-Group${i}-4  $work_dir/01_Guayule-AZ2-Cal3-Group${i}.mums
done

```
**For chromosomes with 3 haplotypes**
```bash
genome_dir="/home/lyu/01_project/04_Guayule/01_Genome"
anno_dir="/home/lyu/01_project/04_Guayule/02_Annotation"
work_dir="/home/lyu/01_project/04_Guayule/03_Nucmer/01_trail"

for i in 09,10,11,14;
do
  ### For AZ-2 genome
  samtools faidx $genome_dir/Parthenium_argentatum_var_AZ2_D.mainGenome.fasta Chr${i}A > $work_dir/AZ2-group${i}.fasta
  samtools faidx $genome_dir/Parthenium_argentatum_var_AZ2_D.mainGenome.fasta Chr${i}B >> $work_dir/AZ2-group${i}.fasta
  samtools faidx $genome_dir/Parthenium_argentatum_var_AZ2_D.mainGenome.fasta Chr${i}C >> $work_dir/AZ2-group${i}.fasta

  ### For Cal-3 genome
  samtools faidx $genome_dir/Parthenium_argentatum_var_Cal3_B2.HAP1.mainGenome.fasta Chr${i} > $work_dir/Cal3-group${i}.fasta
  sed -i 's/Chr01/Chr01-Hap1/g' $work_dir/Cal3-group1.fasta
done
```

### Using nucmer and mummer for plotting
```bash
for i in $(seq -w 01 18);
do

  ### build the lib for comparison
  nucmer -p Group1 AZ2-group1.fasta Cal3-group1.fasta
  ### Show coods inforamtion
  show-coords Group1.delta Group1.coords

  ### Show plots
  mummer -mum -b -c Cal3-group1.fasta AZ2-group1.fasta > Guayule.mums

  mummerplot -postscript -p Group01-mummer Guayule.mums
done 

```

Run for the entire combined genome
```bash
genome_dir="/home/lyu/01_project/04_Guayule/01_Genome"
work_dir="/home/lyu/01_project/04_Guayule/03_Nucmer"

nucmer -p Hap1 $genome_dir/Parthenium_argentatum_var_AZ2_D.mainGenome.fasta $genome_dir/Parthenium_argentatum_var_Cal3_B2.HAP1.mainGenome.fasta
nucmer -p Hap2 $genome_dir/Parthenium_argentatum_var_AZ2_D.mainGenome.fasta $genome_dir/Parthenium_argentatum_var_Cal3_B2.HAP2.mainGenome.fasta

show-coords Hap1.delta Hap1.coords
show-coords Hap2.delta Hap2.coords

mummer -mum -b -c -l 2000 $genome_dir/Parthenium_argentatum_var_AZ2_D.mainGenome.fasta $genome_dir/Parthenium_argentatum_var_Cal3_B2.HAP1.mainGenome.fasta > Guayule-Hap1.mums
mummer -mum -b -c -l 2000 $genome_dir/Parthenium_argentatum_var_AZ2_D.mainGenome.fasta $genome_dir/Parthenium_argentatum_var_Cal3_B2.HAP2.mainGenome.fasta > Guayule-Hap2.mums

mummerplot -postscript --color -p Guayule-Hap1-dotplot Guayule-Hap1.mums
mummerplot -postscript --color -p Guayule-Hap2-dotplot Guayule-Hap2.mums
```
Plot the seperate maps for each Hap to homologs groups

```bash 

for i in 02 04 05 06 07 08 12 13 15 16 17 18;
do
  mummerplot -postscript -q Chr${i}A -r Chr_Hap1_${i} -p Guayule-AZ2-Cal3-Hap1-Group${i}-1 Guayule-Hap1.mums
  mummerplot -postscript -q Chr${i}B -r Chr_Hap1_${i} -p Guayule-AZ2-Cal3-Hap1-Group${i}-2 Guayule-Hap1.mums
  mummerplot -postscript -q Chr${i}C -r Chr_Hap1_${i} -p Guayule-AZ2-Cal3-Hap1-Group${i}-3 Guayule-Hap1.mums
  mummerplot -postscript -q Chr${i}D -r Chr_Hap1_${i} -p Guayule-AZ2-Cal3-Hap1-Group${i}-4 Guayule-Hap1.mums
done

for i in 02 04 05 06 07 08 12 13 15 16 17 18;
do
  mummerplot -postscript -q Chr${i}A -r Chr_Hap2_${i} -p Guayule-AZ2-Cal3-Hap2-Group${i}-1 Guayule-Hap2.mums
  mummerplot -postscript -q Chr${i}B -r Chr_Hap2_${i} -p Guayule-AZ2-Cal3-Hap2-Group${i}-2 Guayule-Hap2.mums
  mummerplot -postscript -q Chr${i}C -r Chr_Hap2_${i} -p Guayule-AZ2-Cal3-Hap2-Group${i}-3 Guayule-Hap2.mums
  mummerplot -postscript -q Chr${i}D -r Chr_Hap2_${i} -p Guayule-AZ2-Cal3-Hap2-Group${i}-4 Guayule-Hap2.mums
done
```

## Test the mapping depth of resequencing data derived from 10 groups of crosses (Diploid * Tetraploid)

### Trim reads
```bash
work_dir="/home/nelsonlab/guayule/01_Genome/03_ReadsDepth/01_AZ2"
fastq_dir="/data/home/nelsonlab/guayule/01_Genome/03_ReadsDepth/00_Reads"
trimmed_dir="/data/home/nelsonlab/guayule/01_Genome/03_ReadsDepth/00_trim"

for i in $(cat ${fastq_dir}/sample.list);
do
  fastp \
    --in1 ${fastq_dir}/${i}_combine_1.fq.gz --out1 ${trimmed_dir}/${i}_trim_R1.fastq.gz \
    --in2 ${fastq_dir}/${i}_combine_2.fq.gz --out2 ${trimmed_dir}/${i}_trim_R2.fastq.gz \
    --thread 5 \
    --trim_front1 5 --trim_tail1 5 \
    --cut_front --cut_front_window_size 3 \
    --cut_front_mean_quality 15 --cut_tail \
    --cut_tail_window_size 3 --cut_tail_mean_quality 15 \
    --n_base_limit 5 --average_qual 10 --length_required 75 \
    --dedup --dup_calc_accuracy 6 --reads_to_process 450000000  &
done

cd  $trimmed_dir

 ls | grep -v "sample" | cut -d "_" -f1 | uniq > sample.list
```
### Mapping reads to reference genome
**AZ2 genome**
```bash
trimmed_dir="/data/home/nelsonlab/guayule/01_Genome/03_ReadsDepth/00_trim"
Ref_dir="/data/home/nelsonlab/guayule/01_Genome/AZ2"
bam_dir="/data/home/nelsonlab/guayule/01_Genome/03_ReadsDepth/01_AZ2/01_bam"

### index the genome
bwa index $Ref_dir/Parthenium_argentatum_var_AZ2_D.mainGenome.fast

### align reads by accessions
for i in $(cat $trimmed_dir/sample.list);
do
  bwa mem -M -R "@RG\tID:$i\tSM:$i\tPL:ILLUMINA" \
    -t 6  \
    $Ref_dir/Parthenium_argentatum_var_AZ2_D.mainGenome.fasta \
    ${trimmed_dir}/${i}_trim_R1.fastq.gz \
    ${trimmed_dir}/${i}_trim_R2.fastq.gz | samtools view -@10 -bS - -o ${bam_dir}/${i}.bam &

    ### check numbers of reads within each bam file
    samtools flagstat ${bam_dir}/${i}.bam > ${bam_dir}/${i}_flagstat.txt
    
    ### summarize mappeds reads into logStats file
    grep 'mapped' ${bam_dir}/${i}_flagstat.txt >> bwaMapReads.stats
done
```
**Combined genome AZ2 and Cal3 genome together**
```bash
trimmed_dir="/data/home/nelsonlab/guayule/01_Genome/03_ReadsDepth/00_trim"
Ref_dir="/data/home/nelsonlab/guayule/01_Genome/AZ2-Cal3"
bam_dir="/home/nelsonlab/guayule/01_Genome/03_ReadsDepth/03_combine/01_bam"

### align reads by accessions
for i in $(cat $trimmed_dir/sample.list);
do
  bwa mem -M -R "@RG\tID:$i\tSM:$i\tPL:ILLUMINA" \
    -t 10  \
    $Ref_dir/AZ2-Cal3-genome.fasta \
    ${trimmed_dir}/${i}_trim_R1.fastq.gz \
    ${trimmed_dir}/${i}_trim_R2.fastq.gz | samtools view -@10 -bS - -o ${bam_dir}/${i}.bam

    ### check numbers of reads within each bam file
    samtools flagstat ${bam_dir}/${i}.bam > ${bam_dir}/${i}_flagstat.txt
    
    ### summarize mappeds reads into logStats file
    grep 'mapped' ${bam_dir}/${i}_flagstat.txt >> bwaMapReads.stats
done
```
### Remove multiple-mapped reads
Description of [**sambamba usage:**](https://github.com/biod/sambamba/wiki/%5Bsambamba-view%5D-Filter-expression-syntax ) 

**AZ2 genome**
```bash
bam_dir="/data/home/nelsonlab/guayule/01_Genome/03_ReadsDepth/01_AZ2/01_bam"
uniq_dir="/data/home/nelsonlab/guayule/01_Genome/03_ReadsDepth/01_AZ2/02_uniq_bam"

for i in $(cat $bam_dir/sample.list);
do
  ### remove duplicates
  sambamba view \
    -t 2 -h -f bam \
    -F "mapping_quality >= 1 and not (unmapped or secondary_alignment or mate_is_unmapped or chimeric) and not ([XA] != null or [SA] != null) and proper_pair" \
    ${bam_dir}/${i}.bam \
    -o ${uniq_dir}/${i}.uniq.bam

  ### check numbers of reads within each bam file
  samtools flagstat ${uniq_dir}/${i}.uniq.bam > ${bam_dir}/${i}.uniq.flagstat.txt

  ### summarize mappeds reads into logStats file 
  grep 'mapped' ${uniq_dir}/${i}.uniq.flagstat.txt >> UniqMapReads.stats

  ### sort bam files
  samtools sort -m 4G ${uniq_dir}/${i}.uniq.bam > ${uniq_dir}/${i}.uniq.sort.bam

  ### check reads number after sort
  samtools flagstat ${uniq_dir}/${i}.uniq.sort.bam > ${uniq_dir}/${i}.uniqSort.flagstat.txt
done
```
**Combined genome**
```bash
bam_dir="/data/home/nelsonlab/guayule/01_Genome/03_ReadsDepth/03_combine/01_bam"
uniq_dir="/data/home/nelsonlab/guayule/01_Genome/03_ReadsDepth/03_combine/02_uniq_bam"

### list all samples
cd $bam_dir
ls | grep 'bam' | cut -d "." -f1 > sample.list

for i in $(cat $bam_dir/sample.list);
do
  ### remove low quality reads
  sambamba view \
    -t 2 -h -f bam \
    -F "mapping_quality >= 50 and not (unmapped or secondary_alignment or mate_is_unmapped or chimeric) and not ([XA] != null or [SA] != null) and proper_pair" \
    ${bam_dir}/${i}.bam \
    -o ${uniq_dir}/${i}.uniq.bam &

  for i in $(cat $bam_dir/sample.list);
  do
    ### remove duplicates
    sambamba markdup \
      -t 10 -r   \
      ${uniq_dir}/${i}.uniq.bam \
      ${uniq_dir}/${i}.clean.bam &
  done

  ### check numbers of reads within each bam file
  samtools flagstat ${uniq_dir}/${i}.uniq.bam > ${bam_dir}/${i}.uniq.flagstat.txt

  ### summarize mappeds reads into logStats file 
  grep 'mapped' ${uniq_dir}/${i}.uniq.flagstat.txt >> UniqMapReads.stats

  ### sort bam files
  samtools sort -m 4G ${uniq_dir}/${i}.uniq.bam > ${uniq_dir}/${i}.uniq.sort.bam

  ### check reads number after sort
  samtools flagstat ${uniq_dir}/${i}.uniq.sort.bam > ${uniq_dir}/${i}.uniqSort.flagstat.txt
done
```

### Calculte the coverage in sliding windows 
#### Manually prepare the bed file for each chromosomes (see example)
```
Chr01 0 100000
Chr01 20000 120000
Chr01 40000 140000
Chr01 60000 160000
Chr01 80000 180000
Chr01 100000 200000
```

```bash
### For AZ2 samples
uniq_dir="/data/home/nelsonlab/guayule/01_Genome/03_ReadsDepth/01_AZ2/02_uniq_bam"
bed_dir="/data/home/nelsonlab/guayule/01_Genome/03_ReadsDepth/01_AZ2/03_Coverage"

for i in $(cat $uniq_dir/sample.list);   
do
  samtools bedcov $bed_dir/AZ2_window.bed ${uniq_dir}/${i}.uniq.sort.bam > ${bed_dir}/${i}.AZ2-Cov.txt
done
```

```bash
### For AZ2-Cal3 combined samples
uniq_dir="/data/home/nelsonlab/guayule/01_Genome/03_ReadsDepth/03_combine/02_uniq_bam"
bed_dir="/data/home/nelsonlab/guayule/01_Genome/03_ReadsDepth/03_combine/03_Coverage"

for i in $(cat $uniq_dir/sample.list);   
do
  samtools bedcov $bed_dir/Genome_window.bed ${uniq_dir}/${i}.clean.sort.bam > ${bed_dir}/${i}.AZ2-Cov.txt
done

```

