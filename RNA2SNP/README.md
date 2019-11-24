## Step1: Generate the STAR index

~~~shell
STAR --runMode genomeGenerate --genomeDir $Your_index_dirctory --genomeFastaFiles $Your_genome_fasta --sjdbGTFfile $Your_genome_gtf --sjdbOverhang $(cat $read_length-1 | bc) --runThreadN $threads
~~~

## Step2: Generate bam file

use the domestic  scripts: STAR_align.sh

