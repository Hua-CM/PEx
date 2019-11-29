usage()
{
	cat <<EOF >&2
	OPTIONS:
	-s sample name
	-i STAT index
	-o output directory. MUST BE ABSOLUTE PATH.
	--tmp temporary directory. MUST BE ABSOLUTE PATH. default:./tmp
        --genome genome fasta path
	--fq1 fastq one pair
        --fq2 fastq the other pair
	-h show help of scripts
        Example:
	template.sh -i input.txt -d database.txt -o result.txt
EOF
}
bin_star=/home/rnaseq/STAR-2.7.3a/source/STAR
tmp_dir=./tmp
eval set -- `getopt -o hi:s:o: -l fq1:,fq2:,tmp:,genome: -- "$@"`
while [ $# -gt 0 ]
do
	case $1 in
		--fq1)
			fq1=$2
			shift 2
			;;
		--fq2)
			fq2=$2
			shift 2
			;;
		--tmp)
			tmp_dir=$2
			shift 2
			;;
		--genome)
			genome=$2
			shift 2
			;;
		-i)
			star_index=$2
			shift 2
			;;
		-s)
			sample=$2
			shift 2
			;;
		-o)
			out_dir=$2
			shift 2
			;;
		-h)
			usage
			exit
			;;
		--)
			shift 1
			break
			;;
		-*)
			echo "Invalid option:$1"
			exit 2
			;;
		*)
			echo "enter the right para"
			break
			;;
	esac
done
if [ -z "$fq1" ] || [ -z "$fq2" ] || [ -z "$star_index" ] || [ -z "$sample" ] || [ -z "$out_dir" ] || [ -z "$genome" ];then
	usage
	exit 2
fi
start=$(date +%s.%N)
if [ ! -d $tmp_dir ]
then
	mkdir $tmp_dir && cd $_
else
	cd $tmp_dir
fi
echo star start `date`
if [ ! -f "$genome.dict" ]; then
	        samtools faidx $genome
fi
$bin_star --runThreadN 30  --genomeDir $star_index  \
	--twopassMode Basic --outReadsUnmapped None --chimSegmentMin 12  \
	--alignIntronMax 100000 --chimSegmentReadGapMax parameter 3  \
	--alignSJstitchMismatchNmax 5 -1 5 5  \
	--readFilesIn $fq1 $fq2 --outFileNamePrefix  ${sample}
mv ${sample}_Aligned.out.sam $sample.sam
samtools sort -o $sample.bam  $sample.sam
samtools index $sample.bam
touch  ok.star.$sample.status
rm  $sample.sam
sambamba markdup --overflow-list-size 600000  --tmpdir='./'  -r  $sample.bam  ${sample}_rmd.bam
if [ ! -f "$genome.dict" ]; then
	gatk CreateSequenceDictionary -R $genome
fi
gatk SplitNCigarReads -R $genome -I ${sample}_rmd.bam -O  ${sample}_rmd_split.bam
mv ${sample}_rmd_split.bam $out_dir/${sample}_rmd_split.bam
echo star  end  `date`
dur=$(echo "$(date +%s.%N) - $start" | bc)
printf "Execution time for star : %.6f seconds" $dur
