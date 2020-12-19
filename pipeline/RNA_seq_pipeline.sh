#!/bin/bash

#set software path used in this pipeline
star_bin=/home/rnaseq/STAR-2.7.3a/bin/Linux_x86_64/STAR

usage()
{
	cat <<EOF >&2
	OPTIONS:
	-n project name
	-f reference genome fasta
	-g reference genome gff
	-i input directory. MUST BE ABSOLUTE PATH.
	-o output directory. MUST BE ABSOLUTE PATH.
	-p threads default:4
	--tmp temporary directory. MUST BE ABSOLUTE PATH. default:./tmp
	-h show help of scripts
        Example:
	RNA_seq_pipeline.sh -i input.txt -d database.txt -o result.txt
EOF
}
tmp_dir=./tmp
threads=4
eval set -- `getopt -o hn:f:g:i:o:p: -l tmp -- "$@"`
while [ $# -gt 0 ]
do
	case $1 in
		-n)
			project_name=$2
			shift 2
			;;
		-f)
			gfasta=$2
			shift 2
			;;
		-g)
			ggff=$2
			shift 2
			;;
		-i)
			in_dir=$2
			shift 2
			;;
		-o)
			out_dir=$2
			shift 2
			;;
		-p)
			threads=$2
			shift 2
			;;
		--tmp)
			tmp_dir=$2
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
if [ -z "$project_name" ] || [ -z "$gfasta" ] || [ -z "$ggff" ] || [ -z "$in_dir" ] || [ -z "$out_dir" ];then
	usage
	exit 2
fi
start=$(date +%s.%N)
if [ ! -d $out_dir ];then
	mkdir $out_dir
fi
if [ ! -d $tmp_dir ];then
	mkdir $tmp_dir
fi
echo start `date`
# create index
gname=${echo "$ggff"  | awk -F "/" '{print $NF}' | awk -F "." '{print $1}'}
if [ ! -d $out_dir ]
then
	mkdir $out_dir/01.supp
fi
cp $gfasta $out_dir/01.supp
cp $ggff $out_dir/01.supp
gfasta=`realpath $out_dir/01.supp/${gfasta##*/}`
ggff=`realpath $out_dir/01.supp/${ggff##*/}`
## hisat2
hisat2-build -p $threads $gfasta $out_dir/01.supp/${gname}_hisat
## RSEM
rsem-prepare-reference --gff  $ggff \
                    --star \
                    --star-path $star_bin \
                    -p $threads \
                    $gfasta \
                    $out_dir/01.supp/${gname}_rsem
## RSeQC
./gtf2bed $out_dir/01.supp/${gname}.gtf > $out_dir/01.supp/$gname.bed

# map and quantify for each sample
## make directory
mkdir $out_dir/02.data
mkdir $out_dir/03.data

for sample_ in $in_dir/*
do
	s_name=$(echo "${sample_}"| awk -F "/" '{print $NF}')
	mkdir $out_dir/$s_name
	echo start "$s_name"
	# If there is more than one SRR* file in a SRS file, run the following two lines.
	if [ ${ls "$sample_" | wc -l } -gt 2 ];then
		cat ${sample_}/*_1.fastq.gz >> $tmp_dir/${s_name}_1.fastq.gz
		cat ${sample_}/*_2.fastq.gz >> $tmp_dir/${s_name}_2.fastq.gz
		 fastp -q 20 -u 40 -l 50 -g -x -r -W 4 -M 20 -w $threads \
                                        -i $tmp_dir/${s_name}_1.fastq.gz -o $out_dir/02.data/${s_name}_clean_1.fastq.gz \
                                        -I $tmp_dir/${s_name}_2.fastq.gz -O $out_dir/02.data/${s_name}_clean_2.fastq.gz \
                                        -h $out_dir/02.data${s_name}_report.html > $out_dir/02.data/${s_name}_fastp_report.txt
	else
		fastp -q 20 -u 40 -l 50 -g -x -r -W 4 -M 20 -w $threads \
                                        -i ${sample_}/*_1.fastq.gz -o $out_dir/02.data/${s_name}_clean_1.fastq.gz \
                                        -I ${sample_}/*_2.fastq.gz -O $out_dir/02.data/${s_name}_clean_2.fastq.gz \
                                        -h $out_dir/02.data${s_name}_report.html > $out_dir/02.data/${s_name}_fastp_report.txt
	fi
	## determine strand
	hisat2 -p $threads -x $out_dir/01.supp/${gname}_hisat \
		-1 $out_dir/02.data/${s_name}_clean_1.fastq.gz \
		-2 $out_dir/02.data/${s_name}_clean_2.fastq.gz \
		-S $out_dir/03.mapping/${s_name}_alignment_unsorted.sam  2> $out_dir/03.mapping/${s_name}_hisat2_Mapping_Rate.txt
	samtools view -b -@ $1 -S $out_dir/03.mapping/${s_name}_alignment_unsorted.sam -o $out_dir/03.mapping/${s_name}_alignment_unsorted.bam
	
	/home/anaconda3/envs/python27/bin/python \
	/home/rnaseq/RSeQC-2.6.4/scripts/infer_experiment.py \
	-r $out_dir/01.supp/$gname.bed -i $out_dir/03.mapping/${s_name}_alignment_unsorted.bam > $out_dir/03.mapping/${s_name}_RseQC.txt
	#6.According to RseQC.txt, confirm the value of "--forward-prob"
	parameter_1=`sed -n 5p $out_dir/03.mapping/${s_name}_RseQC.txt | awk '{print $7}'`
	parameter_2=`sed -n 6p $out_dir/03.mapping/${s_name}_RseQC.txt | awk '{print $7}'`
	decimal=`awk -v x="$parameter_1" -v y="$parameter_2" 'BEGIN{printf "%.0f\n",(x-y)*100}'`
	if [ $decimal -gt 20 ]; then
	    value=1
	elif [ $decimal -gt -20 ]; then
	    value=0.5
	else
	    value=0
	fi
	/home/rnaseq/RSEM-1.3.1/rsem-calculate-expression \
	--keep-intermediate-files --temporary-folder $out_dir/03.mapping/${s_name}_STAR_rsem \
	--paired-end --forward-prob=$value -p $threads --time \
	--star --star-path /home/rnaseq/STAR-2.7.3a/bin/Linux_x86_64/ \
	--append-names --output-genome-bam \
	--sort-bam-by-coordinate \
	--star-gzipped-read-file $out_dir/02.data/${s_name}_clean_1.fastq.gz $out_dir/02.data/${s_name}_clean_2.fastq.gz \
	$out_dir/01.supp/${gname}_rsem \
	$out_dir/03.mapping/${s_name}_star_RSEM
done
# quantative
mkdir $out_dir/04.quan
/home/rnaseq/RSEM-1.3.1/rsem-generate-data-matrix $out_dir/03.mapping/*.genes.results > $out_dir/04.quan${project_name}_GeneMat_rawCounts.txt
/home/rnaseq/RSEM-1.3.1/rsem-generate-data-matrix $out_dir/03.mapping/*.isoforms.results > $out_dir/04.quan/${project_name}_TransMat_rawCounts.txt

/home/rnaseq/RSEM-1.3.1/rsem-generate-data-matrix-FPKM $out_dir/03.mapping/*.genes.results > $out_dir/04.quan/${project_name}_GeneMat_FPKM.txt
/home/rnaseq/RSEM-1.3.1/rsem-generate-data-matrix-FPKM $out_dir/03.mapping/*.isoforms.results > $out_dir/04.quan/${project_name}_TransMat_FPKM.txt

/home/rnaseq/RSEM-1.3.1/rsem-generate-data-matrix-TPM $out_dir/03.mapping/*.genes.results > $out_dir/04.quan/${project_name}_GeneMat_TPM.txt
/home/rnaseq/RSEM-1.3.1/rsem-generate-data-matrix-TPM $out_dir/03.mapping/*.isoforms.results > $out_dir/04.quan/${project_name}_TransMat_TPM.txt