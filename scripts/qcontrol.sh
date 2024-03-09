dir=part_2
outdir=qcontrol

mkdir -p $dir/$outdir

fastqc $dir/*.bam -o $dir/$outdir