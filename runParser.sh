#!/bin/sh
set -euo pipefail
## usage
usage() {
    echo "
usage:   run_parser.sh [options]
options:
  -d    data directory
  -o    output directory
  -c    directory containing read depth information to annotate calls with

  -f    filter
  -s    only operate on somatic tumour variants
  -m    merge
  -a    annotate
  -r    remove false positives anad reannotate
  -e    exclude bed file
  -h    show this message
"
}

filter=0
merge=0
annotate=0
replace=0
out_dir='filtered/'
data_dir='data/'
cnv_dir=
exclude_file='/Users/Nick_curie/Documents/Curie/Data/Genomes/Dmel_v6.12/Mappability/dmel6_unmappable_100.bed'
somatic=0

while getopts 'fmarsho:d:c:e:' flag; do
  case "${flag}" in
    f)  filter=1 ;;
    m)  merge=1 ;;
    a)  annotate=1 ;;
    r)  replace=1 ;;
    s)  somatic=1 ;;
    d)  data_dir="$OPTARG" ;;
    o)  out_dir="$OPTARG" ;;
    c)  cnv_dir="$OPTARG" ;;
    e)  exclude_file="$OPTARG" ;;
    h)  usage
        exit 0 ;;
  esac
done

if [[ $# -eq 0 ]]
then
  usage
  exit 0
fi

# Change to the 'script' dir in svParser
dir=$(dirname "$0")
script_bin="$dir/script"
script_bin="$( cd "$script_bin" ; pwd -P )"

echo "Reading data from '$data_dir'"
echo "Writing data to '$out_dir'"
echo "Exclude file set to '$exclude_file'"

mkdir -p "$out_dir/summary"
s=''
annoS=''
# Run svParser for each type of variant file
if [[ $filter -eq 1 ]]
then
  if [[ $somatic -eq 1 ]]
  then
    s='-f st=1'
    annoS='-s'
  fi

  echo "**************************"
  echo "*** Filtering variants ***"
  echo "**************************"

  for lumpy_file in $data_dir/lumpy/*.vcf
  do
    [ -f $lumpy_file ] || continue
    echo "perl "$script_bin"/svParse.pl -v $lumpy_file -m l -f chr=1 -f su=3 -f dp=10 -f sq=0.1 -f rdr=0.05 $s -e $exclude_file -o $out_dir -p"
    perl "$script_bin"/svParse.pl -v $lumpy_file -m l -f chr=1 -f su=3 -f dp=10 -f sq=0.1 -f rdr=0.05 $s -e $exclude_file -p -o $out_dir
  done

  for delly_file in $data_dir/delly/*.vcf
  do
  [ -f $delly_file ] || continue
    echo "perl "$script_bin"/svParse.pl -v $delly_file -m d -f chr=1 -f su=3 -f dp=10 -f sq=0.1 -f rdr=0.05 $s -e $exclude_file -o $out_dir -p"
    perl "$script_bin"/svParse.pl -v $delly_file -m d -f chr=1 -f su=3 -f dp=10 -f sq=0.1 -f rdr=0.05 $s -e $exclude_file -p -o $out_dir
done

  for novo_file in $data_dir/novobreak/*.vcf
  do
    [ -f $novo_file ] || continue
    echo "perl "$script_bin"/svParse.pl -v $novo_file -m n -f chr=1 -f su=3 -f dp=10 -f sq=0.1 -f rdr=0.05 $s -e $exclude_file -o $out_dir -p"
    perl "$script_bin"/svParse.pl -v $novo_file -m n -f chr=1 -f su=3 -f dp=10 -f sq=0.1 -f rdr=0.05 $s -e $exclude_file -p -o $out_dir
  done

  for freec_file in $data_dir/freec/*filt_cnvs.txt
  do
    [ -f $freec_file ] || continue
    echo "perl "$script_bin"/parseCF.pl -c $freec_file -o $out_dir/summary"
    perl "$script_bin"/parseCF.pl -c $freec_file -o $out_dir/summary
  done

  for cnv_file in $data_dir/cnv/*.txt
  do
    [ -f $cnv_file ] || continue
    echo "perl "$script_bin"/parseCNV.pl -c $cnv_file -o $out_dir/summary"
    perl "$script_bin"/parseCNV.pl -c $cnv_file -o $out_dir/summary
  done

fi

# If CNV-Seq has been run, the cnv directory can be specified with the -c flag
# For each summary file, annotate somatic events with log2(FC) from .cnv file
if [[ -n "$cnv_dir" && $filter -eq 1 ]]
then
  echo "*******************************************************"
  echo "*** Annotating variants with read depth information ***"
  echo "*******************************************************"

  cd $out_dir/summary
  samples+=( $(ls -1 *.filtered.summary.txt | cut -d '.' -f 1 | sort -u ) )

  for ((i=0;i<${#samples[@]};++i))
  do
    if [ ! -f $cnv_dir/${samples[i]}.*.cnv ]
    then
      echo " -> ! No corresponding CNV file for ${samples[i]} in $cnv_dir"
    else
      echo "perl "$script_bin"/findCNV.pl -c $cnv_dir/${samples[i]}.*.cnv -v $out_dir/summary/${samples[i]}*.filtered.summary.txt"
      perl "$script_bin"/findCNV.pl -c $cnv_dir/${samples[i]}.*.cnv -v $out_dir/summary/${samples[i]}*.filtered.summary.txt
    fi
  done
fi

cd $out_dir

if [[ $merge -eq 1 ]]
then
  echo "************************"
  echo "*** Merging variants ***"
  echo "************************"

  mergeVCF=`which mergevcf || true`
  if [[ -z "$mergeVCF" ]]
  then
    usage
    echo -e "Error: mergevcf was not found. Please set in path\n`pip install mergevcf`"
    exit 1
  fi
  echo "perl "$script_bin"/merge_vcf.pl"
  #perl "$script_bin"/merge_vcf.pl
fi

cd $out_dir/summary

if [[ $merge -eq 1 ]]
then
  mkdir -p "$out_dir/summary/merged/"
  if [[ $cnv_dir ]]
  then
    samples+=( $(ls -1 *.summary.cnv.txt | cut -d '.' -f 1 | sort -u ) )

    for ((i=0;i<${#samples[@]};++i))
    do
      echo "perl "$script_bin"/svMerger.pl -f ${samples[i]}.*.summary.cnv.txt"
      perl "$script_bin"/svMerger.pl -f ${samples[i]}.*.summary.cnv.txt -o "$out_dir/summary/merged"
    done
  else
    samples+=( $(ls -1 *.summary.txt | cut -d '.' -f 1 | sort -u ) )

    for ((i=0;i<${#samples[@]};++i))
    do
      echo "perl "$script_bin"/svMerger.pl -f ${samples[i]}.*.summary.txt"
      perl "$script_bin"/svMerger.pl -f ${samples[i]}.*.summary.txt -o "$out_dir/summary/merged"
    done
  fi
fi

cd $out_dir/summary/merged

if [[ $merge -eq 1 ]]
then
  for f in *_merged_SVs.txt
  do
    echo "perl "$script_bin"/svClusters.pl -v $f -d 500"
    perl "$script_bin"/svClusters.pl -v $f -d 500
    rm $f
  done
fi

#features=/Users/Nick/Documents/Curie/Data/Genomes/Dmel_v6.12/Features/dmel-all-r6.12.gtf # home
features=/Users/Nick_curie/Documents/Curie/Data/Genomes/Dmel_v6.12/Features/dmel-all-r6.12.gtf # work
blacklist=${out_dir}/summary/merged/all_samples_blacklist.txt
whitelist=${out_dir}/summary/merged/all_samples_whitelist.txt

if [[ $annotate -eq 1 ]]
then

  echo "***************************"
  echo "*** Annotating variants ***"
  echo "***************************"

  cd $out_dir/summary/merged

  if [ -f "all_genes.txt" ] && [ -f "all_bps.txt" ]
  then
    rm "all_genes.txt"
    rm "all_bps.txt"
  fi

  for clustered_file in *clustered_SVs.txt
  do
    echo "Annotating $clustered_file"
    # Should check both files individually
    if [ -f $blacklist ] && [ -f $whitelist ]
    then
      echo "perl "$script_bin"/sv2gene.pl -f $features -i $clustered_file -b $blacklist -w $whitelist $annoS"
      perl "$script_bin"/sv2gene.pl -f $features -i $clustered_file -b $blacklist -w $whitelist $annoS
    elif [ -f $blacklist ]
    then
      echo "perl "$script_bin"/sv2gene.pl -f $features -i $clustered_file -b $blacklist $annoS"
      perl "$script_bin"/sv2gene.pl -f $features -i $clustered_file -b $blacklist $annoS
    elif [ -f $whitelist ]
    then
      echo "perl "$script_bin"/sv2gene.pl -f $features -i $clustered_file -w $whitelist $annoS"
      perl "$script_bin"/sv2gene.pl -f $features -i $clustered_file -w $whitelist $annoS
    else
      echo "perl "$script_bin"/sv2gene.pl -f $features -i $clustered_file $annoS"
      perl "$script_bin"/sv2gene.pl -f $features -i $clustered_file $annoS
    fi
    rm $clustered_file
  done

fi

if [[ $replace -eq 1 ]]
then
    # echo "Adding any new CNV calls to data/cnv'"
    # for annofile in *_annotated_SVs.txt
    # do
    #   python "$script_bin"/getCNVs.py -f $annofile
    # done

  if [ -f "all_genes_filtered.txt" ] && [ -f "all_bps_filtered.txt" ]
  then
    rm "all_genes_filtered.txt"
    rm "all_bps_filtered.txt"
  fi

  echo "Removing calls marked as false positives in 'all_samples_blacklist.txt'"
  for annofile in *_stitched.txt
  # do
  #   perl "$script_bin"/clean_files.pl -v $annofile -o $out_dir/summary/merged -b $blacklist -w $whitelist
  # done

  #
  # for clean_file in *cleaned_SVs.txt
  do
    # Delete file if empty
    # if [[ ! -s $annofile ]]
    # then
    #   rm $annofile
    # else
      # Annotate un-annotated (manually added) calls
      # Append any new hit genes to 'all_genes.txt'
      perl "$script_bin"/sv2gene.pl -r -f $features -i $annofile -s
      # rm $clean_file
    # fi
  done

  echo "Writing bp info for cleaned, reannotated SV calls to 'all_bps_filtered.txt'"

  # for reanno_file in *reannotated_SVs.txt
  # do
  #   # Grab some of the fields from these newly annotated files, and write them to 'all_bps_cleaned.txt'
  #   python "$script_bin"/getbps.py -f $reanno_file
  # done
  #
  # # This shouldn't be neccessary. All calls in this file are taken from 'reannotated' files, which should have FP removed already...
  # echo "Removing false positives from bp file 'all_bps_cleaned.txt', writing new bp file to 'all_bps_filtered.txt'"
  # python "$script_bin"/update_bps.py
  # rm 'all_bps_cleaned.txt'

  # Merge all samples
  cd $out_dir/summary/merged
  echo "Merging all samples into single file..."
  perl "$script_bin"/merge_files.py
fi


function getBase(){
  stem=$(basename "$1" )
  name=$(echo $stem | cut -d '.' -f 1)
  echo $name
}


exit 0
