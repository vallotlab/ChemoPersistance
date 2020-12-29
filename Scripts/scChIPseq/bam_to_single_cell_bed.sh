## Mapping Research Pipeline
## Copyleft 2018 Institut Curie
## Author(s): Nicolas Servant, Pacôme Prompsy
## Contact:  pacome.prompsy@curie.fr
## This software is distributed without any guarantee under the terms of the CECILL License
## See the LICENCE file for details

## bam_to_sc_bed 
# Generates single-cell BED files (one entry per read) from barcode-tagged BAM file (XB tag)
# Only generate bed for cells with > {count} reads (recommended 1000)
# Creates a directory with the single cell files 
# Increase number of threads for faster conversion (~6 threads ~ 8min/gb for th.=1000 reads) 
# Require samtools 1.3
# Usage :
# bam_to_sc_bed {BAM} ${minimum read per cell} {path to output dir} {path to log dir} {number of threads} 
# Exemple:
# bam_to_sc_bed MM468_DMSO3_day77.bam 1000 ./ 6
bam_to_sc_bed() {
  bam_in=$1
  count=$2
  local odir=$3
  NB_PROC=$4

  local prefix=${odir}/$(basename $bam_in | sed -e 's/.bam$//')
  
  echo -e "Creating sc bed from Mapped & Dedup BAM ..."
  echo
  
  mkdir -p ${odir}
  mkdir -p ${odir}/scBed_${count}
  
  #Get barcode field & read length
  barcode_field=$(samtools view $bam_in  | sed -n "1 s/XB.*//p" | sed 's/[^\t]//g' | wc -c)

  echo -e "Sorting BAM by barcode  ..."
  echo
  #Create header
  samtools view -H $bam_in | sed '/^@HD/ d' > ${prefix}_tmp_header.sam
    
  #Sort by Barcode, Chr, Pos R1 :
  samtools view $bam_in | LC_ALL=C sort --parallel=${NB_PROC} -t $'\t' -k "$barcode_field.8,$barcode_field"n -k 3.4,3g -k 4,4n >> ${prefix}_tmp_header.sam
    
  samtools view -@ ${NB_PROC} -b ${prefix}_tmp_header.sam > ${prefix}_tmp.sorted.bam
  
  echo -e "Create single-cell bed for each cell with > $count reads ..."
  echo
  #Convert to scBed
  samtools view ${prefix}_tmp.sorted.bam | awk -v odir=${odir}/scBed -v bc_field=$barcode_field -v OFS="\t" -v count=$count '
  BEGIN{
    split(count,min_counts,",")
  }
  NR==1{
    lastBC=substr($bc_field,6,15);
    i=1
    chr[i] = $3
    start[i] = $4
    end[i] = $4 +1
  }
  NR>1{
  if(lastBC==substr($bc_field,6,15)){
    i = i +1
    chr[i] = $3
    start[i] = $4
    end[i] = $4 +1
    }
    else{
    for(m=1; m<=length(min_counts);m++){
      if(i > min_counts[m]){
        for (x=1; x<=i; x++){
          out = odir"_"min_counts[m]"/"lastBC".bed"
          print chr[x],start[x],end[x] >> out
        }
      }
    }
    i=0
    }
     lastBC=substr($bc_field,6,15);
}
'

  #Gzip
  echo -e "Gzipping single cell BED files ..."
  echo
  for i in $odir/scBed_${count}/*.bed; do gzip -9 $i; done
  mv $odir/scBed_${count} ${prefix}_scBed_${count}
  
  rm -f ${prefix}_tmp_header.sam ${prefix}_tmp.sorted.bam
 
  echo -e "Done !"
  echo
}

