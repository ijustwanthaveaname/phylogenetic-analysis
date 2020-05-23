#!/bin/bash
ls *.fasta|while read id
do
#TransDecoder.LongOrfs -t ${id}
#TransDecoder.Predict -t ${id}
t_coffee -other_pg seq_reformat -in ${id} -action +translate -output fasta_seq > protein_${id}
mafft protein_${id} > aln_protein_${id}
trimal -automated1 -in aln_protein_${id} -out protein_trim_${id}
pal2nal.pl protein_trim_${id} ${id} -output paml > ${id%.fasta}.paml
iqtree -s protein_trim_HBB.fasta -m MFP
#进行支位点模型
done
ls *.paml|while read id
do 
	cat branch-site.ctl > ${id%.paml}.ctl
	echo "
  	seqfile = ${id}
  	outfile = ${id}_mlc " >>  ${id%.paml}.ctl
	cat branch-site_null.ctl >${id%.paml}_null.ctl 
	echo "
  	seqfile = ${id}
	outfile = ${id}_null_mlc" >> ${id%.paml}_null.ctl  
	codeml ${id%.paml}.ctl
      	codeml ${id%.paml}_null.ctl 
	names=`cat ${id}_mlc ${id}_null_mlc|grep 'lnL'|awk '
	NR==1{alternative=$5}NR==2{null=$5}
	END{if(2*(alternative-null)>3.84){print 1}}'`
	if [ $names -eq 1 ];then
		echo ${id%.paml} >> possitive_gene.txt
	fi
done	
#cat mlc|grep 'lnL'|awk 'NR==1{M0=$5}NR==2{M1=$5}NR==3{M2=$5}NR==4{M3=$5}NR==5{M7=$5}NR==6{M8=$5}END{if (2*(M3-M0)>9.49){print "M3 is right"};if(2*(M2-m1)>5.99){print "M2 is right"};if(2*(M8-M7)>5.99){print "M8 is right"}}'
