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

