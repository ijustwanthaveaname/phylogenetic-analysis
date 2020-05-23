#将fasta格式的序列比对转换为sto格式(利用biopython)
from Bio import SeqIO
inputfile=input("please input your fasta_file:")
outputfile=input("please input your output_file:")
records = SeqIO.parse(inputfile, "fasta")
count = SeqIO.write(records,outputfile, "stockholm")
print("Converted %i records" % count)
