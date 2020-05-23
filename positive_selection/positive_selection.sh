#!/bin/bash
echo "
*seqfile = paml   * sequence data file name
     treefile = $1 * tree structure file name

      *outfile = mlc          * main result file name
        noisy = 3   * 0,1,2,3,9: how much rubbish on the screen
      verbose = 0   * 1: detailed output, 0: concise output
      runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic
                    * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

      seqtype = 1   * 1:codons; 2:AAs; 3:codons-->AAs
    CodonFreq = 2   * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
        clock = 0   * 0: no clock, unrooted tree, 1: clock, rooted tree
       aaDist = 0   * 0:equal, +:geometric; -:linear, {1-5:G1974,Miyata,c,p,v}
        model = 2	* 0:one,1:b,2:2 or more dN/dS ratios for branches

      NSsites = 2
                    * 0:one w; 1:NearlyNeutral; 2:PositiveSelection; 3:discrete;
                    * 4:freqs; 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;10:3normal
        icode = 0   * 0:standard genetic code; 1:mammalian mt; 2-10:see below
        Mgene = 0   * 0:rates, 1:separate; 2:pi, 3:kappa, 4:all

    fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
        kappa = .3   * initial or fixed kappa
    fix_omega = 0   * 1: omega or omega_1 fixed, 0: estimate 
        omega = 2  * initial or fixed omega, for codons or codon-based AAs
        ncatG = 5   * # of categories in the dG or AdG models of rates

        getSE = 0   * 0: don't want them, 1: want S.E.s of estimates
 RateAncestor = 0   * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)

   Small_Diff = .45e-6
    cleandata = 1  * remove sites with ambiguity data (1:yes, 0:no)?
  fix_blength = 0  * 0: ignore, -1: random, 1: initial, 2: fixed, 3: proportional" >branch-stie.ctl
echo "*seqfile = paml  * sequence data file name
     treefile = protein_trim_HBB.fasta.treefile.Labeled.txt  * tree structure file name

      *outfile = mlc_null          * main result file name
        noisy = 3   * 0,1,2,3,9: how much rubbish on the screen
      verbose = 0   * 1: detailed output, 0: concise output
      runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic
                    * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

      seqtype = 1   * 1:codons; 2:AAs; 3:codons-->AAs
    CodonFreq = 2   * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
        clock = 0   * 0: no clock, unrooted tree, 1: clock, rooted tree
       aaDist = 0   * 0:equal, +:geometric; -:linear, {1-5:G1974,Miyata,c,p,v}
        model = 2	* 0:one,1:b,2:2 or more dN/dS ratios for branches

      NSsites = 2
                    * 0:one w; 1:NearlyNeutral; 2:PositiveSelection; 3:discrete;
                    * 4:freqs; 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;10:3normal
        icode = 0   * 0:standard genetic code; 1:mammalian mt; 2-10:see below
        Mgene = 0   * 0:rates, 1:separate; 2:pi, 3:kappa, 4:all

    fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
        kappa = .3   * initial or fixed kappa
    fix_omega = 1   * 1: omega or omega_1 fixed, 0: estimate 
        omega = 1  * initial or fixed omega, for codons or codon-based AAs
        ncatG = 5   * # of categories in the dG or AdG models of rates

        getSE = 0   * 0: don't want them, 1: want S.E.s of estimates
 RateAncestor = 0   * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)

   Small_Diff = .45e-6
    cleandata = 1  * remove sites with ambiguity data (1:yes, 0:no)?
  fix_blength = 0  * 0: ignore, -1: random, 1: initial, 2: fixed, 3: proportiona" > branch-site_null.ctl
#ls *.fasta|while read id
#do
#TransDecoder.LongOrfs -t ${id}
#TransDecoder.Predict -t ${id}
#t_coffee -other_pg seq_reformat -in ${id} -action +translate -output fasta_seq > protein_${id}
#mafft protein_${id} > aln_protein_${id}
#trimal -automated1 -in aln_protein_${id} -out protein_trim_${id}
#pal2nal.pl protein_trim_${id} ${id} -output paml > ${id%.fasta}.paml
#iqtree -s protein_trim_HBB.fasta -m MFP
#进行支位点模型
#done
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
