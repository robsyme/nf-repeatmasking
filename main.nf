#!/usr/bin/env nextflow
VERSION = "1.0.0"

params.reference = 'data/example/genome.fasta'
params.trnaprot = 'http://www.hrt.msu.edu/uploads/535/78637/Tpases020812.gz'
params.trnanuc = 'http://gtrnadb2009.ucsc.edu/download/tRNAs/eukaryotic-tRNAs.fa.gz'
params.outdir = 'output'

trnanuc = file(params.trnanuc)
trnaprot = file(params.trnaprot)
reference = file(params.reference)

log.info ""
log.info "R E P E A T M A S K I N G  ~  version " + VERSION
log.info "reference input    : ${params.reference}"
log.info "output directory   : ${params.outdir}"
log.info ""

process recentLTRs {
  container 'robsyme/nf-repeatmasking'

  input:
  file 'genome.fasta' from reference
  file 'eukaryotic-tRNAs.fasta.gz' from trnanuc

  output:
  set age, 'seqfile.result' into ltrHarvestNew
  set age, 'seqfile.outinner' into ltrInnerSeqNew
  set age, 'CRL_Step2_Passed_Elements.fasta', 'Repeat_down*.fasta', 'Repeat_up*.fasta' into recentLTRs

  script:
  age = 'new'
  """
gt suffixerator -db genome.fasta -indexname genome.fasta -tis -suf -lcp -des -ssp -dna

gt ltrharvest \
 -index genome.fasta \
 -out seqfile.out \
 -outinner seqfile.outinner \
 -gff3 seqfile.gff \
 -minlenltr 100 \
 -maxlenltr 6000 \
 -mindistltr 1500 \
 -maxdistltr 25000 \
 -mintsd 5 \
 -maxtsd 5 \
 -motif tgca \
 -similar 99 \
 -vic 10 \
> seqfile.result

gt gff3 \
 -sort seqfile.gff \
> seqfile.gff.sort

zcat eukaryotic-tRNAs.fasta.gz > eukaryotic-tRNAs.fasta

gt ltrdigest \
 -trnas eukaryotic-tRNAs.fasta \
 seqfile.gff.sort \
 genome.fasta \
> seqfile.gff.dgt

CRL_Step1.pl \
 --gff seqfile.gff.dgt

CRL_Step2.pl \
 --step1 CRL_Step1_Passed_Elements.txt \
 --repeatfile seqfile.out \
 --resultfile seqfile.result \
 --sequencefile genome.fasta \
 --removed_repeats CRL_Step2_Passed_Elements.fasta
  """
}

process olderLTRs {
  container 'robsyme/nf-repeatmasking'

  input:
  file 'genome.fasta' from reference
  file 'eukaryotic-tRNAs.fasta.gz' from trnanuc

  output:
  set age, 'seqfile.result' into ltrHarvestOld
  set age, 'seqfile.outinner' into ltrInnerSeqOld
  set age, 'CRL_Step2_Passed_Elements.fasta', 'Repeat_down*.fasta', 'Repeat_up*.fasta' into olderLTRs

  script:
  age = 'old'
  """
gt suffixerator -db genome.fasta -indexname genome.fasta -tis -suf -lcp -des -ssp -dna

gt ltrharvest \
 -index genome.fasta \
 -out seqfile.out \
 -outinner seqfile.outinner \
 -gff3 seqfile.gff \
 -minlenltr 100 \
 -maxlenltr 6000 \
 -mindistltr 1500 \
 -maxdistltr 25000 \
 -mintsd 5 \
 -maxtsd 5 \
 -vic 10 \
> seqfile.result

gt gff3 \
 -sort seqfile.gff \
> seqfile.gff.sort

zcat eukaryotic-tRNAs.fasta.gz > eukaryotic-tRNAs.fasta

gt ltrdigest \
 -trnas eukaryotic-tRNAs.fasta \
 seqfile.gff.sort \
 genome.fasta \
> seqfile.gff.dgt

CRL_Step1.pl \
 --gff seqfile.gff.dgt

CRL_Step2.pl \
 --step1 CRL_Step1_Passed_Elements.txt \
 --repeatfile seqfile.out \
 --resultfile seqfile.result \
 --sequencefile genome.fasta \
 --removed_repeats CRL_Step2_Passed_Elements.fasta
  """
}

ltrHarvestNew
.tap { ltrHarvestResultsNew }
.set { ltrHarvestResultsForExamplarNew }

ltrInnerSeqNew
.tap { ltrHarvestInnerNew }
.set { outinnerForBlastXNew }

ltrHarvestOld
.tap { ltrHarvestResultsOld }
.set { ltrHarvestResultsForExamplarOld }

ltrInnerSeqOld
.tap { ltrHarvestInnerOld }
.set { outinnerForBlastXOld }

ltrs = recentLTRs.mix(olderLTRs)
ltrHarvestResults = ltrHarvestResultsOld.mix(ltrHarvestResultsNew)
ltrHarvestInner = ltrHarvestInnerOld.mix(ltrHarvestInnerNew)
outinnerForBlastX = outinnerForBlastXOld.mix(outinnerForBlastXNew)
ltrHarvestResultsForExamplar = ltrHarvestResultsForExamplarOld.mix(ltrHarvestResultsForExamplarNew)

process CRL_Step3 {
  container 'robsyme/nf-repeatmasking'
  tag { age }

  input:
  set age, 'CRL_Step2_Passed_Elements.fasta', 'Repeat_down*.fasta', 'Repeat_up*.fasta' from ltrs

  output:
  set age, 'CRL_Step3_Passed_Elements.fasta' into step3Passed
  set age, 'CRL_Step3_Passed_Elements.fasta' into step3PassedForExamplars

  """
CRL_Step3.pl \
 --directory . \
 --step2 CRL_Step2_Passed_Elements.fasta \
 --pidentity 60 \
 --seq_c 25
  """
}

ltrHarvestResults
.combine(step3Passed, by: 0)
.set { nestedInput }

process identifyNestedInsetions {
  container 'robsyme/nf-repeatmasking'
  tag { age }
  input:
  file 'genome.fasta' from reference
  set age, 'seqfile.result', 'CRL_Step3_Passed_Elements.fasta' from nestedInput

  output:
  set age, 'repeats_to_mask_LTR.fasta' into repeatsToMaskLTR

  """
ltr_library.pl \
 --resultfile seqfile.result \
 --step3 CRL_Step3_Passed_Elements.fasta \
 --sequencefile genome.fasta
cat lLTR_Only.lib \
| awk 'BEGIN {RS = ">" ; FS = "\\n" ; ORS = ""} \$2 {print ">"\$0}' \
> repeats_to_mask_LTR.fasta
  """
}

process RepeatMasker1 {
  container 'robsyme/nf-repeatmasking'
  tag { age }

  input:
  set age, 'repeats_to_mask_LTR.fasta', 'seqfile.outinner' from repeatsToMaskLTR.combine(ltrHarvestInner, by: 0)

  output:
  set age, 'seqfile.outinner.out', 'seqfile.outinner.masked' into repeatMasker1Unclean

  """
RepeatMasker \
 -lib repeats_to_mask_LTR.fasta \
 -nolow \
 -no_is \
 -dir . \
 seqfile.outinner

if [ ! -f seqfile.outinner.masked ]; then
  cp seqfile.outinner seqfile.outinner.masked
fi
  """
}

process cleanRM {
  tag { age }

  input:
  set age, 'seqfile.outinner.out', 'seqfile.outinner.masked' from repeatMasker1Unclean

  output:
  set age, 'seqfile.outinner.clean' into repeatMasker1Clean

  """
cleanRM.pl seqfile.outinner.out seqfile.outinner.masked > seqfile.outinner.unmasked
rmshortinner.pl seqfile.outinner.unmasked 50 > seqfile.outinner.clean
  """
}

process blastX {
  container 'robsyme/nf-repeatmasking'
  tag { age }
  cpus 4

  input:
  file 'Tpases020812DNA.fasta.gz' from trnaprot
  set age, 'seqfile.outinner.clean', 'seqfile.outinner' from repeatMasker1Clean.combine(outinnerForBlastX, by: 0)

  output:
  set age, 'passed_outinner_sequence.fasta' into blastxPassed

  """
zcat Tpases020812DNA.fasta.gz > Tpases020812DNA.fasta
makeblastdb -in Tpases020812DNA.fasta -dbtype prot
blastx \
 -query seqfile.outinner.clean \
 -db Tpases020812DNA.fasta \
 -evalue 1e-10 \
 -num_descriptions 10 \
 -num_threads ${task.cpus} \
 -out seqfile.outinner.clean_blastx.out.txt

outinner_blastx_parse.pl \
 --blastx seqfile.outinner.clean_blastx.out.txt \
 --outinner seqfile.outinner

if [ ! -s passed_outinner_sequence.fasta ]; then
  echo -e '>dummy empty sequence\nACTACTAC' > passed_outinner_sequence.fasta
fi
  """
}

blastxPassed
.combine(step3PassedForExamplars, by: 0)
.combine(ltrHarvestResultsForExamplar, by: 0)
.set { forExamplarBuilding }

process buildExemplars {
  container 'robsyme/nf-repeatmasking'
  tag { age }
  cpus 4

  input:
  file 'genome.fasta' from reference
  set age, 'passed_outinner_sequence.fasta', 'CRL_Step3_Passed_Elements.fasta', 'seqfile.result' from forExamplarBuilding

  output:
  set age, 'LTR.lib' into exemplars

  """
CRL_Step4.pl \
 --step3 CRL_Step3_Passed_Elements.fasta \
 --resultfile seqfile.result \
 --innerfile passed_outinner_sequence.fasta \
 --sequencefile genome.fasta

for lib in lLTRs_Seq_For_BLAST.fasta Inner_Seq_For_BLAST.fasta; do
  makeblastdb -in \$lib -dbtype nucl
  blastn \
   -query \${lib} \
   -db \${lib} \
   -evalue 1e-10 \
   -num_threads ${task.cpus} \
   -num_descriptions 1000 \
   -out \${lib}.out
done

CRL_Step5.pl \
 --LTR_blast lLTRs_Seq_For_BLAST.fasta.out \
 --inner_blast Inner_Seq_For_BLAST.fasta.out \
 --step3 CRL_Step3_Passed_Elements.fasta \
 --final LTR.lib \
 --pcoverage 90 \
 --pidentity 80
  """
}

newLTRs = Channel.create()
oldLTRs = Channel.create()

exemplars
.choice( newLTRs, oldLTRs ) { it[0] == "new" ? 0 : 1 }

process removeDuplicates {
  container 'robsyme/nf-repeatmasking'

  input:
  set _, 'ltrs.new.fasta' from newLTRs
  set _, 'ltrs.old.fasta' from oldLTRs

  output:
  set 'ltrs.old.fasta.masked', 'ltrs.new.fasta' into bothLTRforMasking

  "RepeatMasker -lib ltrs.new.fasta -dir . ltrs.old.fasta"
}

process filterOldLTRs {
  container 'robsyme/nf-repeatmasking'

  input:
  set 'ltrs.old.fasta.masked', 'ltrs.new.fasta' from bothLTRforMasking

  output:
  file 'allLTRs.fasta' into allLTR

  """
remove_masked_sequence.pl \
 --masked_elements ltrs.old.fasta.masked \
 --outfile ltrs.old.final.fasta
cat ltrs.new.fasta ltrs.old.final.fasta > allLTRs.fasta
  """
}

allLTR
.splitFasta(record: [id: true, sequence: true ])
.collectFile( name: 'allLTRs.fasta' ) { ">" + it.id + "#LTR\n" + it.sequence }
.tap { allLTR2 }
.set { inputForRM2 }

process RepeatMasker2 {
  container 'robsyme/nf-repeatmasking'
  cpus 10

  input:
  file 'genome.fasta' from reference
  file 'allLTR.lib' from inputForRM2

  output:
  file 'genome.fasta.masked' into genomeLtrMasked

  """
RepeatMasker \
 -no_is \
 -nolow \
 -pa ${task.cpus} \
 -lib allLTR.lib \
 -dir . \
 genome.fasta
  """
}

process RepeatModeler {
  container 'repeats'
  cpus 4

  input:
  file 'genome.masked' from genomeLtrMasked

  output:
  file 'consensi.fa.classified' into rmOutput

  """
rmaskedpart.pl genome.masked 50 > umseqfile
BuildDatabase -name umseqfiledb -engine ncbi umseqfile
RepeatModeler -pa ${task.cpus} -database umseqfiledb >& umseqfile.out
mv RM*/consensi.fa.classified consensi.fa.classified
  """
}

identityUnknown = Channel.create()
identityKnown = Channel.create()

rmOutput
.splitFasta(record: [id: true, text: true])
.choice(identityUnknown, identityKnown) { record -> record.id =~ /#Unknown/ ? 0 : 1 }

repeatmaskerUnknowns = identityUnknown.collectFile() { record -> ['unknown.fasta', record.text] }
repeatmaskerKnowns = identityKnown.collectFile() { record -> ['known.fasta', record.text] }

process transposonBlast {
  container 'robsyme/nf-repeatmasking'
  cpus 4

  input:
  file 'transposases.fasta.gz' from trnaprot
  file 'repeatmodeler_unknowns.fasta' from repeatmaskerUnknowns

  output:
  file 'identified_elements.txt' into identifiedTransposons
  file 'unknown_elements.txt' into unknownElements

  """
zcat transposases.fasta.gz > transposases.fasta
makeblastdb \
 -in transposases.fasta \
 -dbtype prot
blastx \
 -query repeatmodeler_unknowns.fasta \
 -db transposases.fasta \
 -evalue 1e-10 \
 -num_descriptions 10 \
 -num_threads ${task.cpus} \
 -out modelerunknown_blast_results.txt
transposon_blast_parse.pl \
 --blastx modelerunknown_blast_results.txt \
 --modelerunknown repeatmodeler_unknowns.fasta
  """
}

repeatmaskerKnowns
.mix(identifiedTransposons)
.collectFile() { it.text }
.combine(allLTR2)
.set { knownRepeats }

process repeatMaskerKnowns {
  container 'robsyme/nf-repeatmasking'
  publishDir "${params.outdir}/repeatMaskerKnowns", mode: 'copy'

  input:
  file 'reference.fasta' from reference
  set 'knownTransposons.lib', 'allLTRs.lib' from knownRepeats

  output:
  set 'reference.fasta.out', 'reference.fasta.masked' into repeatMaskerKnownsMasked
  file 'reference.fasta.out.gff'

  """
cat *.lib > knownRepeats.fasta
RepeatMasker \
 -lib knownRepeats.fasta \
 -nolow \
 -no_is \
 -dir . \
 -gff \
 -s \
 reference.fasta
  """
}

process octfta {
  container 'robsyme/nf-repeatmasking'

  input:
  file 'reference.fa' from reference
  set 'rm.out', 'rm.masked' from repeatMaskerKnownsMasked

  output:
  file 'summary.tsv' into repeatmaskerSummaryTable

  """
build_dictionary.pl --rm rm.out > ltr.dict
one_code_to_find_them_all.pl --rm rm.out --ltr ltr.dict --fasta reference.fa
echo -e 'Family\\tElement Length\\tFragments\\tCopies\\tSolo_LTR\\tTotal_Bp\\tCover\\tchrname' > summary.tsv
for file in *.copynumber.csv; do
  chrname=`echo \$file | sed -e 's/^rm\\.out_//' -e 's/.copynumber.csv\$//'`
  awk -v chrname=\$chrname 'BEGIN{OFS="\\t"} NR>1 && /^[^#]/ {print(\$0, chrname)}' \$file
done >> summary.tsv
  """
}

process summarise {
  publishDir "${params.outdir}/summarise", mode: 'copy'

  input:
  file 'summary.tsv' from repeatmaskerSummaryTable

  output:
  set 'summary.bycontig.tidy.tsv', 'summary.tidy.tsv'

  """
#!/usr/bin/env Rscript
library(ggplot2)
library(dplyr)
library(tidyr)
library(magrittr)

data <- read.table('summary.tsv', header=TRUE) %>%
        separate(Family, into=c("Family", "Subfamily"), sep="/") %>%
        group_by(chrname, Family, Subfamily) %>%
        summarise(fragment.count = sum(Fragments), length = sum(Total_Bp)) %>%
        unite("Family", Family, Subfamily, sep="/")

write.table(data, file='summary.bycontig.tidy.tsv')

data <- read.table('summary.tsv', header=TRUE) %>%
        separate(Family, into=c("Family", "Subfamily"), sep="/") %>%
        group_by(Family, Subfamily) %>%
        summarise(fragment.count = sum(Fragments), length = sum(Total_Bp)) %>%
        unite("Family", Family, Subfamily, sep="/")

write.table(data, file='summary.tidy.tsv', row.names = FALSE)
  """
}

workflow.onComplete {
	log.info "Pipeline completed at: $workflow.complete"
	log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
