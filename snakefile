#Steph_wildlife_trial snakefile from 13 June 2022
#trial with Oct 2023 data

configfile: "config_40-919527447_ThBa_opos.yaml"
# print (config['samples'])

rule main:
  # input: expand("results/{sample}/{sample}.clusters.blast.taxonomy.log", sample=config["samples"])
  input: expand("results/{sample}/{sample}.{pcr}.clusters.blast.log", sample=config["samples"], pcr=config["primers"].keys())

rule trim:
 output:
  read1 = "results/{sample}/{sample}.trimmed_read1.fastq",
  read2 = "results/{sample}/{sample}.trimmed_read2.fastq",
  singles = "results/{sample}/{sample}.trimmed_singles.fastq",
  log = "results/{sample}/{sample}.trimming_by_sickle.log"
 input:
  read1 =lambda wildcards: "fastq/{}_1.fastq.gz".format(config['reads'][wildcards.sample]),
  read2 =lambda wildcards: "fastq/{}_2.fastq.gz".format(config['reads'][wildcards.sample]),
 params:
  qual_threshold = config["read_quality_threshold"],
  len_threshold = config["read_trim_length_threshold"],
  folder = "results/{sample}"
 shell:
  r"""sickle pe -f {input.read1} -r {input.read2} -o {output.read1} -p {output.read2} -s {output.singles} -t sanger -q {params.qual_threshold} -l {params.len_threshold} > {output.log}
      fastqc -f fastq -o {params.folder} {input.read1}
      fastqc -f fastq -o {params.folder} {input.read2}
   """

rule overlap:
 output: "results/{sample}/{sample}.overlap_by_flash.log"
 input:
  read1 = "results/{sample}/{sample}.trimmed_read1.fastq",
  read2 = "results/{sample}/{sample}.trimmed_read2.fastq",
 params:
  min_overlap = config["minimum_overlap"],
  max_overlap = config["maximum_overlap"],
  directory = "results/{sample}",
  prefix = "{sample}"
 shell:
  r"""flash -m {params.min_overlap} -M {params.max_overlap} -O -o {params.prefix} -d {params.directory} --threads=1 {input.read1} {input.read2} > {output}
      fastqc -f fastq -o {params.directory} {params.directory}/{params.prefix}.extendedFrags.fastq
   """

rule sort:
 output: "results/{sample}/{sample}.{pcr}.sortPrimers.log"
 input: "results/{sample}/{sample}.overlap_by_flash.log"
 params:
  prefix = "{sample}",
  directory = "results/{sample}",
  pcr = "{pcr}",
  primers = lambda wildcards: config['primers'][wildcards.pcr]
 shell:
  "perl scripts/sortPrimers.pl {params.prefix} {params.directory} {params.pcr} {params.primers} > {output}"

rule blast:
 output: "results/{sample}/{sample}.{pcr}.clusters.blast"
 input: 
  log = "results/{sample}/{sample}.{pcr}.sortPrimers.log"
 params:
  reference = lambda wildcards: config['database'][wildcards.pcr],
  query = "results/{sample}/{sample}.{pcr}.clusters.fasta"
 shell:
  "blastn -db {params.reference} -query {params.query} -outfmt '6 qseqid sseqid pident length qlen qstart qend slen sstart send mismatch gapopen evalue bitscore' -out {output}"
   
rule analyse_blast:
 output: "results/{sample}/{sample}.{pcr}.clusters.blast.log"
 input: "results/{sample}/{sample}.{pcr}.clusters.blast"
 params:
  prefix = "{sample}",
  directory = "results/{sample}",
  taxonomy = config['taxonomy'],
  pcr = "{pcr}"
 shell:
  "perl scripts/reportMapping.pl {params.prefix} {params.directory} {input} {params.taxonomy} {params.pcr} > {output}"





