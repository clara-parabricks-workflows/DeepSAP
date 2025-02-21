# Benchmarking RNA-seq Aligners using Baruzzo and SimBA Datasets

## Introduction
This document describes our Nextflow pipeline designed for benchmarking RNA-seq aligners using the Baruzzo and SimBA datasets. 

The benchmarked RNA-seq aligners in this pipeline include:
- DeepSAP with GSNAP TGGA v24.06.2024 and DNABERT MS150
- Illumina’s DRAGEN™ v4.0.3
- novoSplice v0.8.4 (Berakdar et al. 2019)
- STAR v2.7.10a (Dobin et al. 2013)
- HISAT2 v2.2.1 (Kim et al. 2019)
- Subjunc v2.0.1 (Liao et al. 2013)

Results from these aligners are collected in a standardized `results.tsv` file and subsequently visualized using an R Shinyapp [Dashboard](https://rna-seqbenchmark.shinyapps.io/DeepSAP-Dashboard).

## Datasets
### Baruzzo Dataset
The Baruzzo dataset originates from the study *"Simulation-based comprehensive benchmarking of RNA-seq aligners"*, from Human and Malaria with T1, T2, and T3 levels of complexities.

- **Article:** [Simulation-based benchmarking of RNA-seq aligners](https://www.nature.com/articles/nmeth.4106#Sec9)
- **BEER Simulator:** [GitHub Repository](https://github.com/itmat/beers_simulator)
- **Baruzzo datasets:** [Homepage](http://bioinf.itmat.upenn.edu/BEERS/bp1/datasets.html)

### SimBA Dataset
The SimBA dataset was generated using the SimCT simulator from *"A methodology and tools for evaluating the performance of RNA-Seq bioinformatic pipelines"* to mimic the complexities of the Baruzzo datasets with the same genome and annotation files.

- **Article:** [SimBA: A methodology and tools for evaluating the performance of RNA-Seq bioinformatic pipelines](https://pubmed.ncbi.nlm.nih.gov/28969586/)
- **SimCT simulator:** [GitHub Repository](https://github.com/jaudoux/simct)

#### Human Simulations
- **T1 Complexity:**
  ```bash
    simCT --vcf-file /data/references/GRCh37_Baruzzo/00-common_all.vcf.gz \
            --annotations /data/references/GRCh37_Baruzzo/simulator_config_geneinfo_hg19_GTF \
            --genome-dir /data/references/GRCh37_Baruzzo/GRCh37.primary_assembly.chr.fasta/ \
            --output-dir /data/simCT/BarruzoSimBA/SimBA_human_t1/ \
            --disable-error-encoding --substitution-rate 0.001 --insertion-rate 0.0001 --deletion-rate 0.0001 \
            --vcf-ratio 0.95 --nb-fusions 100 --fragment-length 250 --fragment-sd 50 --reads-length 100 \
            --nb-molecules 1000000 --nb-reads 20000000
    ```
- **T2 Complexity:**
  ```bash
    simCT --vcf-file /data/references/GRCh37_Baruzzo/00-common_all.vcf.gz \
        --annotations /data/references/GRCh37_Baruzzo/simulator_config_geneinfo_hg19_GTF \
        --genome-dir /data/references/GRCh37_Baruzzo/GRCh37.primary_assembly.chr.fasta/ \
        --output-dir /data/simCT/BarruzoSimBA/SimBA_human_t2/ \
        --disable-error-encoding --substitution-rate 0.005 --insertion-rate 0.0020 --deletion-rate 0.0020 \
        --vcf-ratio 0.95 --nb-fusions 100 --fragment-length 250 --fragment-sd 50 --reads-length 100 \
        --nb-molecules 1000000 --nb-reads 20000000
    ```

- **T3 Complexity:**
  ```bash
    simCT --vcf-file /data/references/GRCh37_Baruzzo/00-common_all.vcf.gz \
        --annotations /data/references/GRCh37_Baruzzo/simulator_config_geneinfo_hg19_GTF \
        --genome-dir /data/references/GRCh37_Baruzzo/GRCh37.primary_assembly.chr.fasta/ \
        --output-dir /data/simCT/BarruzoSimBA/SimBA_human_t3/ \
        --disable-error-encoding --substitution-rate 0.030 --insertion-rate 0.0050 --deletion-rate 0.0050 \
        --vcf-ratio 0.95 --nb-fusions 100 --fragment-length 250 --fragment-sd 50 --reads-length 100 \
        --nb-molecules 1000000 --nb-reads 20000000
    ```


#### Malaria Simulations
- **T1 Complexity:**
  ```bash
    simCT --annotations /data/references/PFAL_Baruzzo/simulator_config_geneinfo_pfal_GTF \
        --genome-dir /data/references/PFAL_Baruzzo/genome_sequence_pfal.chr.fasta \
        --output-dir /data/simCT/BarruzoSimBA/malaria_t1/ \
        --disable-error-encoding --substitution-rate 0.001 --insertion-rate 0.0001 --deletion-rate 0.0001 \
        --nb-fusions 100 --fragment-length 250 --fragment-sd 50 --reads-length 100 \
        --nb-molecules 1000000 --nb-reads 20000000
    ```
- **T2 Complexity:**
  ```bash
    simCT --annotations /data/references/PFAL_Baruzzo/simulator_config_geneinfo_pfal_GTF \
        --genome-dir /data/references/PFAL_Baruzzo/genome_sequence_pfal.chr.fasta \
        --output-dir /data/simCT/BarruzoSimBA/malaria_t1/ \
        --disable-error-encoding --substitution-rate 0.005 --insertion-rate 0.0020 --deletion-rate 0.0020 \
        --nb-fusions 100 --fragment-length 250 --fragment-sd 50 --reads-length 100 \
        --nb-molecules 1000000 --nb-reads 20000000
    ```

- **T3 Complexity:**
  ```bash
    simCT --annotations /data/references/PFAL_Baruzzo/simulator_config_geneinfo_pfal_GTF \
        --genome-dir /data/references/PFAL_Baruzzo/genome_sequence_pfal.chr.fasta \
        --output-dir /data/simCT/BarruzoSimBA/malaria_t1/ \
        --disable-error-encoding --substitution-rate 0.030 --insertion-rate 0.0050 --deletion-rate 0.0050 \
        --nb-fusions 100 --fragment-length 250 --fragment-sd 50 --reads-length 100 \
        --nb-molecules 1000000 --nb-reads 20000000
    ```

**Note:** The VCF file used for human simulations is obtained from [NCBI](https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/00-common_all.vcf.gz)

## Reference Genomes and Annotations
### Benchmarking
The human reference genome is  [GRCh37](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz), however, the Human annotation and Malaria refernce genome and annotation, can be found on the Baruzzo [homepage](http://bioinf.itmat.upenn.edu/BEERS/bp1/datasets.html) 

(http://bp1.s3.amazonaws.com/human.tar.bz2) and Malaria reference [genome and annaotation](http://bp1.s3.amazonaws.com/malaria.tar.bz2).

### Fine-tuning DNABERT MS150 model

#### Human (Homo sapiens RefSeq)
- **Genome:** [GRCh38.p14 FASTA](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/reference/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz)
- **Annotation:** [GRCh38.p14 GTF](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/reference/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz)

#### Human (Homo sapiens GENCODE)
- **Genome:** [GENCODE v44 FASTA](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz)
- **Annotation:** [GENCODE v44 GTF](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz)

#### Zebrafish (Danio rerio)
- **Genome:** [GRCz11 FASTA](https://ftp.ensembl.org/pub/release-110/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.toplevel.fa.gz)
- **Annotation:** [GRCz11 GTF](https://ftp.ensembl.org/pub/release-110/gtf/danio_rerio/Danio_rerio.GRCz11.110.gtf.gz)

#### Mouse (Mus musculus)
- **Genome:** [GRCm39 FASTA](https://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.toplevel.fa.gz)
- **Annotation:** [GRCm39 GTF](https://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz)

#### Fruitfly (Drosophila melanogaster)
- **Genome:** [BDGP6.46 FASTA](https://ftp.ensembl.org/pub/release-110/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa.gz)
- **Annotation:** [BDGP6.46 GTF](https://ftp.ensembl.org/pub/release-110/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.46.110.gtf.gz)

#### C. elegans (Caenorhabditis elegans)
- **Genome:** [WBcel235 FASTA](https://ftp.ensembl.org/pub/release-110/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz)
- **Annotation:** [WBcel235 GTF](https://ftp.ensembl.org/pub/release-110/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.110.gtf.gz)

#### Rat (Rattus norvegicus)
- **Genome:** [mRatBN7.2 FASTA](https://ftp.ensembl.org/pub/release-110/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz)
- **Annotation:** [mRatBN7.2 GTF](https://ftp.ensembl.org/pub/release-110/gtf/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.110.gtf.gz)

#### Malaria (Plasmodium species)
##### Plasmodium falciparum
- **Genome:** [ASM276v2 FASTA](http://ftp.ensemblgenomes.org/pub/protists/release-57/fasta/plasmodium_falciparum/dna/Plasmodium_falciparum.ASM276v2.dna.toplevel.fa.gz)
- **Annotation:** [ASM276v2 GTF](http://ftp.ensemblgenomes.org/pub/protists/release-57/gtf/plasmodium_falciparum/Plasmodium_falciparum.ASM276v2.57.gtf.gz)

##### Plasmodium vivax
- **Genome:** [ASM241v2 FASTA](http://ftp.ensemblgenomes.org/pub/protists/release-57/fasta/plasmodium_vivax/dna/Plasmodium_vivax.ASM241v2.dna.toplevel.fa.gz)
- **Annotation:** [ASM241v2 GTF](http://ftp.ensemblgenomes.org/pub/protists/release-57/gtf/plasmodium_vivax/Plasmodium_vivax.ASM241v2.57.gtf.gz)

##### Plasmodium berghei
- **Genome:** [PBANKA01 FASTA](http://ftp.ensemblgenomes.org/pub/protists/release-57/fasta/plasmodium_berghei/dna/Plasmodium_berghei.PBANKA01.dna.toplevel.fa.gz)
- **Annotation:** [PBANKA01 GTF](http://ftp.ensemblgenomes.org/pub/protists/release-57/gtf/plasmodium_berghei/Plasmodium_berghei.PBANKA01.57.gtf.gz)

##### Plasmodium knowlesi
- **Genome:** [ASM635v1 FASTA](http://ftp.ensemblgenomes.org/pub/protists/release-57/fasta/plasmodium_knowlesi/dna/Plasmodium_knowlesi.ASM635v1.dna.toplevel.fa.gz)
- **Annotation:** [ASM635v1 GTF](http://ftp.ensemblgenomes.org/pub/protists/release-57/gtf/plasmodium_knowlesi/Plasmodium_knowlesi.ASM635v1.57.gtf.gz)

## Alignment runs
### DeepSAP
```bash
    run_id=DeepSAP__${version}__${description}__${sample_id}

    docker run --gpus all --ipc=host --ulimit memlock=-1 --ulimit stack=67108864 --rm \
        --volume  \${PWD}:/workdir \
        --volume  \${PWD}:/outputdir \
        ${container} \
        --out /outputdir/ \
        --prefix \${run_id}.sorted \
        --sam /workdir/${alignment_sam} \
        --fasta /workdir/${genome_FASTA} \
        --gtf /workdir/${genome_GTF} &> \${run_id}.log.txt
```

### DRAGEN
#### Indexing
```bash
    options="--num-threads ${params.nThreads}"
    ${DRAGEN} --build-hash-table true --ht-build-rna-hashtable true ${options} --output-directory  dragen_${version}_${description}.idx/ --ht-reference ${genome_FASTA} 
```

#### Aligning
```bash
    options= "--num-threads ${params.nThreads}"
    run_id="DRAGEN__${version}__${description}__${sample_id}"

    ${DRAGEN} \
    -r ${dragen_index} --output-directory ./ --output-file-prefix \${run_id} --output-format BAM \
    -1 ${reads_forward} -2 ${reads_reverse} --enable-rna true --RGID baruzzo__${description} --RGSM ${sample_id} --annotation-file ${genome_GTF} --enable-rna-gene-fusion true ${options}

    ${params.condaEnv}/bin/samtools sort -n --threads ${params.nThreads} \${run_id}.bam > \${run_id}.sorted.bam
```
### novoSplice 
```bash
    options= "--nThreads ${params.nThreads}"
    run_id="novoSplice__${version}__${description}__${sample_id}"

    ${novoSplice} ${options} \
    --ignoreErrors --fasta ${genome_FASTA} --gtf ${genome_GTF} \
    -1 ${reads_forward} -2 ${reads_reverse} -o ./ 
    
    ${params.condaEnv}/bin/samtools view -Sb --threads ${params.nThreads} alignment.sam | ${params.condaEnv}/bin/samtools sort -n --threads ${params.nThreads} > \${run_id}.sorted.bam
    rm -rf alignment.sam 
```

### STAR 
#### Indexing
```bash
    options="--runThreadN ${params.nThreads} --sjdbOverhang 99"
    ${STAR} ${options} \
    --runMode genomeGenerate \
    --genomeDir star_${version}_${description}.idx/ \
    --genomeFastaFiles ${genome_FASTA} \
    --sjdbGTFfile ${genome_GTF}
```
#### Aligning
```bash
    options="--runThreadN ${params.nThreads} --twopassMode Basic --outSAMunmapped Within"
    run_id="STAR__${version}__${description}__${sample_id}"

    ${STAR} \
    --runMode alignReads \
    --genomeDir ${star_index} \
    --readFilesIn ${reads_forward} ${reads_reverse} \
    --outSAMtype SAM \
    --outFileNamePrefix ./ \
    ${options}

    ${params.condaEnv}/bin/samtools view -Sb --threads ${params.nThreads} Aligned.out.sam | ${params.condaEnv}/bin/samtools sort -n --threads ${params.nThreads} > \${run_id}.sorted.bam
```
### HISAT2
#### Indexing
```bash
    options="-p ${params.nThreads}"
    hisat2_extract_splice_sites.py ${genome_GTF} > \${index_id}/genome.ss
    bin/hisat2_extract_exons.py ${genome_GTF} > \${index_id}/genome.exon

    ${HISAT2_build} ${options} \
    --ss \${index_id}/genome.ss \
    --exon \${index_id}/genome.exon \
    ${genome_FASTA} \${index_id}/index
```
#### Aligning
```bash
    options="-p ${params.nThreads}"
    fasta_flag=""

    if [[ ${reads_forward} == *.fasta ]]
    then
        fasta_flag="-f"
    fi 
    
    run_id="HISAT2__${version}__${description}__${sample_id}"

    ${HISAT2} \
    \${fasta_flag} ${options} -x ${hisat2_index_prefix} -1 ${reads_forward} -2 ${reads_reverse} > \${run_id}.sam 

    ${params.condaEnv}/bin/samtools view --threads ${params.nThreads} -Sb \${run_id}.sam | ${params.condaEnv}/bin/samtools sort -n --threads ${params.nThreads} > \${run_id}.sorted.bam
```

### Subjunc
#### Indexing
```bash
    options="-B"
    ${Subread_build} ${options} -o subjunc_${version}_${description}.idx/index ${genome_FASTA}
```
#### Aligning
```bash
    options="-a ${genome_GTF} --SAMoutput -T ${params.nThreads} --allJunctions"
    run_id="Subjunc__${version}__${description}__${sample_id}"

    ${Subjunc} \
    ${options} -i ${subjunc_index_prefix} -r ${reads_forward} -R ${reads_reverse} -o \${run_id}.sam 

    ${params.condaEnv}/bin/samtools view -Sb --threads ${params.nThreads} \${run_id}.sam  | ${params.condaEnv}/bin/samtools sort -n --threads ${params.nThreads} > \${run_id}.sorted.bam
```

## Benchmarking
### Baruzzo script
- **Benchmarking Scripts:** [GitHub Repository](https://github.com/khayer/aligner_benchmark)
```bash
    samtools view -h ${sample_bam} > alignment.sam
    
    ruby ${params.fixSAM} alignment.sam > alignment.fixed.sam

    ruby ${params.compareToTruth} ${sample_id}.cig alignment.fixed.sam > ${run_id}.bench
```
### SimBA BenchCT script
- **BenchCT:** [GitHub Repository](https://github.com/jaudoux/benchct)
```bash
    ${params.condaEnv}/bin/samtools view -Sbq 5 ${sample_bam} > alignment.unique.bam

    ${params.condaEnv}/bin/python ${params.writeYML} --dataset_info ${sample_id} --output ./ --run_id ${run_id} --path_to_BAM ${sample_bam} --path_to_uniquq_BAM alignment.unique.bam

    ${params.benchCT} bench_all.yml    &> ${run_id}__all_stats.log
    ${params.benchCT} bench_unique.yml &> ${run_id}__unique_stats.log 
```

in which writeYML script is python script that has writeYAML_unique function:
```python
    def writeYAML_unique(_dataset_info, _output, _run_id, _path2SAM):
        buffer = "---\n"
        buffer += "checker:\n"
        buffer += "    files:\n"
        buffer += "        infos:     " + _dataset_info + ".info.txt\n"
        buffer += "        mutations: " + _dataset_info + ".mutations.vcf\n"
        buffer += "        splices:   " + _dataset_info + ".splices.modified.bed\n"
        buffer += "        chimeras:  " + _dataset_info + ".chimeras.tsv\n"
        buffer += "    thresholds:\n"
        buffer += "        MAPPING:   5\n"
        buffer += "        SNP:       5\n"
        buffer += "        INSERTION: 5\n"
        buffer += "        DELETION:  5\n"
        buffer += "        CHIMERA:   20\n"
        buffer += "        ERROR:     5\n"
        buffer += "        SPLICE:    5\n"

        buffer += "softwares:\n"
        buffer += "    - name: " + _run_id + "\n"
        buffer += "      files:\n"
        buffer += "        - name: " + _path2SAM + "\n"
        buffer += "          false_positives: " + _output + _run_id + "__unique_fp\n"
        buffer += "          true_positives: " + _output + _run_id + "__unique_tp\n"
        
        buffer += "          type: SAM::Crac\n"
        buffer += "          check:\n"
        buffer += "            - mapping\n"
        buffer += "            - splice\n"
        buffer += "            - chimera\n"
        buffer += "            - error\n"

        buffer += "output:\n"
        buffer += "    statistics:\n"
        buffer += "      - sensitivity\n"
        buffer += "      - accuracy\n"
        buffer += "      - true-positives\n"
        buffer += "      - false-positives\n"
        buffer += "      - nb-elements\n"
        buffer += "      - false-negatives\n"
        buffer += "    nb_decimals: 4\n"
        
        with open(_output + "/bench_unique.yml", "w") as file_parm:
                file_parm.write(buffer)
```

## Visualization
We hosted the benchmarking results on a ShinyApp Dashboard: [DeepSAP Dashboard](https://rna-seqbenchmark.shinyapps.io/DeepSAP-Dashboard).
