#!/usr/bin/env nextflow

params.acc = "/data/acc.txt"
List words=[]
file(params.acc).eachLine { line ->
    words.add(line)
}
inputChannel=Channel
.fromList(words)


process gen_ref{
    publishDir "/data/indices", mode: 'symlink'
        input:
            val url from "ftp://ftp.ensemblgenomes.org/pub/plants/release-40/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"
        output:
            path "*.ht2" into indices    

    """
    wget $url 
    gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
    hisat2-build Arabidopsis_thaliana.TAIR10.dna.toplevel.fa TAIR10
    """
}

process download_annotation{
    publishDir "/data/indices", mode: 'symlink'
    input:
        val url from "ftp://ftp.ensemblgenomes.org/pub/plants/release-50/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.50.gff3.gz"
    output:
        path "*.gff3" into annotations    

    """
    wget $url 
    gunzip Arabidopsis_Thaliana.TAIR10.50.gff3.gz 
    """    
}

process download_sra{

    input:
        val acc_id from inputChannel

    output:
        path "${acc_id}.sra"  into srafiles

    """
    prefetch -p --output-directory ./ $acc_id
    mv $acc_id/${acc_id}.sra .
    """
}

process convert{
    input:
        path sra_file from srafiles
    output:
        path "*.fastq" into fastqs

    """
    fasterq-dump -e 8 $sra_file 
    """
    
}

process quality_control{
    input:
        path fastq from fastqs
    output:
        path "./trimmed.fastq" into trimmed

    """
    fastp -i $fastq -q 30  -o ./trimmed.fastq -w 4
    """
}

process mapping{
    input:
        path fastq from trimmed
        path index from indices  
    output:
        path "output.sam" into sams

    """
    hisat2 -q -p 16 -x /data/indices/TAIR10 -S ./output.sam -U ${fastq}
    """
}

process create_bam{
    input:
        path sam from sams
    output:
        path "output.bam" into bams
    
    """
    samtools sort -@ 16 -O bam -o output.bam ${sam}
    samtools index output.bam
    """    
}

process calc_status{
    input:
        path ann from annotations
        path bam from bams 
    output:
        path "output.gtf" into gtfs
    """
    stringtie $bam -o output.gtf -G $ann -A ex.tab
    """              
}



gtfs.view{
    "directory content: $it" 
}
