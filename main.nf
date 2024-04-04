#!/usr/bin/env nextflow
RES_DIR = params.resultsDir

process count_samples {
    label 'bcftools'
    def id = "01_count_samples"
    publishDir "$RES_DIR/$id", mode: 'copy', overwrite: true

    input:
        path(vcf)

    output:
        path "*.samples" 

    script:
    """
    bcftools query -l ${vcf} > ${vcf.getFileName()}.${task.index}.samples
    """
}

process plot_allelefreqs {
    cache false   
    def id = "02_plot_allelefreqs"
    publishDir "$RES_DIR/$id", mode: 'copy', overwrite: true

    input:
        tuple path(vcf), path(notebook)

    output:
        path "*.*" 

    script:
    """    
    papermill ${notebook} \
        ${id}-${vcf}-NB-0.ipynb \
         --report-mode \
         --prepare-only \
         -p vcf_file ${vcf} 
    jupyter nbconvert ${notebook} \
        --to html \
        --stdout > ${id}-${vcf}-NB-B.html
    """
}

process make_report {
    cache false   

    publishDir "$RES_DIR", mode: 'copy', overwrite: true
    
    input:
        tuple val(samplefiles), path(html_report), path(javascript)

    output:
        tuple path("*.html"), path("index.js")
        

    script:
    """
    make_report.py --samples ${samplefiles.join(" ")} --output sample_memberships.html
    """
}


process reflect {
    input:
        val subject

    // output:

    println "Project : $workflow.projectDir"

    exec:
        println "REFLECT"        
        println subject
        
}


process countlines {
    def id = "03_countlines"
    def id2 = "${workflow.projectDir}"
    // id3 { "$task.process" }

    publishDir "$RES_DIR/$id", mode: 'copy', overwrite: true

    input:
        path(datafile)

    output:
        path "counts.csv" 

    println "Project : $workflow.projectDir"
    // println "task.process : "
    
    script:
    """
    wc -l ${datafile} > counts.csv
    echo ${task.process} >> counts.csv
    echo ${id2} >> counts.csv
    echo ${task.process} >> counts.csv
    """
}

workflow {
    if (params.input == null) {
        // exit 1, helpMessage("ERROR: no --input specified")
        VCFMASK = './input/test/o*cohort*.vcf.gz'
        ch_vcfs = Channel.fromPath(VCFMASK)
    } else {
        ch_vcfs = Channel.fromPath(params.input.tokenize(' '))
    }
    
    ch_samples = count_samples(ch_vcfs)

    // ch_template = Channel.fromPath('assets/notebooks/plot_allelefreqs.ipynb')
    ch_template = Channel.fromPath("${workflow.projectDir}/assets/notebooks/plot_allelefreqs.ipynb")//.view()
    
    
    // ch_vcfs.combine(ch_template).view()


    ch_plots = plot_allelefreqs(ch_vcfs.combine(ch_template))

    ch_fancychannel = ch_samples.map{it -> ["sample",it]}.groupTuple().map{it-> [it[1]]}
        .combine([["${workflow.projectDir}/assets/html/report.html",
        "${workflow.projectDir}/assets/js/index.js"]]).view()

    //reflect(ch_samples.toSortedList())
    
    
    // make_report(ch_samples.toSortedList())
     make_report(ch_fancychannel)



    // ch_notebook = Channel.fromPath('./analyses/01_generate_data.ipynb')
    // generate_data(ch_notebook) | countlines | view { "directory contents: $it" }
}
