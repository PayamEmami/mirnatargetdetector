#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/mirnatargetdetector
========================================================================================
 nf-core/mirnatargetdetector Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/mirnatargetdetector
----------------------------------------------------------------------------------------
*/

log.info Headers.nf_core(workflow, params.monochrome_logs)

////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////+
def json_schema = "$projectDir/nextflow_schema.json"
if (params.help) {
    def command = "nextflow run nf-core/mirnatargetdetector --input '*_R{1,2}.fastq.gz' -profile docker"
    log.info NfcoreSchema.params_help(workflow, params, json_schema, command)
    exit 0
}

////////////////////////////////////////////////////
/* --         VALIDATE PARAMETERS              -- */
////////////////////////////////////////////////////+
if (params.validate_params) {
    NfcoreSchema.validateParameters(params, json_schema, log)
}

////////////////////////////////////////////////////
/* --     Collect configuration parameters     -- */
////////////////////////////////////////////////////



// Check AWS batch settings
if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, 'Specify correct --awsqueue and --awsregion parameters on AWSBatch!'
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, 'Outdir not on S3 - specify S3 Bucket to run on AWSBatch!'
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (params.tracedir.startsWith('s3:')) exit 1, 'Specify a local tracedir or run without trace! S3 cannot be used for tracefiles.'
}

// Stage config files
ch_multiqc_config = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs = file("$projectDir/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("$projectDir/docs/images/", checkIfExists: true)

/*
 * Create a channel for input read files
 */


Channel.fromPath(params.input).ifEmpty { exit 1, 'params.input_mirna was empty - no input files supplied' }.into{targets}
Channel.fromPath(params.input_utr).ifEmpty { exit 1, 'params.input_utr was empty - no input files supplied' }.into{reference}




////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////
log.info NfcoreSchema.params_summary_log(workflow, params, json_schema)

// Header log info
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = workflow.runName
// TODO nf-core: Report custom parameters here
summary['Input']            = params.input
summary['UTRs']        = params.input_utr
summary['kmer']        = params.kmer
summary['number_of_chunks_mirna']        = params.number_of_chunks_mirna
summary['number_of_chunks_utr']        = params.number_of_chunks_utr
summary['skip_miranda']        = params.skip_miranda
summary['skip_rnahybrid']        = params.skip_rnahybrid
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if (workflow.profile.contains('awsbatch')) {
    summary['AWS Region']   = params.awsregion
    summary['AWS Queue']    = params.awsqueue
    summary['AWS CLI']      = params.awscli
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Profile Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Profile Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config Profile URL']         = params.config_profile_url
summary['Config Files'] = workflow.configFiles.join(', ')
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
    summary['MultiQC maxsize']   = params.max_multiqc_email_size
}

// Check the hostnames against configured profiles
checkHostname()

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-mirnatargetdetector-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/mirnatargetdetector Workflow Summary'
    section_href: 'https://github.com/nf-core/mirnatargetdetector'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }

/*
 * Parse software version numbers
 */
process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.indexOf('.csv') > 0) filename
                      else null
        }

    output:
    file 'software_versions_mqc.yaml' into ch_software_versions_yaml
    file 'software_versions.csv'

    script:
    // TODO nf-core: Get all tools to print their version number here
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    #fastqc --version > v_fastqc.txt
    #multiqc --version > v_multiqc.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}



/*
 * STEP 1 - Split utrs
 */

process split_utrs{
  tag "$ref"
  label 'process_medium'
  publishDir "${params.outdir}/split_UTRs", mode: params.publish_dir_mode

input:
file ref from reference
output:
file "*.fasta" into shuffle_process

"""
pyfasta split -n $params.number_of_chunks_utr $ref
"""

}


/*
 * STEP 2 - Split mirna
 */

process split_mirna{
  tag "$ref"
  label 'process_medium'
  publishDir "${params.outdir}/split_mirna", mode: params.publish_dir_mode



input:
file ref from targets
output:
file "*.fasta" into shuffle_process_f1

"""
pyfasta split -n $params.number_of_chunks_mirna $ref
"""

}


shuffle_process_fl_target=shuffle_process_f1.flatten()

shuffle_process_fl=shuffle_process.flatten()


/*
 * STEP 3 - shuffle utrs
 */

process shuffle_utrs{
 tag "$ref"
 label 'process_medium'
 publishDir "${params.outdir}/shuffle_utrs", mode: params.publish_dir_mode


input:
file ref from shuffle_process_fl
output:
file "out/${ref.baseName}_shuffled_combined.fasta" into miranda_process, rnahybrid_process

"""
mkdir out
mkdir tmp
fasta-shuffle-letters -kmer $params.kmer -seed 100 $ref tmp/$ref
cat $ref tmp/$ref > out/${ref.baseName}_shuffled_combined.fasta
"""

}



/*
 * STEP 4 - shuffle mirna
 */

process shuffle_mirna{
  tag "$ref"
  label 'process_medium'
  publishDir "${params.outdir}/shuffle_mirna", mode: params.publish_dir_mode


input:
file ref from shuffle_process_fl_target
output:
file "out/${ref.baseName}_shuffled_combined.fasta" into miranda_process_targets, rnahybrid_process_targets

"""
mkdir out
mkdir tmp
fasta-shuffle-letters -kmer $params.kmer -seed 100 $ref tmp/$ref
cat $ref tmp/$ref > out/${ref.baseName}_shuffled_combined.fasta
"""
}


/*
 * STEP 4 - perform miranda
 */


process miranda{
tag "$query $ref"
label 'process_medium'
publishDir "${params.outdir}/miranda", mode: params.publish_dir_mode

when:
params.skip_miranda == false

input:
each file(query) from miranda_process_targets
each file(ref) from miranda_process

output:
file "out/${ref.baseName}_${query.baseName}_results.txt" into combine

"""
mkdir out
miranda $query $ref  -quiet -strict -out out/${ref.baseName}_${query.baseName}_results.txt
"""

}


/*
 * STEP 5 - perform rnahybrid
 */


process rnahybrid{
 publishDir 'results/rnahybrid'
tag "$query $ref"
label 'process_medium'
publishDir "${params.outdir}/rnahybrid", mode: params.publish_dir_mode

when:
params.skip_rnahybrid == false

input:
each file(query) from rnahybrid_process_targets
each file(ref) from rnahybrid_process

output:
file "out/${ref.baseName}_${query.baseName}_results_rnahybrid.txt" into combine_hybrid

"""
mkdir out
RNAhybrid -s 3utr_fly -m 9000 -n 100 -q $query -t $ref > out/${ref.baseName}_${query.baseName}_results_rnahybrid.txt
"""

}


/*
 * STEP 6 - merge miranda
 */


process miranda_merge{
  label 'process_medium'
  publishDir "${params.outdir}/miranda_merge", mode: params.publish_dir_mode

  when:
  params.skip_miranda == false

input:
file query from combine.collect()


output:
file "miranda_merged_results.txt" into combineOut

"""
cat *.txt > miranda_merged_results.txt

"""

}


process rnahybrid_merge{
  label 'process_medium'
  publishDir "${params.outdir}/rnahybrid_merge", mode: params.publish_dir_mode

  when:
  params.skip_rnahybrid == false

input:
file query from combine_hybrid.collect()

output:
file "rnahybrid_merged_results.txt" into combineOutRnahybrid

"""
cat *.txt > rnahybrid_merged_results.txt

"""
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/mirnatargetdetector] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[nf-core/mirnatargetdetector] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp



    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$projectDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$projectDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, projectDir: "$projectDir", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
    def sf = new File("$projectDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[nf-core/mirnatargetdetector] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            def mail_cmd = [ 'mail', '-s', subject, '--content-type=text/html', email_address ]
            if ( mqc_report.size() <= params.max_multiqc_email_size.toBytes() ) {
              mail_cmd += [ '-A', mqc_report ]
            }
            mail_cmd.execute() << email_html
            log.info "[nf-core/mirnatargetdetector] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[nf-core/mirnatargetdetector]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[nf-core/mirnatargetdetector]${c_red} Pipeline completed with errors${c_reset}-"
    }

}

workflow.onError {
    // Print unexpected parameters - easiest is to just rerun validation
    NfcoreSchema.validateParameters(params, json_schema, log)
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = 'hostname'.execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "${c_red}====================================================${c_reset}\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "${c_red}====================================================${c_reset}\n"
                }
            }
        }
    }
}
