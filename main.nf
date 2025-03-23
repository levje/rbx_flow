#!/usr/bin/env nextflow

if(params.help) {
    usage = file("$baseDir/USAGE")
    cpu_count = Runtime.runtime.availableProcessors()

    bindings = ["atlas_directory":"$params.atlas_directory",
                "minimal_vote_ratio":"$params.minimal_vote_ratio",
                "seed":"$params.seed",
                "outlier_alpha":"$params.outlier_alpha",
                "register_processes":"$params.register_processes",
                "rbx_processes":"$params.rbx_processes",
                "single_dataset_size_GB":"$params.single_dataset_size_GB",
                "cpu_count":"$cpu_count"]

    engine = new groovy.text.SimpleTemplateEngine()
    template = engine.createTemplate(usage.text).make(bindings)
    print template.toString()
    return
}

log.info "SCIL RecobundlesX pipeline"
log.info "=========================="
log.info ""
log.info "Start time: $workflow.start"
log.info ""

log.debug "[Command-line]"
log.debug "$workflow.commandLine"
log.debug ""

log.info "[Git Info]"
log.info "$workflow.repository - $workflow.revision [$workflow.commitId]"
log.info ""

log.info "Options"
log.info "======="
log.info ""
log.info "[Atlas]"
log.info "Atlas Directory: $params.atlas_directory"
log.info ""
log.info "[Recobundles options]"
log.info "Minimal Vote Percentage: $params.minimal_vote_ratio"
log.info "Random Seed: $params.seed"
log.info "Outlier Removal Alpha: $params.outlier_alpha"
log.info ""
log.info ""

log.info "Input: $params.input"
root = file(params.input)
/* Watch out, files are ordered alphabetically in channel */
tractogram_for_recognition = Channel
     .fromFilePairs("$root/**/{*tracking*.*,}",
                    size: -1,
                    maxDepth:1) {it.parent.name}

Channel
    .fromPath("$root/**/*fa.nii.gz",
                    maxDepth:1)
    .map{[it.parent.name, it]}
    .into{anat_for_registration;anat_for_reference_bundles}

atlas_directory = Channel.fromPath("$params.atlas_directory/atlas")

Channel.fromPath("$params.atlas_directory/mni_masked.nii.gz")
    .into{atlas_anat;atlas_anat_for_average}
atlas_config = Channel.fromPath("$params.atlas_directory/config_fss_1.json")
centroids_dir = Channel.fromPath("$params.atlas_directory/centroids/",
                                   type: 'dir')
centroids_dir
    .into{atlas_centroids;atlas_centroids_for_average}

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}

anat_for_registration
    .combine(atlas_anat)
    .set{anats_for_registration}
process Register_Anat {
    cpus params.register_processes
    memory '2 GB'

    input:
    set sid, file(native_anat), file(atlas) from anats_for_registration

    output:
    set sid, "${sid}__output0GenericAffine.mat" into transformation_for_recognition
    file "${sid}__outputWarped.nii.gz"
    file "${sid}__native_anat.nii.gz"

    script:
    """
    export ANTS_RANDOM_SEED=1234
    antsRegistrationSyNQuick.sh -d 3 -f ${native_anat} -m ${atlas} -n ${params.register_processes} -o ${sid}__output -t a
    cp ${native_anat} ${sid}__native_anat.nii.gz
    """
}

tractogram_for_recognition
    .join(anat_for_reference_bundles)
    .join(transformation_for_recognition)
    .combine(atlas_config)
    .combine(atlas_directory)
    .set{tractogram_and_transformation}
process Recognize_Bundles {
    cpus params.rbx_processes
    memory { params.single_dataset_size_GB.GB * params.rbx_processes }

    input:
    set sid, file(tractograms), file(refenrence), file(transfo), file(config), file(directory) from tractogram_and_transformation

    output:
    set sid, "*.trk" into bundles_for_cleaning
    file "results.json"
    file "logfile.txt"

    script:
    """
    mkdir tmp/
    scil_recognize_multi_bundles.py ${tractograms} ${config} ${directory}/ ${transfo} --inverse --out_dir tmp/ \
        --log_level DEBUG --minimal_vote_ratio $params.minimal_vote_ratio \
        --seed $params.seed --processes $params.rbx_processes
    mv tmp/* ./
    """
}

bundles_for_cleaning
    .combine(transformation_for_average, by:0)
    .combine(atlas_anat_for_average)
    .set{all_bundles_transfo_for_clean_average}
process Clean_Bundles {
    input:
    set sid, file(bundles), file(transfo), file(atlas) from all_bundles_transfo_for_clean_average

    output:
    set sid, "${sid}__*_cleaned.trk" into bundle_for_count

    script:
    String bundles_list = bundles.join(", ").replace(',', '')
    """
    for bundle in $bundles_list;
        do if [[ \$bundle == *"__"* ]]; then
            pos=\$((\$(echo \$bundle | grep -b -o __ | cut -d: -f1)+2))
            bname=\${bundle:\$pos}
            bname=\$(basename \$bname .trk)
        else
            bname=\$(basename \$bundle .trk)
        fi

        scil_outlier_rejection.py \${bundle} "${sid}__\${bname}_cleaned.trk" \
            --alpha $params.outlier_alpha
            
        if [ -s "${sid}__\${bname}_cleaned.trk" ]; then 
            echo "Bundle \${bname} cleaned."
        else
            echo "After cleaning \${bundle} all streamlines were outliers."
        fi
    done
    """
}

process Count_Streamlines {
    input:
    set sid, file(bundles) from bundle_for_count

    output:
    set sid, file("${sid}__*_count.txt") into bundle_counts

    script:
    String bundles_list = bundles.join(", ").replace(',', '')
    """
    total_count=0
    for bundle in $bundles_list; do
        bname=\$(basename \${bundle} .trk)
        scil_count_streamlines.py \${bundle} --print_count_alone > ${sid}__\${bname}_count.txt
        count=\$(< ${sid}__\${bname}_count.txt)
        total_count=\$((total_count + count))
    done
    echo "Total streamlines counted: \$total_count" > ${sid}__total_count.txt
    """
}
