#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { HANDLE_USAGE_AND_HEADER } from './inputs.nf'
include { RBX_CORE } from './rbx_flow.nf'

workflow {
    // Print usage and header information.
    // This will also check if the required parameters are set.
    HANDLE_USAGE_AND_HEADER()

    // Directory holding the input data.
    root = file(params.input)

    anat_channel = Channel.fromPath("$root/**/*fa.nii.gz", maxDepth:1)
        .map{[it.parent.name, it]}

    input_tractograms = Channel.fromFilePairs("$root/**/{*tracking*.*,}", size: -1, maxDepth:1) {it.parent.name}

    atlas_directory = Channel.fromPath("$params.atlas_directory/atlas")
    atlas_config = Channel.fromPath("$params.atlas_directory/config_fss_1.json")
    atlas_anat = Channel.fromPath("$params.atlas_directory/mni_masked.nii.gz")

    Channel.fromPath("$params.atlas_directory/config_fss_1.json").set{atlas_config}

    RBX_CORE(anat_channel, atlas_directory, atlas_config, atlas_anat, input_tractograms, params.outlier_alpha)
}
