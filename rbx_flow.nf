
process Register_Anat {
    cpus params.register_processes
    memory '2 GB'

    input:
    tuple val(sid), path(native_anat), path(atlas)

    output:
    tuple val(sid), path("${sid}__output0GenericAffine.mat"), emit: transformation
    path "${sid}__outputWarped.nii.gz"
    path "${sid}__native_anat.nii.gz"

    script:
    """
    export ANTS_RANDOM_SEED=1234
    antsRegistrationSyNQuick.sh -d 3 -f ${native_anat} -m ${atlas} -n ${params.register_processes} -o ${sid}__output -t a
    cp ${native_anat} ${sid}__native_anat.nii.gz
    """
}

process Recognize_Bundles {
    cpus params.rbx_processes
    memory { params.single_dataset_size_GB.GB * params.rbx_processes }

    input:
    tuple val(sid), path(tractograms), path(reference), path(transfo), path(config), path(directory)

    output:
    tuple val(sid), path ("*.trk"), path("results.json"), emit: bundles_for_cleaning
    path "logfile.txt"

    script:
    """
    mkdir tmp/

    # Support for both old and new scripts.
    if command -v scil_bundle_reject_outliers.py >/dev/null 2>&1; then
        SCRIPT="scil_tractogram_segment_with_bundleseg.py"
        VERBOSE_FLAG="-v"
    else
        SCRIPT="scil_recognize_multi_bundles.py"
        VERBOSE_FLAG="--log_level"
    fi

    \${SCRIPT} ${tractograms} ${config} ${directory}/ ${transfo} --inverse --out_dir tmp/ \
        \${VERBOSE_FLAG} DEBUG --minimal_vote_ratio $params.minimal_vote_ratio \
        --seed $params.seed --processes $params.rbx_processes
    mv tmp/* ./
    """
}

process Clean_Bundles {
    cpus 1
    memory '2 GB'

    input:
    tuple val(sid), path(bundles), path(results), path(transfo), path(atlas)
    val(outlier_alpha)

    output:
    tuple val(sid), path("${sid}__*_cleaned.trk"), emit: cleaned_bundles
    path "${sid}__results_indices.json"

    script:
    String bundles_list = bundles.join(", ").replace(',', '')
    """
    # Check if the newer outlier rejection script is available.
    if command -v scil_bundle_reject_outliers.py >/dev/null 2>&1; then
        OUTLIER_SCRIPT="scil_bundle_reject_outliers.py"
    else
        OUTLIER_SCRIPT="scil_outlier_rejection.py"
    fi
    echo "Cleaning bundle with \${OUTLIER_SCRIPT}"

    cp ${results} ${sid}__results_indices.json

    for bundle in $bundles_list;
        do if [[ \$bundle == *"__"* ]]; then
            pos=\$((\$(echo \$bundle | grep -b -o __ | cut -d: -f1)+2))
            bname=\${bundle:\$pos}
            bname=\$(basename \$bname .trk)
        else
            bname=\$(basename \$bundle .trk)
        fi

        \${OUTLIER_SCRIPT} \${bundle} "${sid}__\${bname}_cleaned.trk" \
            --alpha $outlier_alpha --results_json ${sid}__results_indices.json --bname \${bname}
            
        if [ -s "${sid}__\${bname}_cleaned.trk" ]; then 
            echo "Bundle \${bname} cleaned."
        else
            echo "After cleaning \${bundle} all streamlines were outliers."
        fi
    done
    """
}

workflow RBX_CORE {
    take:
    anat_channel
    atlas_directory
    atlas_config
    atlas_anat
    input_tractograms
    outlier_alpha

    main:
    anat_for_registration = anat_channel.combine(atlas_anat)
    Register_Anat(anat_for_registration)
    
    tractogram_and_transformation = input_tractograms.join(anat_channel)
        .join(Register_Anat.out.transformation)
        .combine(atlas_config)
        .combine(atlas_directory)
        
    Recognize_Bundles(tractogram_and_transformation)

    all_bundles_transfo_for_clean_average = Recognize_Bundles.out.bundles_for_cleaning
        .combine(Register_Anat.out.transformation, by:0)
        .combine(atlas_anat)

    Clean_Bundles(all_bundles_transfo_for_clean_average, outlier_alpha)
}