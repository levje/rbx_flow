def check_required_params(param_names) {
    // Loop through each parameter name and check if it exists in params
    // We need to accumulate errors to report them all at once
    def missing_params = []
    param_names.each { param ->
        if (!params.containsKey(param) || params[param] == false || params[param] == '' || params[param] == null) {
            missing_params << param
        }
    }

    if (missing_params) {
        throw new Exception("Missing required parameters: ${missing_params.join(', ')}")
    }
}

def check_nb_cpus() {
    if(params.processes) {
        if(params.processes > Runtime.runtime.availableProcessors()) {
            throw new RuntimeException("Number of processes higher than available CPUs.")
        }
        else if(params.processes < 1) {
            throw new RuntimeException("When set, number of processes must be >= 1 " +
                                    "and smaller or equal to the number of CPUs.")
        }
    }
}

workflow HANDLE_USAGE_AND_HEADER {
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

    if (params.help) {
        engine = new groovy.text.SimpleTemplateEngine()
        template = engine.createTemplate(usage.text).make(bindings)
        print template.toString()
        System.exit(0)
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

    required_params = ["input", "atlas_directory"]
    check_required_params(required_params)

    log.info "Options"
    log.info "======="
    log.info "[Input] Atlas directory: $params.atlas_directory"
    log.info "[Input] Input: $params.input"
    log.info ""
    log.info "[RBX] Minimal Vote Percentage: $params.minimal_vote_ratio"
    log.info "[RBX] Random Seed: $params.seed"
    log.info "[RBX] Outlier Removal Alpha: $params.outlier_alpha"

    check_nb_cpus()
}
