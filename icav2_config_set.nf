#!/usr/bin/env nextflow
nextflow.enable.dsl=2
params.help = false

/*
 * Configure ICAv2 from the current environment where `icav2` is on PATH while
 * storing the config in a host-mounted credentials directory.  If you launch
 * this workflow inside the ICAv2 CLI container, the same mounted directory is
 * visible on both the host and inside the container.
 *
 * Required inputs:
 *   --server_url   ICA server URL (for example, https://ica.illumina.com)
 *   --x_api_key    ICA API key
 *
 * Optional inputs:
 *   --colormode    CLI color mode (default: none)
 *   --table_format CLI output table format (default: table)
 *   --icav2_container   Path to ICAv2 CLI .sif (default: bin/icav2_cli_v2.43.0.sif)
 *   --icav2_creds_dir   Host credentials directory (default: ~/.icav2)
 */

params.server_url      = null
params.x_api_key       = null
params.colormode       = 'none'
params.table_format    = 'table'
params.icav2_container = "$projectDir/bin/icav2_cli_v2.43.0.sif"
params.icav2_creds_dir = System.getProperty('user.home') + '/.icav2'

workflow {
    def effective_container = params.icav2_container == '/path/to/icav2_cli.sif' \
        ? "$projectDir/bin/icav2_cli_v2.43.0.sif" \
        : params.icav2_container

    if (!params.server_url) {
        throw new IllegalArgumentException('Missing required parameter: --server_url')
    }
    if (!params.x_api_key) {
        throw new IllegalArgumentException('Missing required parameter: --x_api_key')
    }
    if (!effective_container) {
        throw new IllegalArgumentException('Missing required parameter: --icav2_container')
    }
    def container_file = file(effective_container)
    if (!container_file.exists()) {
        throw new IllegalArgumentException("ICAv2 container does not exist: ${effective_container}")
    }
    if (!container_file.name.toLowerCase().endsWith('.sif')) {
        throw new IllegalArgumentException("ICAv2 container must be a .sif file: ${effective_container}")
    }

    ICAV2_CONFIG_SET(
        Channel.value(params.server_url),
        Channel.value(params.x_api_key),
        Channel.value(params.colormode),
        Channel.value(params.table_format),
        Channel.value(effective_container),
        Channel.value(params.icav2_creds_dir)
    )
}

process ICAV2_CONFIG_SET {
    cpus 1
    errorStrategy 'terminate'
    maxRetries 1

    input:
    val server_url
    val x_api_key
    val colormode
    val table_format
    val icav2_container
    val icav2_creds_dir

    script:
    """
    set -euo pipefail

        if ! command -v apptainer >/dev/null 2>&1; then
            echo "apptainer not found in PATH." >&2
            echo "On this cluster, run this workflow from a worker shell (e.g. srun --pty bash)." >&2
            exit 127
        fi
    mkdir -p "${icav2_creds_dir}"

    printf '%s\n%s\n%s\n%s\n' \
      '${server_url}' \
      '${x_api_key}' \
      '${colormode}' \
      '${table_format}' | \
      apptainer exec -B "${icav2_creds_dir}:/root/.icav2" "${icav2_container}" icav2 config set

    echo
    echo "ICAv2 config stored in container path: /root/.icav2"
    echo "Host path bound to this location: ${icav2_creds_dir}"
    apptainer exec -B "${icav2_creds_dir}:/root/.icav2" "${icav2_container}" icav2 config get
    """

    stub:
    """
    set -euo pipefail
    mkdir -p "${icav2_creds_dir}"
    echo "Stub run: would run '${icav2_container}' and configure ICAv2 in ${icav2_creds_dir}"
    """
}