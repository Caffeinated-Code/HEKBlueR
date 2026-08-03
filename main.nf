nextflow.enable.dsl = 2

params.raw = 'data/simulated/raw_plate_reader.csv'
params.plate_map = 'data/simulated/plate_map.csv'
params.metadata = 'data/simulated/run_metadata.csv'
params.samplesheet = ''
params.outdir = 'results/nextflow_demo'

process VALIDATE_INPUTS {
  tag "validate"
  publishDir "${params.outdir}/validation", mode: 'copy'

  input:
  tuple val(run_id), path(raw), path(plate_map), path(metadata)

  output:
  path "validated_inputs.txt"
  tuple val(run_id), path(raw), path(plate_map), path(metadata), emit: validated_files

  script:
  """
  test -s ${raw}
  test -s ${plate_map}
  test -s ${metadata}
  printf "run_id=%s\\nraw=%s\\nplate_map=%s\\nmetadata=%s\\n" "${run_id}" "${raw}" "${plate_map}" "${metadata}" > validated_inputs.txt
  """
}

process RUN_HEKBLUE_ANALYSIS {
  tag "hekblue-analysis"
  publishDir "${params.outdir}", mode: 'copy'
  cpus 1
  memory '2 GB'
  time '1h'

  input:
  tuple val(run_id), path(raw), path(plate_map), path(metadata)

  output:
  path "${run_id}_analysis_package"

  script:
  """
  mkdir -p ${run_id}_analysis_package
  export HEKBLUER_HOME="${projectDir}"
  Rscript "${projectDir}/scripts/run_pipeline_cli.R" \\
    --raw ${raw} \\
    --plate-map ${plate_map} \\
    --metadata ${metadata} \\
    --out ${run_id}_analysis_package
  """
}

workflow {
  if (params.samplesheet && params.samplesheet.toString().trim()) {
    input_ch = Channel
      .fromPath(params.samplesheet, checkIfExists: true)
      .splitCsv(header: false)
      .map { row -> row.collect { it == null ? null : it.toString().replaceAll('"', '').trim() } }
      .filter { row -> row[0] != 'run_id' }
      .map { row -> tuple(row[0] ?: 'hekblue_run', file(row[1]), file(row[2]), file(row[3])) }
  } else {
    input_ch = Channel.of(tuple('hekblue_run', file(params.raw), file(params.plate_map), file(params.metadata)))
  }

  validated = VALIDATE_INPUTS(input_ch)
  RUN_HEKBLUE_ANALYSIS(validated.validated_files)
}
