nextflow.enable.dsl = 2

params.raw = 'data/simulated/raw_plate_reader.csv'
params.plate_map = 'data/simulated/plate_map.csv'
params.metadata = 'data/simulated/run_metadata.csv'
params.outdir = 'results/nextflow_demo'

process VALIDATE_INPUTS {
  tag "validate"
  publishDir "${params.outdir}/validation", mode: 'copy'

  input:
  path raw
  path plate_map
  path metadata

  output:
  path "validated_inputs.txt"
  path raw, emit: raw_out
  path plate_map, emit: plate_map_out
  path metadata, emit: metadata_out

  script:
  """
  test -s ${raw}
  test -s ${plate_map}
  test -s ${metadata}
  printf "raw=%s\\nplate_map=%s\\nmetadata=%s\\n" "${raw}" "${plate_map}" "${metadata}" > validated_inputs.txt
  """
}

process RUN_HEKBLUE_ANALYSIS {
  tag "hekblue-analysis"
  publishDir "${params.outdir}", mode: 'copy'
  cpus 1
  memory '2 GB'
  time '1h'

  input:
  path raw
  path plate_map
  path metadata

  output:
  path "analysis_package"

  script:
  """
  mkdir -p analysis_package
  export HEKBLUER_HOME="${projectDir}"
  Rscript "${projectDir}/scripts/run_pipeline_cli.R" \\
    --raw ${raw} \\
    --plate-map ${plate_map} \\
    --metadata ${metadata} \\
    --out analysis_package
  """
}

workflow {
  raw_ch = Channel.fromPath(params.raw, checkIfExists: true)
  plate_map_ch = Channel.fromPath(params.plate_map, checkIfExists: true)
  metadata_ch = Channel.fromPath(params.metadata, checkIfExists: true)

  validated = VALIDATE_INPUTS(raw_ch, plate_map_ch, metadata_ch)
  RUN_HEKBLUE_ANALYSIS(validated.raw_out, validated.plate_map_out, validated.metadata_out)
}
