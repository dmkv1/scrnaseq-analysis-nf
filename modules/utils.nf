process MATRIX_TO_SCE {
    tag "${sample_id}"

    input:
    tuple val(sample_id), val(patient_id), val(label), val(replicate), path(fbmtx)

    output:
    tuple val(sample_id), path("${sample_id}.sce"), emit: sce

    script:
    """
    matrix_to_sce.R --fbmtx_path '${fbmtx}' \
        --sample_id '${sample_id}' \
        --patient_id '${patient_id}' \
        --label '${label}' \
        --replicate '${replicate}'
    """
}