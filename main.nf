#!/usr/bin/env nextflow

// =============================================================================
// mRNA Design Pipeline - Simplified Workflow
// =============================================================================

nextflow.enable.dsl = 2

// Parameters
params {
    input = "${projectDir}/test_input.fasta"
    sample_seqs = 50
}

log.info "Starting mRNA Design Pipeline..."

// Get input file - use absolute path
def inputFile = params.input ? file(params.input) : file("test_input.fasta")
log.info "Input: ${inputFile}"

workflow {
    // Step 1: Parse input
    log.info "Step 1: Parsing input..."
    aa_seq_file = PARSE_INPUT(inputFile)
    
    // Step 2: Generate candidates
    log.info "Step 2: Generating candidate sequences..."
    candidates_file = REVERSE_TRANSLATE(aa_seq_file, params.sample_seqs)
    
    // Step 3: Evaluate all
    log.info "Step 3: Evaluating candidates..."
    all_scores = EVALUATE_ALL(candidates_file)
    
    // Step 4: Rank and generate report
    log.info "Step 4: Ranking and generating report..."
    report_file = RANK_AND_REPORT(all_scores)
    
    log.info "Complete! Report: ${report_file}"
}

process PARSE_INPUT {
    input:
        path input_file
    
    output:
        path "parsed.fasta"
    
    """
    python3 ${projectDir}/bin/parse_input.py ${input_file} > parsed.fasta
    """
}

process REVERSE_TRANSLATE {
    input:
        path aa_seq
        val num_seqs
    
    output:
        path "candidates.fasta"
    
    """
    python3 ${projectDir}/bin/reverse_translate.py < ${aa_seq} ${num_seqs} > candidates.fasta
    """
}

process EVALUATE_ALL {
    input:
        path candidates
    
    output:
        path "all_scores.json"
    
    """
    python3 ${projectDir}/bin/evaluate_all.py < ${candidates} > all_scores.json
    """
}

process RANK_AND_REPORT {
    input:
        path scores_file
    
    output:
        path "report.html"
    
    """
    python3 ${projectDir}/bin/rank_and_report.py ${scores_file} > report.html
    """
}
