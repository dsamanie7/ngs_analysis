#!/bin/bash

run_star_alignment() {
    # Function to generate a SLURM job for STAR alignment for a single sample
    # Arguments:
    #   --sample_name: Sample name (e.g., "Sample1")
    #   --input_dir: Directory containing input fastq files (default: predefined path)
    #   --genome_dir: STAR genome directory (default: predefined path)
    #   --output_dir: Directory for output files (default: predefined path)

    # Default values for directories
    local GENOME_DIR="/projects/p31767/Bioinformatics/mouse/star/star_50/"
    local INPUT_DIR="/projects/b1042/Shukla_lab/Daniela/projects/ngs_default_analysis/input_files"
    local OUTPUT_DIR="/projects/b1042/Shukla_lab/Daniela/projects/ngs_default_analysis/output_analysis"
    local SAMPLE_NAME=""
    local PAIRED_END=true
    local SUBMIT_JOB=false

    # Parse arguments
    while [[ "$#" -gt 0 ]]; do
        case $1 in
            --sample_name)
                SAMPLE_NAME="$2"
                shift 2
                ;;
            --input_dir)
                INPUT_DIR="$2"
                shift 2
                ;;
            --genome_dir)
                GENOME_DIR="$2"
                shift 2
                ;;
            --output_dir)
                OUTPUT_DIR="$2"
                shift 2
                ;;
            --single_end)
                PAIRED_END=false
                shift
                ;;
            --submit)
                SUBMIT_JOB=true
                shift
                ;;
            *)
                echo "Unknown parameter: $1"
                exit 1
                ;;
        esac
    done

    if [ -z "$SAMPLE_NAME" ]; then
        echo -e "\nUsage: run_star_alignment --sample_name <SAMPLE_NAME> [--input_dir <INPUT_DIR>] [--genome_dir <GENOME_DIR>] [--output_dir <OUTPUT_DIR>] [--single_end] [--submit]\n"
        echo "Requirements:"
        echo "  - STAR aligner installed and accessible from PATH"
        echo "  - Input files located in the specified input directory or default: $INPUT_DIR"
        echo "  - Genome directory specified or default: $GENOME_DIR"
        echo "  - Output directory specified or default: $OUTPUT_DIR"
        echo "  - Specify single-end with '--single_end'. Default assumes Paired-End data"
        echo "  - Use '--submit' to submit the SLURM job after creation"
        printf "\nExample: run_star_alignment.sh --sample_name Sample1 --submit\n\n"
        return 1
    fi

    # Ensure output directory exists
    mkdir -p "$OUTPUT_DIR/slurm_logs"

    # Define input files based on single or paired-end
    local READ_FILES=""
    if [ "$PAIRED_END" = true ]; then
        READ_FILES="$INPUT_DIR/${SAMPLE_NAME}_R1.fastq $INPUT_DIR/${SAMPLE_NAME}_R2.fastq"
    else
        READ_FILES="$INPUT_DIR/${SAMPLE_NAME}_R1.fastq"
    fi

    # Make SLURM script unique by appending a timestamp
    local TIMESTAMP=$(date +%Y%m%d)

    # SLURM script content
    local SLURM_SCRIPT="$OUTPUT_DIR/slurm_logs/star_mapping_${TIMESTAMP}_${SAMPLE_NAME}.slurm"
    cat <<EOT > "$SLURM_SCRIPT"
#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --job-name=star_mapping_${SAMPLE_NAME}
#SBATCH --output=${OUTPUT_DIR}/slurm_logs/star_mapping_${TIMESTAMP}_${SAMPLE_NAME}.out
#SBATCH --error=${OUTPUT_DIR}/slurm_logs/star_mapping_${TIMESTAMP}_${SAMPLE_NAME}.err

STAR --runThreadN 16 \
     --genomeDir "$GENOME_DIR" \
     --readFilesIn $READ_FILES \
     --genomeLoad LoadAndRemove \
     --outFilterType BySJout \
     --outFilterMultimapNmax 10 \
     --alignSJoverhangMin 8 \
     --alignSJDBoverhangMin 1 \
     --outFilterMismatchNmax 4 \
     --alignIntronMin 20 \
     --alignIntronMax 1000000 \
     --alignMatesGapMax 100000 \
     --genomeFileSizes "$GENOME_DIR/chrLength.txt" \
     --outFileNamePrefix "$OUTPUT_DIR/${SAMPLE_NAME}-"
EOT

    echo "SLURM job script generated: $SLURM_SCRIPT"

    if [ "$SUBMIT_JOB" = true ]; then
        sbatch "$SLURM_SCRIPT"
    fi
}

if [[ "$0" == "$BASH_SOURCE" ]]; then
    run_star_alignment "$@"
fi
