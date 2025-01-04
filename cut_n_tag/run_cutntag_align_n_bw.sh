#!/bin/bash
set -e  # This ensures the script exits on errors
set -u  # Treat unset variables as an error

run_cutntag_align_n_bw() {
    # Function to generate a SLURM job for BOWTIE alignment for a single sample
    # Arguments:
    #   --sample_name: Sample name (e.g., "Sample1")
    #   --input_dir: Directory containing input fastq files (default: predefined path)
    #   --bowtie_idx: STAR genome directory (default: predefined path)
    #   --output_dir: Directory for output files (default: predefined path)
    #   --single_end: Flag Indicating data is not paired-end (default: assumes Paired-End data)
    #   --cpus: How many CPUs to allocate for job (default: 16)

    # Default values for directories
    local BOWTIE_IDX="/projects/p31767/Bioinformatics/mouse/mm10/mm10_bowtie_index/mm10"
    local INPUT_DIR="/projects/b1042/Shukla_lab/Daniela/projects/ngs_default_analysis/input_files"
    local OUTPUT_DIR="/projects/b1042/Shukla_lab/Daniela/projects/ngs_default_analysis/output_analysis"
    local SAMPLE_NAME=""
    local PAIRED_END=true
    local SUBMIT_JOB=false
    local CPUs=16

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
            --bowtie_idx)
                BOWTIE_IDX="$2"
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
            --cpus)
                CPUs="$2"
                shift 2
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
        echo -e "\nUsage: run_cutntag_align_n_bw.sh --sample_name <SAMPLE_NAME> [--bowtie_idx <BOWTIE_IDX>] [--input_dir <INPUT_DIR>] [--output_dir <OUTPUT_DIR>] [--single_end] [--submit]\n"
        echo "Requirements:"
        echo "  - Preffix of sample name within input directory. e.g. INPUT_DIR/\${SAMPLE_NAME}_R1_001.fastq"
        echo "Arguments:"
        echo "  - Bowtie index directory specified or default: $BOWTIE_IDX"
        echo "  - Input files located in the specified input directory or default: $INPUT_DIR"
        echo "  - Output directory specified or default: $OUTPUT_DIR"
        echo "  - Specify single-end with '--single_end'. Default assumes Paired-End data"
        echo "  - Use '--submit' to submit the SLURM job after creation"
        printf "\nExample: run_cutntag_align_n_bw.sh --sample_name Sample1 --submit\n\n"
        return 1
    fi

    # Ensure output directory exists
    mkdir -p "$OUTPUT_DIR/slurm_logs"

    # Check if input file exists
    if [ ! -f "${INPUT_DIR}/${SAMPLE_NAME}_R1_001.fastq" ]; then
        echo "Error: Input file ${INPUT_DIR}/${SAMPLE_NAME}_R1_001.fastq not found."
        exit 1
    fi
    if [ "$PAIRED_END" = true ]; then
        if [ ! -f "${INPUT_DIR}/${SAMPLE_NAME}_R2_001.fastq" ]; then
            echo "Error: Input file ${INPUT_DIR}/${SAMPLE_NAME}_R2_001.fastq not found."
            exit 1
        fi
    fi


    # Define input files based on single or paired-end
    local READ_FILES=""
    if [ "$PAIRED_END" = true ]; then
        READ_FILES="$INPUT_DIR/${SAMPLE_NAME}_R1_001.fastq -2 $INPUT_DIR/${SAMPLE_NAME}_R2_001.fastq"
    else
        READ_FILES="$INPUT_DIR/${SAMPLE_NAME}_R1_001.fastq"
    fi

    # Make SLURM script unique by appending a timestamp
    local TIMESTAMP=$(date +%Y%m%d)

    # SLURM script content
    local SLURM_SCRIPT="$OUTPUT_DIR/slurm_logs/cutntag_align_n_bw_${TIMESTAMP}_${SAMPLE_NAME}.slurm"
    cat <<EOT > "$SLURM_SCRIPT"
#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --nodes=1
#SBATCH --cpus-per-task=${CPUs}
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --job-name=cutntag_align_n_bw_${SAMPLE_NAME}
#SBATCH --output=${OUTPUT_DIR}/slurm_logs/cutntag_align_n_bw_${TIMESTAMP}_${SAMPLE_NAME}.out
#SBATCH --error=${OUTPUT_DIR}/slurm_logs/cutntag_align_n_bw_${TIMESTAMP}_${SAMPLE_NAME}.err

set -e  # This ensures the script exits on errors
set -u  # Treat unset variables as an error

SAMPLE_HANDLE=${OUTPUT_DIR}/${SAMPLE_NAME}

# Map to reference genome using Bowtie
module load bowtie
bowtie -p ${CPUs} \
    -S \
    --fr \
    --chunkmbs 1024 \
    ${BOWTIE_IDX} \
    -1 ${READ_FILES} \
    -S ${SAMPLE_HANDLE}.mapped.sam


# Binarize, sort and index alignment
module load samtools/1.10.1
samtools view \
    -bS ${SAMPLE_HANDLE}.mapped.sam \
    > ${SAMPLE_HANDLE}.mapped.bam
samtools sort \
    -@ ${CPUs} \
    -o ${SAMPLE_HANDLE}.sorted.bam \
    ${SAMPLE_HANDLE}.mapped.bam
samtools index ${SAMPLE_HANDLE}.sorted.bam


# Mark and Remove PCR duplicates:
module load java/jdk1.8.0_191
module load picard/2.21.4
picard MarkDuplicates \
    I=${SAMPLE_HANDLE}.sorted.bam \
    O=${SAMPLE_HANDLE}.sorted_noPCR.bam \
    M=${OUTPUT_DIR}/metrics_${SAMPLE_NAME}.txt \
    REMOVE_DUPLICATES=true
samtools index ${SAMPLE_HANDLE}.sorted_noPCR.bam


# Calculate coverage 
# # with PCRs
module load deeptools
bamCoverage \
    --numberOfProcessors ${CPUs} \
    -b ${SAMPLE_HANDLE}.sorted.bam \
    --normalizeUsing CPM \
    -o ${SAMPLE_HANDLE}.sorted.bw

# # withOut PCRs
bamCoverage \
    --numberOfProcessors ${CPUs} \
    -b ${SAMPLE_HANDLE}.sorted_noPCR.bam \
    --normalizeUsing CPM \
    -o ${SAMPLE_HANDLE}.sorted_noPCR.bw


# Remove raw files:
rm ${SAMPLE_HANDLE}.mapped.sam
rm ${SAMPLE_HANDLE}.mapped.bam


EOT

    echo "SLURM job script generated: $SLURM_SCRIPT"

    if [ "$SUBMIT_JOB" = true ]; then
        sbatch "$SLURM_SCRIPT"
    fi
}

if [[ "$0" == "$BASH_SOURCE" ]]; then
    run_cutntag_align_n_bw "$@"
fi
