import os
from concurrent.futures import ThreadPoolExecutor
import argparse

def get_options():
    description = 'Concatenate single aligned genes into a msa'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('input',
                        help='Directory containing the single aligned gene files with the .aln.fa extension')
    parser.add_argument('output',
                        help='Output file name')
    parser.add_argument('--threads',
                        default=2,
                        type=int,
                        help='Number of cores')

    return parser.parse_args()

if __name__ == "__main__":
    options = get_options()

    def read_alignment_file(file_path):
        """Reads an alignment file and returns a dictionary of sample IDs and gene sequences."""
        sequences = {}
        with open(file_path, 'r') as file:
            lines = file.readlines()
            for line in lines:
                if line.startswith('>'):
                    header = line.strip()[1:]
                    sample_id = header.split(';')[0].strip()
                    sequences[sample_id] = sequences.get(sample_id, '') + lines[lines.index(line) + 1].strip()
        return sequences

    def concatenate_alignments(directory, max_workers=None):
        """Concatenates gene sequences based on sample IDs from multiple alignment files in a directory."""
        all_sequences = {}
        file_paths = [os.path.join(directory, filename) for filename in os.listdir(directory) if filename.endswith('.aln.fa')]

        print(f"Processing {len(file_paths)} files with a maximum of {options.threads} cores...")

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            individual_sequences = list(executor.map(read_alignment_file, file_paths))

        print("Combining individual sequences...")

        for sequences in individual_sequences:
            for sample_id, sequence in sequences.items():
                all_sequences[sample_id] = all_sequences.get(sample_id, '') + sequence

        print("Writing concatenated sequences to the output file...")

        return all_sequences

    def write_concatenated_sequences(output_file, sequences):
        """Writes concatenated gene sequences to an output file."""
        with open(output_file, 'w') as file:
            for sample_id, sequence in sequences.items():
                file.write(f'>{sample_id}\n{sequence}\n')

    print(f"Concatenated sequences written to {options.output}")

    max_workers = options.threads
    concatenated_sequences = concatenate_alignments(options.input)
    write_concatenated_sequences(options.output, concatenated_sequences)
