import os
import pickle
from collections import defaultdict

def read_data(file_path, max_lines=500000):
    sequence_counts = defaultdict(int)
    sequence_function_counts = defaultdict(lambda: defaultdict(int))

    with open(file_path, 'r') as file:
        for line_number, line in enumerate(file, 1):
            # Stop processing after reading the first `max_lines` lines
            if line_number > max_lines:
                break

            parts = line.strip().split('|')
            sequences = parts[0].split()
            functions = parts[1].split()

            for sequence in sequences:
                sequence_counts[sequence] += 1
                for function in functions:
                    sequence_function_counts[sequence][function] += 1

    return sequence_counts, sequence_function_counts

def calculate_probabilities(sequence_counts, sequence_function_counts):
    # Use a regular dictionary instead of defaultdict for pickling compatibility
    probabilities = defaultdict(dict)

    for sequence in sequence_counts:
        for function in sequence_function_counts[sequence]:
            probabilities[sequence][function] = sequence_function_counts[sequence][function] / sequence_counts[sequence]

    return probabilities

def save_probabilities(probabilities, filename="probabilities.pkl"):
    with open(filename, 'wb') as file:
        pickle.dump(probabilities, file)

def main():
    # 1. Load the first 5000 lines of the data and calculate probabilities
    sequence_counts, sequence_function_counts = read_data('../data/New_filtered_balanced_protein_language_data_34567_top2000')
    probabilities = calculate_probabilities(sequence_counts, sequence_function_counts)

    # 2. Save the probabilities to a file
    save_probabilities(probabilities)
    print("Probabilities saved to 'probabilities.pkl'")

if __name__ == "__main__":
    main()
