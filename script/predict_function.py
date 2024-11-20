import os
import pickle
import argparse
from collections import defaultdict

# 1. Load pre-calculated probabilities
def load_probabilities(filename="probabilities.pkl"):
    if os.path.exists(filename):
        with open(filename, 'rb') as file:
            return pickle.load(file)
    else:
        print(f"Error: {filename} does not exist. Please calculate the probabilities first.")
        return None

# 2. Read FASTA file (input protein sequence)
def read_fasta(file_path):
    with open(file_path, 'r') as file:
        return file.read().replace("\n", "").strip()

# 3. Predict protein function based on input sequence
def predict_function(input_sequence, sequence_words, probabilities):
    function_probabilities = defaultdict(float)

    for sequence in sequence_words:
        if sequence in input_sequence:
            for function in probabilities[sequence]:
                function_probabilities[function] += probabilities[sequence][function]

    # Normalize the probabilities
    total_prob = sum(function_probabilities.values())
    if total_prob > 0:
        for function in function_probabilities:
            function_probabilities[function] /= total_prob

    return function_probabilities

# 4. Load GO term mappings
def load_go_mapping(file_path):
    go_mapping = {}
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split()
            go_mapping[parts[1]] = parts[0]  # Mapping function to GO term

    return go_mapping

# 5. Main Prediction Workflow
def main(fasta_file):
    # Load pre-calculated probabilities from the file
    probabilities = load_probabilities()

    # If probabilities are not loaded, stop the script
    if probabilities is None:
        return

    # Load the sequence words from the file
    with open('../data/New_Balanced_FragDATABASE_34567_top2000.txt', 'r') as file:
        sequence_words = file.read().splitlines()

    # Read the input protein sequence (FASTA format)
    input_sequence = read_fasta(fasta_file)

    # Get function probabilities based on the input sequence
    function_probabilities = predict_function(input_sequence, sequence_words, probabilities)

    # Load GO term mappings
    go_mapping = load_go_mapping('../data/F_3_GO_table.DAT')

    # Find the top three functions with the highest probabilities
    if function_probabilities:
    # Sort functions by probability in descending order
        sorted_functions = sorted(function_probabilities.items(), key=lambda x: x[1], reverse=True)

    # Get the top three functions (or fewer if there are less than three)
        top_three = sorted_functions[:3]

        print("Top Predicted Functions:")
        for function, probability in top_three:
            go_term = go_mapping.get(function, 'Unknown GO term')
            print(f"Function: {function}, GO Term: {go_term}, Probability: {round(probability, 3)}")
    else:
        print("No valid functions predicted.")


if __name__ == "__main__":
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Predict protein function from FASTA file.")
    parser.add_argument('fasta_file', type=str, help="Path to the input FASTA file.")
    
    # Parse the arguments
    args = parser.parse_args()

    # Run the main function with the provided FASTA file
    main(args.fasta_file)
