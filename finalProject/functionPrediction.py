from sklearn.linear_model import LogisticRegression
from sklearn.multioutput import MultiOutputClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report
from tqdm import tqdm
import joblib
import numpy as np
from scipy.sparse import lil_matrix

# Paths to data files
sequence_file = '../data/New_filtered_balanced_protein_language_data_34567_top2000'
go_table_file = '../data/F_3_GO_table.DAT'
vocab_file = '../data/New_Balanced_FragDATABASE_34567_top2000.txt'

# Load GO terms
print("Mapping function word to GO term...")
go_terms = {}
with open(go_table_file, 'r') as f:
    for line in f:
        go_term, func_word = line.strip().split('\t')
        go_terms[func_word] = go_term

# Load vocabulary
print("Building vocabulary and index mapping...")
with open(vocab_file, 'r') as f:
    vocab = [line.strip() for line in f.readlines()]
vocab_dict = {word: idx for idx, word in enumerate(vocab)}

func_words = list(go_terms.keys())
func_dict = {func_word: idx for idx, func_word in enumerate(func_words)}

# Determine number of sequences
print("Determining the number of sequences...")
with open(sequence_file, 'r') as f:
    total_sequences = sum(1 for line in f)

num_sequences = total_sequences
num_features = len(vocab)
num_labels = len(func_words)

# Initialize feature and label matrices
X = lil_matrix((num_sequences, num_features), dtype=np.int8)
y = lil_matrix((num_sequences, num_labels), dtype=np.int8)

# Convert sequence to feature vector
def sequence_to_features(sequence, vocab_dict, num_features):
    features = np.zeros(num_features, dtype=np.int8)
    for word in sequence:
        if word in vocab_dict:
            features[vocab_dict[word]] = 1
    return features

# Process data
print("Processing sequences and functions...")
with open(sequence_file, 'r') as f:
    for i, line in enumerate(tqdm(f, desc="Processing sequences", total=num_sequences)):
        if i >= num_sequences:
            break
        sequence, func = line.strip().split('|')
        sequence = sequence.split()
        func_list = func.split()
        feature_vector = sequence_to_features(sequence, vocab_dict, num_features)
        X[i, :] = feature_vector
        for func_word in func_list:
            if func_word in func_dict:
                y[i, func_dict[func_word]] = 1

X = X.tocsr()
y = y.tocsr()

# Split into training and testing sets
print("Splitting into training and testing sets...")
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Ensure that X is in dense format
X_train_dense = X_train.toarray() if hasattr(X_train, 'toarray') else X_train
X_test_dense = X_test.toarray() if hasattr(X_test, 'toarray') else X_test

# Convert y to a dense numpy array
y_train_dense = y_train.toarray()
y_test_dense = y_test.toarray()

# Train model using LogisticRegression with MultiOutputClassifier
print("Training model...")
model = LogisticRegression(max_iter=100, random_state=42)
multi_label_model = MultiOutputClassifier(model)

# Fit the model with dense data
multi_label_model.fit(X_train_dense, y_train_dense)

# Save the model
print("Saving model...")
model_file = 'trained_model.pkl'
joblib.dump(multi_label_model, model_file)
print(f"Model saved to {model_file}")

# Evaluate the model
print("Evaluating model...")
y_pred = multi_label_model.predict(X_test_dense)
print(classification_report(y_test_dense, y_pred, target_names=func_words))
