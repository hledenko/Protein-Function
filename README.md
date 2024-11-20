# Protein-Function
### Bioinformatics Final Project

**calculate-probabilties.py** :
- for each sequence, this function calculates the probability of a sequence being associated with a particular function by dividing the count of occurrences of that sequence-function pair by the total count of the sequence in the data.

**predict_function.py** :
- This script is responsible for predicting the function of a protein based on its sequence using the pre-calculated probabilities from the first script.

### To make a prediction:
1. download data files from (http://cs.plu.edu/~caora/cs487/Materials/functionResource.zip) and place in data directory

2. run **calculate_probabities.py** to generate probability dictionary

3. run **predict_function** to generate predicted GO terms like this: 
```bash
python predict_function.py ../fastaFiles/8v12.fasta
```
