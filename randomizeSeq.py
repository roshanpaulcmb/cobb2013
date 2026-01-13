#### IMPORTS ####
import pandas as pd
import random


#### FUNCTIONS ####

def readData(file):
    df = pd.read_csv(file)
    return df


def createRandomSeq(df):
    
    aas = df["code"]
    probs = df["prob"]
    
    # generate random amino acid string based on likelihood of amino acid occuring in eukaryotes
    randomSeq = random.choices(aas, weights = probs, k=100)
    randomSeq = "".join(randomSeq)
    return randomSeq

def writeFasta(fileName, seq, header = None):
    with open(fileName, "w") as f:
        f.write(f">{header}\n")
        f.write(seq + "\n")
    return None

#### RUN ####

def runAll(file):
    df = readData(file)
    randomSeqs = []
    for i in range(10):
        randomSeqs.append(createRandomSeq(df))
        header = "random" + str(i)
        fileName = header + ".fasta"
        writeFasta(fileName, randomSeqs[i], header = header)
        print(i, randomSeqs[i])
    return randomSeqs

if __name__ == "__main__":
    file = "aa_likeliness.csv"
    runAll(file)