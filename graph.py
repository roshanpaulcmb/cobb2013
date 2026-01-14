#### IMPORTS ####
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#### FUNCTIONS ####

def readData(file):
    df = pd.read_csv(file)
    return df

def compareBarplot(csv1, csv2, feature, label1="CSV 1", label2="CSV 2"):
    df1 = readData(csv1)
    df2 = readData(csv2)

    if feature not in df1.columns:
        raise ValueError(f"{feature} not found in {csv1}")
    if feature not in df2.columns:
        raise ValueError(f"{feature} not found in {csv2}")

    x = np.arange(len(df1))
    width = 0.4

    plt.figure(figsize=(12, 5))
    plt.bar(x - width/2, df1[feature], width, label=label1)
    plt.bar(x + width/2, df2[feature], width, label=label2)

    plt.xticks(x, df1["file"], rotation=90)
    plt.ylabel(feature)
    plt.title(f"{feature} comparison")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{feature}.png")
    plt.close()
    
    return None


#### RUN ####

def runAll():
    # alphafold
    afGr = "AlphaFold/gyrationRadii.csv"
    afSasa = "AlphaFold/sasaLeeRichards.csv"
    
    # esm
    esmGr = "ESM/gyrationRadii.csv"
    esmSasa = "ESM/sasaLeeRichards.csv"
    
    compareBarplot(afGr, esmGr, feature = "gyrationRadius",
                   label1 = "AlphaFold", label2 = "ESM")
    compareBarplot(afSasa, esmSasa, feature = "Polar",
                   label1 = "AlphaFold", label2 = "ESM")
    compareBarplot(afSasa, esmSasa, feature = "Apolar",
                   label1 = "AlphaFold", label2 = "ESM")
    
    return None

if __name__ == "__main__":
    runAll()