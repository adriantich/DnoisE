
import numpy as np
import pandas as pd


def entropy(x):
    freqs = pd.Series(x).value_counts()
    freqs = freqs / sum(freqs)
    shannon_entropy = -sum(freqs[freqs > 0] * np.log(freqs[freqs > 0]))
    return shannon_entropy


def mean_entropy(data):
    x = data.sequence.str.split('', expand=True, )
    entropy_values = x.apply(entropy, axis=0)
    first = np.mean(entropy_values.loc[range(1, len(entropy_values), 3)])
    second = np.mean(entropy_values.loc[range(2, len(entropy_values), 3)])
    third = np.mean(entropy_values.loc[range(3, len(entropy_values), 3)])
    return first, second, third
