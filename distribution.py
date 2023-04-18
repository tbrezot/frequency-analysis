import random


def zipf(n):
    """
    Computes the theoretical Zipf distribution with `n` elements.
    """

    C_n = sum([1/(n + 1) for n in range(n)])
    zipf = [1/((i + 1) * C_n) for i in range(n)]

    # The distribution is already normalized.
    return zipf


def sample_distr(distr, N):
    """
    Draws N samples of data following the given distribution. Creates the
    sampled distribution.

    **Note**: a distribution is the list of the normalized probabilities for
    each element in the universe.
    """

    # Build the cumulative distribution.
    cumulative_distr = [distr[0]]
    for i in range(1, len(distr)):
        cumulative_distr.append(distr[i] + cumulative_distr[i-1])

    # Randomly draw `N` samples from the Zipf.
    data = [0 for _ in cumulative_distr]
    for _ in range(N):
        x = random.random()
        for (i, cumulative_freq) in enumerate(cumulative_distr):
            if x <= cumulative_freq:
                data[i] += 1
                break

    # Normalize the sampled distribution.
    res = [f/N for f in data]

    return res
