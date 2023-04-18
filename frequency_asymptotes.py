import math
import matplotlib.pyplot as plt
import distribution as distr

from sage.all import binomial, RealNumber


def compute_zipf_k_snN(s, n, N):
    """
    Computes the safety asymptotes given `s` the security parameter (the
    least number of occurrences that ensures safety of the data).

    The theoretical asymptote is computed by solving the following equation:

        ```
            N (f_k - f_{k + s}) < 1
        <=> k**2 + s * k - C_nN < 0
        <=> k = root_1 or root_2
            and root_{1,2} = (-s +/- sqrt(s**2 + 4 * s * C_nN)) / 2
        ```

    where `N` is the number of samples drawn and `f_k` is the theoretical
    frequency of the `k`th word, and `C_nN = N / sum(1/(i))`.
    """
    C_nN = N / sum([1 / (i + 1) for i in range(n)])
    return math.ceil((-s + math.sqrt(s**2 + 4 * s * C_nN) / 2))


def get_inversion_probability(freq, ifreq, k1, k2, N, coefficients):
    """
    Returns the probability that the `k1`th and `k2`th elements of the
    distribution are reversed in a sample of N elements.
    """
    p = 0
    for t1 in range(N+1):
        p1 = coefficients[t1] * freq[k1][t1] * ifreq[k1][N - t1]
        p2 = 0
        for t2 in range(t1, N+1):
            p2 += coefficients[t2] * freq[k2][t2] * ifreq[k2][N - t2]
        p += p1 * p2
    print(f"{k1} vs {k2} inversion: {p}")
    return p


def zipf_binomial_study(s, n,  N):
    # Compute `k_snN`.
    k_snN = compute_zipf_k_snN(s, n, N)
    print(f"k_snN = {k_snN}")

    # Build the Zipf distribution for `n` elements.
    # Use the Sage `RealNumber` to have enough precision.
    zipf = [RealNumber(f) for f in distr.zipf(n)]

    # Build the matrix of frequency powers.
    freq = []
    ifreq = []
    for i in range(n):
        f = [1]
        inv_f = [1]
        for j in range(N):
            f.append(f[j] * zipf[i])
            inv_f.append(inv_f[j] * (1-zipf[i]))
        freq.append(f)
        ifreq.append(inv_f)

    # Compute the binomial coefficicents once and for all.
    coeffs = [binomial(N, i) for i in range(N+1)]

    # For all `k <= k_snN`, compute the probability of inversion with its
    # neighbour.
    inversions = [get_inversion_probability(
        freq, ifreq, k, k+1, N, coeffs) for k in range(k_snN)]

    assert len(inversions) == k_snN
    score = sum([(1-f) for f in inversions])/k_snN

    plt.vlines([k_snN], 0, inversions[len(inversions) - 1])
    plt.plot([(i+1) for i in range(k_snN)], inversions, "x")
    plt.title(f"Inversion probabilities (s={s}, n={n}, N={N}, score={score})")
    plt.legend([f"Esperance of the first safe keyword (k_snN={k_snN})",
                f"Probability of an inversion with the next keyword"])
    plt.show()


def get_first_safe_index(samples, s, N):
    """
    Finds the index of the first keyword of the first cluster which size is at
    least `s`.
    """
    samples = sorted(samples, reverse=True)
    k = len(samples) - 1
    count = 1
    while k > 0:
        if math.ceil(N * samples[k-1]) == math.ceil(N * samples[k]):
            count += 1
        elif count < s:
            # Keywords `(k-1)` and `k` are not in the same cluster and the
            # cluster of `k` is of size `count < s`. Since the distribution is
            # sorted and `k` ranges from `n - 1` to `0`, `k` is the first
            # keyword of the last unsafe cluster.
            break
        else:
            count = 1
        k -= 1

    # Return the index of the first keyword of the first safe cluster.
    return k + count


def draw_k_snN(s, n):
    C_sn = s / sum([1 / (i + 1) for i in range(n)])

    # The formula is valid for N in `[0; max_N]`.
    max_N = (n - s) * n / C_sn

    # Use the adimensionned value of N.
    adim_N_max = math.ceil(max_N / n)

    # Build the axis values.
    adim_N = [N for N in range(1, adim_N_max)]
    s_2 = s**2
    k_s = [(s + math.sqrt(s_2 + 4 * (n * adim_N) * C_sn)) / (2 * n)
           for adim_N in range(1, adim_N_max)]

    if plot_fig:
        plt.plot([x for x in adim_N], k_s)
        plt.ylabel("sensitive portion of the index (k_snN/n)")
        plt.xlabel("normialized number of index requests (N/n)")
        plt.title(
            "Unsecure portition of the index given the number of requests (s={}, n ={})".format(s, n))
        plt.show()


def draw_curve(s, n, N):
    zipf = distr.zipf(n)
    data = distr.sample_distr(zipf, N)

    k_snN_th = compute_zipf_k_snN(s, n, N)
    # Add one to converts array index to keyword indice.
    k_snN_emp = 1 + get_first_safe_index(data, s, N)

    if plot_fig:
        # Plot asymptotes.
        plt.vlines([k_snN_th/n], 0, max(data[0], zipf[0]),
                   linestyles='dotted')
        plt.vlines([k_snN_emp/n], 0,
                   max(data[0], zipf[0]), linestyles='dashed')

        # Plot Zipf and sampled distributions.
        x = [(i + 1) / n for i in range(n)]
        plt.plot(x, zipf, "x")
        plt.plot(x, data, "x")

        plt.legend(["last safe keyword (theoretical), k = {}".format(k_snN_th),
                    "last safe keyword (empirical), k = {}".format(k_snN_emp),
                   "distribution of the indexed keywords", "distribution of the index requests"])
        plt.ylabel("frequency (f_k)")
        plt.xlabel("normalized keyword index (k/n)")
        plt.title("Zipf distribution analysis (s={}, n={}, N={})".format(s, n, N))
        plt.show()

        plt.vlines([k_snN_th/n], 0, max(data[0], zipf[0]),
                   linestyles='dotted')
        plt.vlines([k_snN_emp/n], 0, max(data[0],
                   zipf[0]), linestyles='dashed')
        x = [(i + 1) / n for i in range(n)]
        plt.plot(x, sorted(data, reverse=True), "x")
        # plt.plot(x, [math.ceil(N * f) / N for f in zipf], "x")
        plt.ylabel("frequency (f_k)")
        plt.xlabel("normalized keyword index (k/n)")
        plt.title(
            "Sorted samples from the Zipf distribution (s={}, n={}, N={})".format(s, n, N))
        plt.legend(["last safe keyword (theoretical), k = {}".format(k_snN_th),
                    "last safe keyword (empirical), k = {}".format(k_snN_emp),
                    "distribution of the index requests (sorted by frequency)"])
        plt.show()


def main(get_score: bool):
    # Arbitrary value, the minimum number of words with the same theoretical
    # frequency for which an attack is considered benign.
    s = 2

    # Cardinal of the universe (nb of words indexed).
    n = 1_000

    # Number of called to the index.
    N = 2 * n

    # Do not run it each time since it is quite long.
    if get_score:
        zipf_binomial_study(s, n, N)

    # Draw the graph of `k_snN`.
    draw_k_snN(s, n)

    # Draw the theoretical and empirical curves.
    draw_curve(s, n, N)


if __name__ == "__main__":
    plot_fig = True

    # Set to `True` to compute the score of the attacker on `[0; k_snN]`.
    get_score = False

    main(get_score)
