import math
import matplotlib.pyplot as plt
import distribution as distr


def compute_zipf_k_snN(s, n, N):
    """
    Computes the safety asymptotes given `s` the security parameter (the
    least number of occurrences that ensures safety of the data).

    The theoretical asymptote is computed by solving the following equation:

        ```
            N (f_i - f_{i + s}) = 1
        <=> i**2 + s * i - C_nN = 0
        <=> i = root_1 or root_2, root_1 < root_2
            and root_{1,2} = (s +/- sqrt(s**2 + 4 * s * C_nN)) / 2
        ```

    where `N` is the number of samples drawn and `f_i` is the theoretical
    frequency of the `i`th word, and `C_nN = N * s / sum(1/(i))`.

    The empirical asymptote is computed by finding the last group of samples
    that is smaller than `s`.
    """

    C_nN = N / sum([1 / (i + 1) for i in range(n)])
    return math.ceil(0.5 * (s + math.sqrt(s**2 + 4 * s * C_nN)))


def get_last_unsafe_index(sampled_distr, s, N):
    """
    Find the last unsafe index in the given samples.

    An index is defined as safe iff it is drawn the same number of times as at
    least (s - 1) other samples.
    """
    # Local sorted copy of the samples.
    n = len(sampled_distr)
    k = n - 1
    count = 1
    while k > 1:
        if math.ceil(N * sampled_distr[k-1]) == math.ceil(N * sampled_distr[k]):
            count += 1
        elif count < s:
            break
        else:
            count = 1
        k -= 1
    return k


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
    k_snN_emp = get_last_unsafe_index(data, s, N)

    if plot_fig:
        # Plot asymptotes.
        plt.vlines([k_snN_th/n], 0, max(data[0], zipf[0]), linestyles='dotted')
        plt.vlines([k_snN_emp/n], 0, max(data[0],
                   zipf[0]), linestyles='dashed')

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

        plt.vlines([k_snN_th/n], 0, max(data[0], zipf[0]), linestyles='dotted')
        x = [(i + 1) / n for i in range(n)]
        plt.plot(x, sorted(data, reverse=True), "x")
        plt.ylabel("frequency (f_k)")
        plt.xlabel("normalized keyword index (k/n)")
        plt.title(
            "Sorted samples from the Zipf distribution (s={}, n={}, N={})".format(s, n, N))
        plt.legend(["last safe keyword (theoretical), k = {}".format(k_snN_th),
                    "distribution of the index requests (sorted by frequency)"])
        plt.show()


def main():
    # Arbitrary value, the minimum number of words with the same theoretical
    # frequency for which an attack is considered benign.
    s = 2

    # Cardinal of the universe (nb of words indexed).
    n = 1_000

    # Draw the graph of `k_snN`.
    draw_k_snN(s, n)

    # Number of called to the index.
    N = 2 * n

    # Draw the theoretical and empirical curves.
    draw_curve(s, n, N)


if __name__ == "__main__":
    plot_fig = True
    main()