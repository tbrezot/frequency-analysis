import math
import distribution as distr
import matplotlib.pyplot as plt

# Length of a Findex block.
L_b = 16
# Length of a Findex UID.
L_uid = 32
# AEAD encryption overhead.
L_enc = 28
# Size of the `K_wi`.
L_kwi = 16
# Size of a keyword hash.
L_hash = 32


def compute_actual_chain_table_storage(indexed_values, B):
    """
    Computes the Chain Table storage for the given volumes and `B` parameter.
    The `indexed_values` parameter lists the sizes of all value stored in all chains.
    """

    line_size = L_uid + L_enc + math.ceil(B/8) + B * (1 + L_b)

    n_lines = 0
    for v_i in indexed_values:
        n_blocks = 0
        for v_ij in v_i:
            n_blocks += math.ceil(v_ij/L_b)
        n_lines += math.ceil(n_blocks/B)

    return line_size * n_lines


def compute_ideal_chain_table_storage(indexed_values):
    """
    Computes the ideal Chain Table storage given the indexed values.
    The `indexed_values` parameter lists the sizes of all value stored in all chains.
    """
    return sum([L_uid + L_enc + sum([v_ij for v_ij in v_i]) for v_i in indexed_values])


def plot_padding_cost(n, max_B):
    """
    Plots the normalized padding for the given `n` and different indexed
    volumes `V`. The repartition of the indexed volume among the indexed
    keywords follows the Zipf distribution. Each indexed value fits in a single
    block.
    """

    zipf = distr.zipf(n)
    x = [i + 1 for i in range(max_B)]

    legend = []
    for i in range(4):
        legend.append(f"V = {10 ** i} * n")
        # Indexed volume.
        V = 10 ** i * n
        # Each indexed value fits in a single block.
        volumes = [[L_b for _ in range(math.ceil(V * f))] for f in zipf]
        size_opt = compute_ideal_chain_table_storage(volumes)
        y = []
        for B in range(1, max_B + 1):
            size_B = compute_actual_chain_table_storage(volumes, B)
            y.append((size_B - size_opt) / size_B)
        plt.plot(x, y)

    plt.title(
        f"Portion of the Chain Table used for padding given B for indexed values following a Zipf distribution (n = {n})")
    plt.ylabel("padding")
    plt.xlabel("B")
    plt.legend(legend)
    plt.show()


def get_first_secure_keyword(distr, V, B, s):
    """
    Computes the indice of the first keyword in the distribution that has been
    drawn the same number of times as (s-1) others.
    """
    i = 0
    for (i, f) in enumerate(distr[:len(distr) - (s - 1)]):
        if math.ceil(V * f / B) == math.ceil(V * distr[i + (s - 1)] / B):
            break
    return i


def plot_storage(distr, V, B, s):
    """
    Plots the given distribution and its padded version. Normalizes them for
    comparison purpose.

    Ars:
        - `V` (int)             : number of indexed values
        - `B` (int)             : padding number
        - `distr` (List(float)) : distribution
    """

    # Padded volume per keyword.
    line_size = L_uid + L_enc + math.ceil(B/8) + B * (1 + L_b)
    padded_storage = [line_size * math.ceil(V * f / B) for f in distr]

    k_safe = get_first_secure_keyword(distr, V, B, s) + 1
    plt.vlines(k_safe, 0, 1)

    x = [(i + 1) for i in range(len(distr))]
    plt.plot(x, [f / distr[0] for f in distr], "x")
    plt.plot(x, [v_wi / padded_storage[0] for v_wi in padded_storage], "x")

    plt.title(
        f"Comparison between theoretical and padded Chain Table storages (B={B}, V={V})")
    plt.ylabel(f"normalized storage")
    plt.xlabel(f"keyword")
    plt.legend([f"first secure keyword (k={k_safe}, s={s})",
               "theoretical storage", f"padded storage (B={B}, V={V})"])
    plt.show()


def plot_first_secure_index_given_B(distr, V, s):
    """
    Plots the evolution of the first safe keyword given B.
    """

    B_range = [B for B in range(1, 101)]
    k_safe_B = [get_first_secure_keyword(distr, V, B, s) for B in B_range]

    plt.title(f"Secure portion of the index given B (s={s}, n={n}, V={V})")
    plt.ylabel(f"secure portion of the index ((n - k)/n)")
    plt.xlabel(f"B")
    plt.plot(B_range, [(len(distr) - y) / len(distr) for y in k_safe_B], "x")
    plt.show()


def plot_first_secure_index_given_V(distr, B, s):
    """
    Plots the evolution of the first safe keyword given V.
    """

    V_range = [len(distr) * i for i in range(1, 1001)]
    k_safe_V = [get_first_secure_keyword(distr, V, B, s) for V in V_range]

    plt.title(f"Secure portion of the index given V (s={s}, n={n}, B={B})")
    plt.ylabel(f"secure portion of the index ((n - k)/n)")
    plt.xlabel(f"V/n")
    plt.plot([x / len(distr) for x in V_range],
             [(len(distr) - y) / len(distr) for y in k_safe_V],
             "x")
    plt.show()


if __name__ == "__main__":
    # Number of keywords indexed.
    n = 1_000

    # The number of values indexed per keyword follows the Zipf distribution.
    zipf = distr.zipf(n)

    plot_storage(zipf, 10 * n, 6, 2)

    plot_first_secure_index_given_B(zipf, 10 * n, 2)
    plot_first_secure_index_given_V(zipf, 6, 2)

    # Maximum value for B
    max_B = 100

    plot_padding_cost(n, max_B)
