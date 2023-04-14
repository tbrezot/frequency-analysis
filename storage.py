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


def plot_padding_space(n, max_B):
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
    i = 0
    for (i, f) in enumerate(distr[:len(distr) - (s - 1)]):
        if math.ceil(V * f / B) == math.ceil(V * distr[i + (s - 1)] / B):
            break
    return i


def plot_padded_distribution(distr, V, B, s):
    """
    Plot the given distribution and its padded version. Normalize them for
    comparison purpose.

    Ars:
        - `V` (int)             : number of indexed values
        - `B` (int)             : padding number
        - `distr` (List(float)) : distribution
    """

    # Padded volume per keyword.
    line_size = L_uid + L_enc + math.ceil(B/8) + B * (1 + L_b)
    padded_storage = [line_size * math.ceil(V * f / B) for f in distr]

    x = [(i + 1) for i in range(len(distr))]
    plt.vlines(get_first_secure_keyword(distr, V, B, s) + 1, 0, 1)
    plt.plot(x, [f / distr[0] for f in distr])
    plt.plot(x, [v_wi / padded_storage[0] for v_wi in padded_storage], "x")
    plt.plot()
    plt.show()


def plot_secure_index(distr, V, B, s):

    B_range = [B for B in range(1, 101)]
    k_safe_B = [get_first_secure_keyword(distr, V, B, s) for B in B_range]

    plt.title(f"Secure portion of the index given B (s={s}, n={n}, V={V})")
    plt.ylabel(f"normalized indice of the first secure keyword (k/n)")
    plt.xlabel(f"B")
    plt.plot(B_range, [(len(distr) - y) / len(distr) for y in k_safe_B])
    plt.show()

    V_range = [n * i for i in range(1, 1001)]
    k_safe_V = [get_first_secure_keyword(
        distr, i * len(distr), B, s) for i in V_range]

    plt.title(f"Secure portion of the index given V (s={s}, n={n}, B={B})")
    plt.ylabel(f"normalized indice of the first secure keyword (k/n)")
    plt.xlabel(f"V/n")
    plt.plot([x / len(distr) for x in V_range],
             [(len(distr) - y) / len(distr) for y in k_safe_V])
    plt.show()


if __name__ == "__main__":
    # Number of keywords indexed.
    n = 1_000

    # Maximum value for B
    max_B = 100

    plot_padded_distribution(distr.zipf(n), 10 * n, 6, 2)
    plot_secure_index(distr.zipf(n), 10 * n, 6, 2)
    plot_padding_space(n, max_B)