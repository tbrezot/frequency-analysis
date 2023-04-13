import math
from numpy import infty
from frequency_analysis import get_last_unsafe_index


def make_a_guess(distr, sampled_distr, guess, score, used_indices, res):
    """
    Make a guess on the distribution that has not been done before and compute
    its score as follows:
        `score = sqrt(sum([distr[i] - guess[i]]))`
    """

    for i in range(len(sampled_distr)):
        if i not in used_indices:
            new_guess = guess.copy()
            new_guess.append(i)

            new_used_indices = used_indices.copy()
            new_used_indices.add(i)

            score += (sampled_distr[i] - distr[len(guess) + 1])**2

            if len(new_guess) == len(sampled_distr):
                res.append((new_guess, math.sqrt(score)))
            elif len(new_guess) < len(sampled_distr):
                make_a_guess(distr, sampled_distr, new_guess,
                             score, new_used_indices, res)


def compute_best_guess(distr, sampled_distr, s, N):
    k_snN = get_last_unsafe_index(sampled_distr, s, N)
    print("k_snN =", k_snN)

    res = []
    make_a_guess(distr, sampled_distr[:k_snN], [], 0, set(), res)
    assert len(res), math.factorial(k_snN)

    best_guess = []
    best_score = infty
    for (guess, score) in res:
        if score < best_score:
            best_score = score
            best_guess = guess
    print("Best score found is:", best_score)
    for (i, j) in enumerate(best_guess):
        print("{} -> {} ({} vs {})".format(i, j, sampled_distr[j], distr[i]))
    return best_guess, best_score


