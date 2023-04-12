# Combination algorithm

I order to compute the best guess of a distribution:
- for each sample compute its distance with each point $$(x, y)$$ in the distribution;
- for each $$x \in \Omega$$, select a sample and add the distances to the score;
- the guess with the best score is selected.

## Problems

1. How many samples it is necessary to draw to be able to make a guess?
