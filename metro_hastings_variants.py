import numpy as np
from ecoli_seq_variants import all_aminos
from funtrp_blosum_predict import (
    fragment_replace_mtx
    )

num_aminos = len(all_aminos)

def int_encode_variant(variant):
    return tuple((all_aminos.index(c) for c in variant))

class MetroHastingsVariants:
    def __init__(self, fragment, funtrp_triples, result_set_size):
        self.fragment = fragment
        self.num_frag_aas = len(fragment)
        self.variant_scores = {}
        self.orig_coded_frag = self.int_encode_variant(fragment)
        self.replace_mtx = fragment_replace_mtx(fragment, funtrp_triples)
        self.score_frag(self.orig_coded_frag) # initialize known variant scores with original fragment
        self.log_replace_mtx = np.log(self.replace_mtx)
        self.best_variants = set()
        self.hist = np.zeros_like(self.replace_mtx)
        self.result_set_size = result_set_size

    def int_encode_variant(self, variant):
        return int_encode_variant(variant)

    def score_frag(self, coded_frag):
        if coded_frag in self.variant_scores:
            return self.variant_scores[coded_frag]
        score = sum((self.replace_mtx[i,j] for j, i in enumerate(coded_frag)))
        self.variant_scores[coded_frag] = score
        return score

    def decode_variant(self, coded_var):
        return ''.join((all_aminos[j] for j in coded_var))

    # Metropolis-Hastings from https://towardsdatascience.com/from-scratch-bayesian-inference-markov-chain-monte-carlo-and-metropolis-hastings-in-python-ef21a29e25a
    # The transition model defines how to move from the current fragment to a new fragment
    def transition_model(self, x, # randomly mutate one amino acid.
            idx_choices=[], # make fewer objects by creating local variables
            amino_choices=np.array(range(num_aminos))):
        if not idx_choices:
            idx_choices.append(np.array(range(self.num_frag_aas)))
        idx = np.random.choice(idx_choices[0])
        newx = list(x[:idx])
        newx.append(np.random.choice(amino_choices))
        newx.extend(x[idx + 1:])
        return tuple(newx)

    #Computes the likelihood of the data given a fragment (new or current) according to equation (2)
    def likelihood_computer(self, x):
        return sum((self.log_replace_mtx[i,j] for j, i in enumerate(x)))

    #Defines whether to accept or reject the new sample
    def acceptance_rule(self, x, x_new):
        if x_new > x:
            return True
        else:
            accept = np.random.uniform(0,1)
            # Since we did a log likelihood, we need to exponentiate in order to compare to the random number
            # less likely x_new are less likely to be accepted
            return accept < np.exp(x_new - x)
            #return accept*x < x_new

    def end_condition(self):
        return len(self.best_variants) >= self.result_set_size

    def run(self, min_interval=10):
        # likelihood_computer(x): returns the likelihood of this sample
        # transition_model(x): a function that draws a sample from a symmetric distribution and returns it
        # param_init: a starting sample
        # iterations: number of accepted to generated
        # acceptance_rule(x,x_new): decides whether to accept or reject the new sample
        if self.best_variants:
            raise RuntimeError("MetroHastingsVariants can only be run once")
        x = self.orig_coded_frag
        accepted = []
        rejected = []
        steps_since_added = 0
        while not self.end_condition():
            x_new =  self.transition_model(x)
            x_lik = self.likelihood_computer(x)
            x_new_lik = self.likelihood_computer(x_new)
            steps_since_added += 1
            if self.acceptance_rule(x_lik, x_new_lik):
                x = x_new
                if steps_since_added < min_interval:
                    continue
                decoded_var = self.decode_variant(x_new)
                if decoded_var not in self.best_variants:
                    accepted.append(x_new)
                    steps_since_added = 0
                    for j, i in enumerate(x_new):
                        self.hist[i, j] += 1
                self.best_variants.add(decoded_var)
            else:
                rejected.append(x_new)
        return np.array(accepted), np.array(rejected) # primarily to debug if needed

    def result_set(self):
        if not self.best_variants:
            self.run()
        return self.best_variants


if __name__ == '__main__':
    accepted, rejected = metropolis_hastings(manual_log_like_normal, transition_model,
            orig_coded_frag, acceptance, end_condition)

    print([decode_variant(f) for f in accepted[:100]])
    print('accepted fraction', len(accepted)/(len(accepted) + len(rejected)), 'of proposals')
    #top = np.array(top)
    top = replace_mtx.ravel('F')
    scale = .01
    top = scale*top/(scale+top) # flatten out values that are way bigger than most others
    bottom = np.zeros_like(top)
    width = .5
    depth = 1

    # This import registers the 3D projection, but is otherwise unused.
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

    # setup the figure and axes
    _x = np.arange(num_aminos)
    _y = np.arange(num_frag_aas)
    _xx, _yy = np.meshgrid(_x, _y)
    x, y = _xx.ravel(), _yy.ravel()
    fig = plt.figure(figsize=(12, 8))
    ax1 = fig.add_subplot(111, projection='3d')
    ax1.bar3d(x, y, bottom, width, depth, top, shade=True, alpha=1)
    top = (hist/sum(hist)).ravel('F')/25
    top = scale*top/(scale+top) # flatten out values that are way bigger than most others
    ax1.bar3d(x, y + 25, bottom, width, depth, top, shade=True, alpha=1, color='green')
    ax1.set_title(''.join(all_aminos))
    plt.xticks(_x+.5, all_aminos)
    plt.yticks(_y+.5, fragment)


    plt.show()
