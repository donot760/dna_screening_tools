import numpy as np
import matplotlib.pyplot as plt
from ecoli_seq_variants import all_aminos
from metro_hastings_variants import MetroHastingsVariants
from funtrp_blosum_predict import (
    fragment_replace_mtx,
    parse_funtrp_line,
    )

fragment = 'YANYEGCLWNATGVVVCTG'

num_aminos = len(all_aminos)
num_frag_aas = len(fragment)

dummy_funtrp = [((1 if y_coord<=5 else 0), (1 if 5<y_coord<=12 else 0), (1 if 12<y_coord else 0)) for y_coord in range(num_frag_aas)]
yanyeg_n = [0.04, 0.02, 0.05, 0.2, 0.35, 0.06, 0.01, 0.11, 0.12, 0.2, 0.16, 0.09, 0.09, 0.13, 0.12, 0.11, 0.11, 0.15, 0.37]
yanyeg_t = [0.55, 0.63, 0.58, 0.14, 0.21, 0.31, 0.77, 0.04, 0.17, 0.22, 0.34, 0.25, 0.17, 0.23, 0.07, 0.08, 0.26, 0.08, 0.1]
yanyeg_r = [0.41, 0.35, 0.37, 0.66, 0.44, 0.63, 0.22, 0.85, 0.71, 0.58, 0.5, 0.66, 0.74, 0.64, 0.81, 0.81, 0.63, 0.77, 0.53]
real_funtrp = list(zip(yanyeg_n, yanyeg_t, yanyeg_r))
replace_mtx = fragment_replace_mtx(fragment, real_funtrp)

variant_scores = {}
def score_frag(coded_frag):
    if coded_frag in variant_scores:
        return variant_scores[coded_frag]
    score = sum((replace_mtx[i,j] for j, i in enumerate(coded_frag)))
    variant_scores[coded_frag] = score
    return score

def int_encode_variant(variant):
    return tuple((all_aminos.index(c) for c in variant))

orig_coded_frag = int_encode_variant(fragment)
score_frag(orig_coded_frag) # initialize known variant scores with original fragment

def decode_variant(coded_var):
    return ''.join((all_aminos[j] for j in coded_var))

# Metropolis-Hastings from https://towardsdatascience.com/from-scratch-bayesian-inference-markov-chain-monte-carlo-and-metropolis-hastings-in-python-ef21a29e25a
# The transition model defines how to move from the current fragment to a new fragment
def transition_model(x): # randomly mutate one amino acid.
    idx = np.random.choice(range(num_frag_aas))
    return x[:idx] + (np.random.choice(range(num_aminos)),) + x[idx + 1:]

#def transition_model(x): # randomly mutate all amino acids.
#    return _transition_model(_transition_model(_transition_model(_transition_model(x))))

#Computes the likelihood of the data given a fragment (new or current) according to equation (2)
log_replace_mtx = np.log(replace_mtx)
def manual_log_like_normal(x):
    return sum((log_replace_mtx[i,j] for j, i in enumerate(x)))


#Defines whether to accept or reject the new sample
def acceptance(x, x_new):
    if x_new>x:
        return True
    else:
        accept=np.random.uniform(0,1)
        # Since we did a log likelihood, we need to exponentiate in order to compare to the random number
        # less likely x_new are less likely to be accepted
        return (accept < (np.exp(x_new-x)))
        #return accept*x < x_new

best_variants = set()
hist = np.zeros_like(replace_mtx)
def metropolis_hastings(likelihood_computer, transition_model, param_init,
        acceptance_rule, end_condition, min_interval=10):
    # likelihood_computer(x): returns the likelihood of this sample
    # transition_model(x): a function that draws a sample from a symmetric distribution and returns it
    # param_init: a starting sample
    # iterations: number of accepted to generated
    # acceptance_rule(x,x_new): decides whether to accept or reject the new sample
    x = param_init
    accepted = []
    rejected = []
    steps_since_added = 0
    while not end_condition():
        x_new =  transition_model(x)
        x_lik = likelihood_computer(x)
        x_new_lik = likelihood_computer(x_new)
        steps_since_added += 1
        if acceptance_rule(x_lik, x_new_lik):
            x = x_new
            if steps_since_added < min_interval:
                continue
            if x_new not in best_variants:
                accepted.append(x_new)
                steps_since_added = 0
                for j, i in enumerate(x_new):
                    hist[i, j] += 1
            best_variants.add(x_new)
        else:
            rejected.append(x_new)

    return np.array(accepted), np.array(rejected)

def end_condition():
    return len(best_variants) >= 100000


#best_variants = set([orig_coded_frag])
#def next_best_variant():
#    best_score = 0
#    best_variant = None
#    existing_variants = list(best_variants)
#    for variant in existing_variants:
#        variant_l = list(variant)
#        for position in range(num_frag_aas):
#            orig_val = variant[position]
#            for j in range(num_aminos):
#                variant_l[position] = j
#                test_variant = tuple(variant_l)
#                if test_variant in best_variants:
#                    continue
#                new_score = score_frag(test_variant)
#                variant_l[position] = orig_val
#                if new_score > best_score:
#                    best_score = new_score
#                    best_variant = test_variant
#    best_variants.add(best_variant)
#    return best_variant
#
#print(decode_variant(orig_coded_frag))
#for _ in range(100):
#    if _%10 == 0:
#        print(decode_variant(next_best_variant()))
#        print(len(variant_scores))
#    else:
#        next_best_variant()

    #top = [replace_score(all_aminos[x_coord], fragment[y_coord], ((1 if y_coord<=5 else 0), (1 if 5<y_coord<=12 else 0), (1 if 12<y_coord else 0))) for x_coord, y_coord in zip(x, y)]

if __name__ == '__main__':
    fragment = 'YANYEGCLWNATGVVVCTG'

    num_aminos = len(all_aminos)
    num_frag_aas = len(fragment)

    dummy_funtrp = [((1 if y_coord<=5 else 0), (1 if 5<y_coord<=12 else 0), (1 if 12<y_coord else 0)) for y_coord in range(num_frag_aas)]
    yanyeg_n = [0.04, 0.02, 0.05, 0.2, 0.35, 0.06, 0.01, 0.11, 0.12, 0.2, 0.16, 0.09, 0.09, 0.13, 0.12, 0.11, 0.11, 0.15, 0.37]
    yanyeg_t = [0.55, 0.63, 0.58, 0.14, 0.21, 0.31, 0.77, 0.04, 0.17, 0.22, 0.34, 0.25, 0.17, 0.23, 0.07, 0.08, 0.26, 0.08, 0.1]
    yanyeg_r = [0.41, 0.35, 0.37, 0.66, 0.44, 0.63, 0.22, 0.85, 0.71, 0.58, 0.5, 0.66, 0.74, 0.64, 0.81, 0.81, 0.63, 0.77, 0.53]
    real_funtrp = list(zip(yanyeg_n, yanyeg_t, yanyeg_r))

    variant_computer = MetroHastingsVariants(fragment, real_funtrp, 10000)
    accepted, rejected = variant_computer.run()

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
    top = (variant_computer.hist/sum(variant_computer.hist)).ravel('F')/25
    top = scale*top/(scale+top) # flatten out values that are way bigger than most others
    ax1.bar3d(x, y + 25, bottom, width, depth, top, shade=True, alpha=1, color='green')
    ax1.set_title(''.join(all_aminos))
    plt.xticks(_x+.5, all_aminos)
    plt.yticks(_y+.5, fragment)


    plt.show()
