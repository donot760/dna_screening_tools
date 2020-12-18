import numpy as np
import matplotlib.pyplot as plt
from ecoli_seq_variants import all_aminos
from metro_hastings_variants import MetroHastingsVariants
from funtrp_blosum_predict import (
    parse_funtrp_line,
    )

if __name__ == '__main__':
    fragment = 'YANYEGCLWNATGVVVCTG'

    num_aminos = len(all_aminos)
    num_frag_aas = len(fragment)

    dummy_funtrp = [((1 if y_coord<=5 else 0), (1 if 5<y_coord<=12 else 0), (1 if 12<y_coord else 0)) for y_coord in range(num_frag_aas)]
    yanyeg_n = [0.04, 0.02, 0.05, 0.2, 0.35, 0.06, 0.01, 0.11, 0.12, 0.2, 0.16, 0.09, 0.09, 0.13, 0.12, 0.11, 0.11, 0.15, 0.37]
    yanyeg_t = [0.55, 0.63, 0.58, 0.14, 0.21, 0.31, 0.77, 0.04, 0.17, 0.22, 0.34, 0.25, 0.17, 0.23, 0.07, 0.08, 0.26, 0.08, 0.1]
    yanyeg_r = [0.41, 0.35, 0.37, 0.66, 0.44, 0.63, 0.22, 0.85, 0.71, 0.58, 0.5, 0.66, 0.74, 0.64, 0.81, 0.81, 0.63, 0.77, 0.53]
    real_funtrp = list(zip(yanyeg_n, yanyeg_t, yanyeg_r))

    variant_computer = MetroHastingsVariants(fragment, real_funtrp, 100)
    accepted, rejected = variant_computer.run()

    print(accepted[:100])
    print([v for _, v in zip(range(100), variant_computer.best_variants)])
    print('accepted fraction', len(accepted)/(len(accepted) + len(rejected)), 'of proposals')
    #top = np.array(top)
    top = variant_computer.replace_mtx.ravel('F')
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
