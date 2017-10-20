from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
import sys
import free_energy_base as feb
# from free_energy_base import Estimator, FreeEnergyCalculation, TI

plt.rcParams['figure.subplot.bottom'] = 0.3


class TI_decomposed(feb.Estimator):
    """Estimate free differences using ThermodynamicIntegration. The class
    performs individual calculations for each component of the system
    interaction and internal energies. Results from each component
    will sum to give the same result as the total energy. The
    decomposition of the total free energy is not well-defined however
    can be indicative of the dominant terms in a free energy change.
    """
    def __init__(self, lambdas):
        self.estimators = defaultdict(lambda: feb.TI(lambdas))

    def add_data(self, series):
        """Save data from a SnapshotResults.series object
        for later calculation"""
        dlam = series.lamf[0]-series.lamb[0]
        for name, term in series.interaction_energies.iteritems():
            for energy in term[:-1]:
                self.estimators["%s_%s" % (name, energy.type)].data.append(
                    (energy.forw-energy.back) / dlam)

        for name, term in series.internal_energies.iteritems():
            for energy in term[:-1]:
                self.estimators["%s_%s" % (name, energy.type)].data.append(
                    (energy.forw-energy.back) / dlam)

    def calculate(self):
        """Calculate the free energy difference for each energy component and
        return a dictionary of PMF objects."""
        return {term: self.estimators[term].calculate()
                for term in self.estimators}

    def apply_slice(self, slc):
        for term in self.estimators:
            self.estimators[term].apply_slice(slc)

def get_arg_parser():
    parser = feb.get_arg_parser()
    return parser


if __name__ == "__main__":
    args = get_arg_parser().parse_args()
    calc = feb.FreeEnergyCalculation(args.directories, 
                                     estimators=[feb.TI, TI_decomposed])
    results = calc.calculate()

    width = 0.8
    N = len(results[TI_decomposed])
    colors = ('b', 'g', 'r', 'y')
    fig, ax = plt.subplots()
    for i, repeat in enumerate(results[TI_decomposed]):
        xs = np.arange(len(repeat)) * N + width * i
        ys, labels = [], []
        for term in sorted(repeat):
            ys.append(repeat[term].dG)
            labels.append(term)

        ax.bar(xs, ys, color=colors[i])
        ax.set_xticks(xs)
        ax.set_xticklabels(labels, 
                           rotation="vertical")

    for ti, ti_decomp in zip(results[feb.TI], results[TI_decomposed]):
        print "FDTI:", ti.dG,
        print "sum:", np.sum(ti_decomp[term].dG for term in ti_decomp)

    plt.legend(calc.root_paths)
    plt.savefig('decomposed.pdf')
