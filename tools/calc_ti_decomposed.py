from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
import sys
from free_energy_base import TI, FreeEnergyCalculation

plt.rcParams['figure.subplot.bottom'] = 0.3


class TI_decomposed(object):
    """Estimate free differences using ThermodynamicIntegration. The class
    performs individual calculations for each component of the system
    interaction and internal energies. Results from each component
    will sum to give the same result as the total energy. The
    decomposition of the total free energy is not well-defined however
    can be indicative of the dominant terms in a free energy change.
    """
    def __init__(self):
        self.estimators = defaultdict(lambda: TI())

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

    def calculate(self, lambdas):
        """Calculate the free energy difference for each energy component and
        return a dictionary of PMF objects."""
        return {estimator: self.estimators[estimator].calculate(lambdas)
                for estimator in self.estimators}


if __name__ == "__main__":
    calc = FreeEnergyCalculation(sys.argv[1:], estimators=[TI_decomposed])
    results = calc.calculate()

    width = 0.8
    N = len(results[TI_decomposed])
    colours = ('b', 'g', 'r', 'y')
    fig, ax = plt.subplots()
    for i, repeat in enumerate(results[TI_decomposed]):
        xs = np.arange(len(repeat)) * N + width * i
        ys = [repeat[term].dG for term in repeat]
        labels = [term for term in repeat]
        ax.bar(xs, ys, color=colours[i])
        ax.set_xticks(xs)
        ax.set_xticklabels(labels, rotation="vertical")

    plt.legend(calc.root_paths)
    plt.savefig('decomposed.pdf')
