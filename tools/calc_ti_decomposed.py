import argparse
from collections import defaultdict
from itertools import chain
import numpy as np
import matplotlib.pyplot as plt
import free_energy_base as feb

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
        self.lambdas = lambdas

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

    def __getitem__(self, val):
        new_est = self.__class__(self.lambdas)
        for term in self.estimators:
            new_est.estimators[term] = self.estimators[term][val]
        return new_est


def get_result_means(repeats):
    """Take a list of estimator results and average all the terms.
    Return a dictionary mapping energy terms to single values."""
    return {term: np.mean([repeat[term].dG 
                           for repeat in repeats])
            for term in repeats[0]}
 

def get_arg_parser():
    """Add custom options for this script"""
    parser = argparse.ArgumentParser(
        description="Calculate individual contributions of different terms "
                    "to the total free energy difference. Although terms are "
                    "guaranteed to be additive with TI, the decomposition "
                    "is not strictly well defined. That said, it can be "
                    "illustrative to consider the dominant contributions of "
                    "a calculation.",
        parents=[feb.get_arg_parser()])
    parser.add_argument("-b", "--bound", nargs="+",
                        help="Output directories of bound phase calculations. "
                             "If present, it is assumed that -d "
                             "provides solvent phase simulation data.")
    parser.add_argument("-g", "--gas", nargs="+",
                        help="Output directories of gas phase calculations. "
                             "If present, it is assumed that -d "
                             "provides solvent phase simulation data.")
    parser.add_argument("-o", "--out", type=str, default="decomposed.pdf",
                        help="Filename of output graph.")
    return parser


if __name__ == "__main__":
    args = get_arg_parser().parse_args()
    calc = feb.FreeEnergyCalculation(args.directories,
                                     estimators=[feb.TI, TI_decomposed])
    # store results in a default dict for later simplicity of arithmetic
    results = calc.calculate()
    fdti_mean = np.mean([rep.dG for rep in results[feb.TI]])
    decomp_means = defaultdict(lambda: 0., 
                               get_result_means(results[TI_decomposed]))

    if (args.bound and args.gas) is not None:
        raise Exception("Cannot supply both bound and gas phase data.")
    elif (args.bound or args.gas) is not None:
        calc2 = feb.FreeEnergyCalculation(args.bound or args.gas,
                                          estimators=[feb.TI, TI_decomposed])
        results2 = defaultdict(lambda: 0., calc2.calculate())
        fdti_mean2 = np.mean([rep.dG for rep in results2[feb.TI]])
        decomp_means2 = defaultdict(lambda: 0., 
                                    get_result_means(results2[TI_decomposed]))

        # iterate over the unique keys from both calc and calc2
        for term in set(chain(decomp_means, decomp_means2)):
            decomp_means[term] -= decomp_means2[term]
            # swap the sign if this is a binding calculation
            if args.bound is not None:
                decomp_means[term] *= -1

        if args.gas is not None:
            fdti_mean -= fdti_mean2
        else:
            fdti_mean = fdti_mean2 - fdti_mean
            
    
    width = 0.8
    N = 1

    fig, ax = plt.subplots()
    xs = np.arange(len(decomp_means)) * N
    
    ys, labels = [], []
    for term in sorted(decomp_means):
        ys.append(decomp_means[term])
        labels.append(term)

    ax.bar(xs, ys)
    ax.set_xticks(xs)
    ax.set_xticklabels(labels, rotation="vertical")

    print "FDTI:", fdti_mean
    print "sum of terms:", np.sum(decomp_means.values())

    # plt.legend(calc.root_paths)
    plt.savefig(args.out)
