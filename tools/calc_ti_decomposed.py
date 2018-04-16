from collections import defaultdict
from itertools import chain
import matplotlib
import numpy as np
import os
import free_energy_base as feb
from table import Table

if "DISPLAY" not in os.environ or os.environ["DISPLAY"] == "":
    matplotlib.use('Agg')
import matplotlib.pyplot as plt


class TI_decomposed(feb.Estimator):
    """Estimate free differences using ThermodynamicIntegration. The class
    performs individual calculations for each component of the system
    interaction and internal energies. Results from each component
    will sum to give the same result as the total energy. The
    decomposition of the total free energy is not well-defined however
    can be indicative of the dominant terms in a free energy change.
    """
    def __init__(self, lambdas, results_name='results',
                 subdir_glob='./', **kwargs):
        self.estimators = defaultdict(
            lambda: feb.TI(lambdas, results_name, subdir_glob))
        self.lambdas = lambdas
        self.subdir_glob = subdir_glob
        self.results_name = results_name

    def add_data(self, series):
        """Save data from a SnapshotResults.series object
        for later calculation"""
        dlam = series.lamf[0]-series.lamb[0]
        min_len = 10E20
        for name, term in series.interaction_energies.iteritems():
            for energy in term[:-1]:
                dat = (energy.forw-energy.back) / dlam
                min_len = len(dat) if len(dat) < min_len else min_len
                self.estimators["%s_%s" % (name, energy.type)].data.append(dat)

        for name, term in series.internal_energies.iteritems():
            for energy in term[:-1]:
                dat = (energy.forw-energy.back) / dlam
                min_len = len(dat) if len(dat) < min_len else min_len
                self.estimators["%s_%s" % (name, energy.type)].data.append(dat)

        return min_len

    def calculate(self, subset=(0., 1., 1)):
        """Calculate the free energy difference for each energy component and
        return a dictionary of PMF objects."""
        return {term: self.estimators[term].calculate(subset)
                for term in self.estimators}

    def __getitem__(self, val):
        new_est = self.__class__(self.lambdas)
        for term in self.estimators:
            new_est.estimators[term] = self.estimators[term][val]
        return new_est

    def __len__(self):
        lengths = [dat.shape[-1] for term in self.estimators
                   for dat in self.estimators[term].data]
        if len(lengths) == 0:
            return 0
        if len(set(lengths)) > 1:
            raise Exception("Found data entries of different lengths.")
        return lengths[0]


class DecomposedCalculation(feb.FreeEnergyCalculation):
    def __init__(self, root_paths, subdir=''):
        # temperature can be given as zero as TI estimators do
        # not make use of it
        feb.FreeEnergyCalculation.__init__(
            self,
            root_paths=root_paths,
            temperature=0.,
            subdir=subdir,
            estimators={feb.TI, TI_decomposed})

    def calculate(self, subset=(0., 1., 1.)):
        results = {}

        leg_result = feb.Result()
        for leg in self.estimators[feb.TI]:
            leg_result += feb.Result(
                [rep.subset(*subset).calculate(self.temperature)
                 for rep in leg])
        results[feb.TI] = leg_result

        # dealing with the decomposed estimator is a little more difficult
        # repeats for the same term need to be grouped together into a
        # single Results object
        # start by pulling out all the data
        data = [[rep.subset(*subset).calculate(self.temperature)
                 for rep in leg]
                for leg in self.estimators[TI_decomposed]]

        # now add into results in a compatible order and combine legs such that
        # results[TI_decomposed][term] = Results object
        results[TI_decomposed] = {}
        for term in data[0][0]:
            result = feb.Result()
            for leg in data:
                result += feb.Result([rep[term] for rep in leg])
            results[TI_decomposed][term] = result

        return results

    def _body(self, args):
        subset = (args.lower_bound, args.upper_bound)
        results = self.calculate(subset=subset)
        decomp = results[TI_decomposed]

        if (args.bound or args.gas) is not None:
            calc2 = DecomposedCalculation(args.bound or args.gas)
            results2 = calc2.calculate(subset=subset)
            decomp2 = results2[TI_decomposed]

            # iterate over the unique keys from both calc and calc2
            for term in set(chain(decomp, decomp2)):
                try:
                    decomp[term] -= decomp2[term]
                except KeyError:
                    try:
                        decomp[term] = -decomp2[term]
                    except KeyError:
                        # term is not in decomp2, no furthor action needed
                        pass
                if args.bound is not None:
                    decomp[term] = -decomp[term]

            # update standard TI result as well
            results[feb.TI] -= results2[feb.TI]
            # swap the sign if this is a binding calculation
            if args.bound is not None:
                results[feb.TI] = -results[feb.TI]

        if args.dualtopology:
            consolidate_terms(decomp)
        self.figures['decomposed'] = plot_terms(decomp)

        if args.pmf:
            self.figures['decomposed_pmfs'] = plot_pmfs(decomp)

        table = Table('', fmts=["%s:", "%.3f"])
        table.add_row(["FDTI", results[feb.TI].dG])
        table.add_blank_row()
        for term in sorted(decomp):
            table.add_row([term, decomp[term].dG])
        table.add_row(['sum of terms', np.sum(decomp.values()).dG])
        self.tables.append(table)

        return results


def consolidate_terms(data):
    # this is a little better than the previous iteration, copes with
    # gcsolutes at least, set operations are still not ideal though
    components = {'protein1', 'solvent', 'GCS'}
    for inter in ('COU', 'LJ'):
        for comp in components:
            interaction_terms = []
            for term in data:
                if inter not in term or comp not in term:
                    continue

                # get the one or more parts of the energy term
                parts = term.split('_')[0].split('-')
                # if it has two parts and exactly one is a non-solute
                # component then we can consolidate
                if len(parts) == 2 and len(set(parts) - components) == 1:
                    interaction_terms.append(term)

            if len(interaction_terms) > 1:
                data['lig-%s_%s' % (comp, inter)] = np.sum(
                    [data.pop(term) for term in interaction_terms])


def plot_terms(data):
    fig, ax = plt.subplots()
    fig.subplots_adjust(bottom=0.3)

    ys, errs, labels = [], [], []
    for term in sorted(data):
        dG = data[term].dG
        # only include terms with non-zero differences
        if dG.value != 0.:
            ys.append(dG.value)
            errs.append(dG.error)
            labels.append(term)

    xs = np.arange(len(ys))
    ax.bar(xs, ys, yerr=errs, ecolor='black')
    ax.set_xticks(xs)
    ax.set_xticklabels(labels, rotation="vertical")
    ax.set_ylabel('Free Energy (kcal/mol)')
    return fig


def plot_pmfs(data):
    fig, ax = plt.subplots()
    for term in data:
        if data[term].dG.value != 0.:
            data[term].pmf.plot(ax, label=term)
    ax.legend(loc='best')
    ax.set_xlabel('Lambda Value')
    ax.set_ylabel('Free Energy (kcal/mol)')
    return fig


def get_arg_parser():
    """Add custom options for this script"""
    parser = feb.FEArgumentParser(
        description="Calculate individual contributions of different terms "
                    "to the total free energy difference. Although terms are "
                    "guaranteed to be additive with TI, the decomposition "
                    "is not strictly well defined. That said, it can be "
                    "illustrative to consider the dominant contributions of "
                    "a calculation.",
        parents=[feb.get_arg_parser()])
    parser.add_argument(
        "-b", "--bound", nargs="+", action='append',
        help="Output directory(s) of additional bound phase calculation(s). "
             "Using this flag causes data loaded via -d to be considered as "
             "solvent phase data. All data is then combined to provide a "
             "decomposition of the binding free energy. Behaves identically "
             "to -d in treatment of repeats and calculation legs.",
        clashes=['gas'])
    parser.add_argument(
        "-g", "--gas", nargs="+", action='append',
        help="As -b except data loaded via this flag is treated as gas phase "
             "data to provide to provide a decomposed solvation free energy.",
        clashes=['bound'])
    parser.add_argument(
        "--dualtopology", action='store_true', default=False,
        help="Indicates provided data is from a dual topology calculation. "
             "Attempts to consolidate terms, for clarity, from ligands that "
             "can have opposite signs and large magnitudes. Please note that "
             "standard errors calculated with this approach are no longer "
             "rigorous and can be spuriously large.")
    parser.add_argument("--pmf", action='store_true', default=False,
                        help="Plot the Potential of Mean Force for all terms.")
    return parser


if __name__ == "__main__":
    args = get_arg_parser().parse_args()
    calc = DecomposedCalculation(args.directories, subdir=args.subdir)
    calc.run(args)
