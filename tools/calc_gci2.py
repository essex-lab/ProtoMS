import matplotlib
import numpy as np
import os
import free_energy_base as feb
import calc_gci as gci

if "DISPLAY" not in os.environ or os.environ["DISPLAY"] == "":
    matplotlib.use('Agg')
import matplotlib.pyplot as plt


# this redundancy needs sorting out, preferably by generalising the
# pmf object to a more series
class GCIResult():
    def __init__(self, lambdas, values, model, pmf):
        self.lambdas = lambdas
        self.values = values
        self.model = model
        self.pmf = feb.PMF(range(len(pmf)), pmf)

    @property
    def dG(self):
        """Return the free energy difference at the PMF end points."""
        return self.values[-1] - self.values[0]

    def __neg__(self):
        return PMF(self.lambdas, [-val for val in self.values])

    def __iter__(self):
        return iter(self.values)


class GCI(feb.Estimator):
    def __init__(self, B_values):
        self.B_values = B_values
        self.data = []

    def add_data(self, series):
        self.data.append(series.solventson)
        return self.data[-1].shape[-1]

    def calculate(self, subset=(0., 1., 1)):
        # model = linear_model.Perceptron()
        # print np.array(self.data).mean(axis=1)
        Ns = np.array(self.data).mean(axis=1)
        model = gci.fit_ensemble(x=np.array(self.B_values), y=Ns, size=2,
                                 verbose=False)[0]

        return GCIResult(
            self.B_values,
            Ns,
            model,
            gci.insertion_pmf(np.array([0, 1, 2]), model, 30.)
        )
        # model.fit(np.array(self.B_values).reshape((len(self.B_values), 1)),
        #           )

    def __getitem__(self, val):
        """Return a new class instance with series[val] applied to each
        individual data series.
        """
        new_est = self.__class__(self.B_values)
        # add data series to the new estimator that have been sliced by val
        # want to always apply slice to last dimension, so transpose array
        # apply slice to first dimension and then transpose back
        for dat in self.data:
            reordered_dat = dat.T[val]
            new_est.data.append(reordered_dat.T)
        return new_est


class TitrationCalculation(feb.FreeEnergyCalculation):
    def __init__(self, root_paths, temperature, subdir=''):
        feb.FreeEnergyCalculation.__init__(
            self,
            root_paths=root_paths,
            temperature=temperature,
            estimators=[GCI],
            subdir=subdir)

    def _path_constructor(self, root_path):
        return os.path.join(root_path, "b_*", self.subdir)

    def _get_lambda(self, path):
        return float(path.split('/')[-2].split('_')[1])

    def calculate(self, *args, **kwargs):
        return feb.FreeEnergyCalculation.calculate(self, *args, **kwargs)[GCI]

    def _body(self, args):
        results = self.calculate(subset=(args.lower_bound, args.upper_bound))

        self.figures['titration'] = plot_titration(results)
        fig, table = insertion_pmf(results)
        self.tables.append(table)
        self.figures['insertion_pmf'] = fig

        return results


def plot_titration(results):
    fig, ax = plt.subplots()
    model_ys = []
    for rep in results.data[0]:
        ax.scatter(rep.lambdas, rep.values)
        rep.model.x = np.linspace(min(rep.lambdas), max(rep.lambdas), 100)
        rep.model.forward()
        model_ys.append(rep.model.predicted)

    ys = np.mean(model_ys, axis=0)
    errs = np.std(model_ys, axis=0) / len(results.data)**0.5

    ax.plot(rep.model.x, ys, 'r')
    ax.fill_between(rep.model.x, ys + 2*errs, ys - 2*errs,
                    facecolor='gray', alpha=0.3,
                    interpolate=True, linewidth=0.)
    ax.fill_between(rep.model.x, ys + errs, ys - errs, facecolor='yellow',
                    alpha=0.3, interpolate=True, linewidth=0.)
    return fig


def insertion_pmf(results):
    table = feb.Table('', fmts=['%d', '%.3f'],
                      headers=['Number of Waters', 'Binding Free Energy'])

    pmf = feb.PMF(*[rep.pmf for rep in results.data[0]])
    for i, fe in enumerate(pmf):
        table.add_row([i, fe])

    fig, ax = plt.subplots()
    pmf.plot(ax)
    return fig, table


def get_arg_parser():
    """Add custom options for this script"""
    parser = feb.FEArgumentParser(
        description="Calculate water binding free energies using Grand "
                    "Canonical Integration.",
        parents=[feb.get_arg_parser()], conflict_handler='resolve')
    parser.add_argument(
        '-d', '--directories', nargs='+', required=True,
        help="Location of folders containing ProtoMS output subdirectories. "
             "Multiple directories can be supplied to this flag and indicate "
             "repeats of the same calculation.")
    return parser


if __name__ == '__main__':
    args = get_arg_parser().parse_args()

    tc = TitrationCalculation(
        [args.directories], 300., subdir=args.subdir)
    tc.run(args)
    plt.show()
