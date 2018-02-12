import matplotlib
import numpy as np
import os
import free_energy_base as feb
import calc_gci as gci
from table import Table

if "DISPLAY" not in os.environ or os.environ["DISPLAY"] == "":
    matplotlib.use('Agg')
import matplotlib.pyplot as plt


class GCIPMF(feb.Series):
    def __init__(self, coordinate, values, model, pmf):
        self.coordinate = coordinate
        self.values = values
        self.model = model
        self.pmf = feb.PMF(range(len(pmf)), pmf)


class GCIResult(feb.BaseResult):
    def __init(self, *args):
        if len(args) > 1:
            raise TypeError(
                'GCIResult objects do not support multiple data instances')
        feb.BaseResult.__init__(self, *args)

    @property
    def B_values(self):
        return feb.Result.lambdas.fget(self)

    @property
    def occupancies(self):
        return feb.Series(self.B_values, *self.data[0])

    @property
    def insertion_pmf(self):
        return feb.PMF(self.data[0][0].pmf.coordinate,
                       *[dat.pmf for dat in self.data[0]])

    @property
    def model(self):
        model_ys = []
        for rep in self.data[0]:
            rep.model.x = np.linspace(
                min(rep.coordinate), max(rep.coordinate), 100)
            rep.model.forward()
            model_ys.append(rep.model.predicted)
        return feb.Series(rep.model.x, *model_ys)


class GCI(feb.Estimator):
    def __init__(self, B_values):
        self.B_values = B_values
        self.data = []

    def add_data(self, series):
        self.data.append(series.solventson)
        return self.data[-1].shape[-1]

    def calculate(self, temp, volume, steps=None):
        Ns = np.array(self.data).mean(axis=1)
        if steps is None:
            steps = int(max(Ns))
            if max(Ns) - steps > 0.9:
                steps += 1

        model = gci.fit_ensemble(x=np.array(self.B_values), y=Ns, size=steps,
                                 verbose=False)[0]

        return GCIPMF(
            self.B_values,
            Ns,
            model,
            gci.insertion_pmf(np.arange(steps+1), model, volume)
        )

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
    def __init__(self, root_paths, temperature, volume, steps=None, subdir=''):
        feb.FreeEnergyCalculation.__init__(
            self,
            root_paths=root_paths,
            temperature=temperature,
            estimators=[GCI],
            subdir=subdir)
        self.volume = volume
        self.steps = steps

    def _path_constructor(self, root_path):
        return os.path.join(root_path, "b_*", self.subdir)

    def _get_lambda(self, path):
        return float(path.split('/')[-2].split('_')[1])

    def calculate(self, subset=(0., 1., 1)):
        """For each estimator return the evaluated potential of mean force.

        Returns
        -------
        dict:
          results of calculation stored in the below structure:
          return_val[estimator_class][i] = PMF instance
          where i is an list index indicating a repeat
        """
        results = {}
        for est, legs in self.estimators.items():
            leg_result = GCIResult()
            for i, leg in enumerate(legs):
                leg_result += GCIResult(
                    [rep.subset(*subset).calculate(self.temperature,
                                                   self.volume,
                                                   self.steps)
                     for rep in leg])
            results[est] = leg_result
        return results[GCI]

    def _body(self, args):
        results = self.calculate(subset=(args.lower_bound, args.upper_bound))

        fig, ax = plt.subplots()
        plot_titration(results, ax)
        self.figures['titration'] = fig

        fig, table = insertion_pmf(results)
        self.tables.append(table)
        self.figures['insertion_pmf'] = fig

        return results


def plot_titration(results, ax, dot_fmt='b'):
    for rep in results.data[0]:
        ax.plot(rep.coordinate, rep.values, 'o')
    results.model.plot(ax, xlabel='B Value', ylabel='Occupancy', color='black')


def insertion_pmf(results, title=''):
    table = Table(title, fmts=['%d', '%.3f'],
                  headers=['Number of Waters', 'Binding Free Energy'])

    steps = results.data[0][0].pmf.coordinate
    pmf = feb.PMF(steps, *[rep.pmf for rep in results.data[0]])
    for i, fe in enumerate(pmf):
        table.add_row([i, fe])

    fig, ax = plt.subplots()
    pmf.plot(ax, xlabel="Occupancy")
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
    parser.add_argument(
        '-v', '--volume', required=True, type=float,
        help="Volume of the calculations GCMC region.")
    parser.add_argument(
        '-n', '--nsteps', type=int,
        help='Override automatic guessing of maximum number of waters for '
             'titration curve fitting.')
    return parser


if __name__ == '__main__':
    args = get_arg_parser().parse_args()
    tc = TitrationCalculation(
        [args.directories], args.temperature, args.volume,
        args.nsteps, subdir=args.subdir)
    tc.run(args)
    plt.show()
