import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from operator import add
import pymbar
from scipy import optimize
from scipy import integrate
import free_energy_base as feb

thermal_wavelength = 1.00778365325  # of water, in angstroms
kB = 0.0019872  # Boltzmann's constant, in kcal.mol^-1.K^-1
standard_volume = 30.  # Standard volume of a water in bulk water
tip4p_excess = -6.2  # Excess chemical potential of tip4p in bulk tip4p


class GCMCPMF(feb.Series):
    def __init__(self, coordinate, values, volume, temperature, model, pmf):
        self.coordinate = coordinate
        self.values = values
        self.volume = volume
        self.temperature = temperature
        self.model = model
        self.pmf = feb.PMF(range(len(pmf)), pmf)

    @property
    def equilibrium_B(self):
        """Get the value of B equating to equilibrium with bulk water."""
        beta = 1 / (self.temperature * kB)
        return beta * tip4p_excess + np.log(self.volume / standard_volume)


class GCMCPMF2D(GCMCPMF):
    def __init__(self, lambdas, Bs, volume, temperature, *args):
        self.coordinate = lambdas
        self.coordinate2 = Bs
        self.volume = volume
        self.temperature = temperature
        self.values = []
        for dat in zip(*args):
            self.values.append(feb.PMF(self.coordinate2, *dat))

    @property
    def dG(self):
        eqb_B = self.equilibrium_B
        closest = np.argmin([abs(B - eqb_B) for B in self.coordinate2])
        return self.values[-1].values[closest] - self.values[0].values[closest]

    def plot(self, axes, show_error=True,  **kwargs):
        lambdas, Bs = np.meshgrid(self.coordinate2, self.coordinate)
        values = np.array([[col.value for col in row] for row in self.values])
        axes.plot_surface(lambdas, Bs, values, cstride=1, rstride=1)
        return axes


class GCMCResult(feb.BaseResult):
    def __init__(self, *args):
        if len(args) > 1:
            raise TypeError(
                'GCMCResult objects do not support multiple data instances')
        feb.BaseResult.__init__(self, *args)

    @property
    def B_values(self):
        return feb.Result.lambdas.fget(self)

    @property
    def equilibrium_B(self):
        return self.data[0][0].equilibrium_B

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


class GCMCResult2D(feb.Result):
    @property
    def pmf(self):
        pmfs = [GCMCPMF2D(self.lambdas, dat[0].coordinate2,
                          dat[0].volume, dat[0].temperature, *dat)
                for dat in self.data]
        return reduce(add, pmfs)


class GCI(feb.Estimator):
    def __init__(self, B_values, results_name='results', **kwargs):
        self.B_values = B_values
        self.data = []
        self.subdir_glob = ''
        self.results_name = results_name

    def add_data(self, series):
        self.data.append(series.solventson)
        return self.data[-1].shape[-1]

    def calculate(self, temp, volume, steps=None):
        Ns = np.array(self.data).mean(axis=1)
        if steps is None:
            steps = int(max(Ns))
            if max(Ns) - steps > 0.9:
                steps += 1

        model_steps = steps if steps != 0 else 1
        model = fit_ensemble(
            x=np.array(self.B_values), y=Ns, size=model_steps,
            verbose=False)[0]

        return GCMCPMF(self.B_values, Ns, volume, temp, model,
                       insertion_pmf(np.arange(steps + 1), model, volume))

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


# Functions specific to grand canonical integration and neural network fitting
def pseudohuber(r, c):
    """
    The pseudo Huber loss/cost function. Smooth version of Huber loss function.

    Parameters
    ----------
    r : float
      residual, i.e. difference between predicted value and its target
    c : float
      parameter that determines roughly where residuals will have less weight

    Return
    ------
    float
      The pseudo Huber cost for a residual r and parameter c

    """
    return (np.sqrt(1 + (r / c)**2) - 1) * c**2


def logistic(params, x):
    """
    Logistic function

    Parameters
    ----------
    params : list or numpy array
      the three parameters of the logistic function
    x : numpy array
      the explanatory variable

    Return
    ------
    numpy array
      the output of the logistic function

    """
    return params[2] / (1 + np.exp(-params[0] - x * params[1]))


def integrated_logistic(params, Bi, Bf):
    """
    The integral of the logistic function, for use in the function
    "insertion_pmf". The problem is that it is prone to numerical instabilities
    when the exponents are large and positive. While one could fix the problem,
    I've opted for numerical integration instead in "insertion_pmf".
    This function has been used to verify the numerical integration.

    Parameters
    ----------
    params : list or numpy array
      the parameters of the logistic function
    Bi : float
      the lower value of the interval of integration
    Bf : float
      the higher value of the interval of integration

    Return
    ------
    float
      the area under the logistic curve between Bi and Bf

    """
    return -(params[2]/params[1]) * \
        (np.log(1 + np.exp(Bf * params[1]+params[0])) -
         np.log(1 + np.exp(Bi * params[1]+params[0])))


def simplex_sample(size):
    """
    Samples uniformly over a simplex, such that the output sums to 1.

    Parameters
    ----------
    size : integer
      the number of random numbers you want

    Return
    ------
    numpy array
      vector of uniform random numbers that will sum to 1.

    """
    samp = np.random.uniform(0.0, 1.0, size - 1)
    samp = np.append([0.0, 1.0], samp)
    samp = np.sort(samp)
    return np.diff(samp)


class Slp(object):
    """
    Class for fitting and storing a monotonically increasing single layer
    perceptron with one input and one output.

    Fitting is achieved by iteratively optimising one unit at a time,
    rather than the more typical back-propigation algorithm, to make use of
    SciPy's monotonic optimisation routines.

    Attributes
    ----------
    size : integer
      the number logsistic terms, equilivalent to the number of
      'hidden units' in an artificial neural network
    x : numpy array
      the explanatory variable, which is the Adams value or
      chemical potential in GCMC
    y : numpy array
      the response variable, which is the average number of
      inserted grand canonical molecules
    weights : numpy array
      the matrix of coefficients of the logistic terms
    offset : numpy array
      the intercept of the artificial neural network. Grouped
      together with the weights when fitting
    streams : list
      contains the output of each logistic function or 'unit'
    predicted : numpy array
      the predicted output of the artificial neural network
    error: float
      the error (either mean-squared, Huber, or absolute loss) of the
      artificial neural network with respect to the response variable.
    """

    def __init__(self, x=np.array(0), y=np.array(0), size=1):
        self.size = size
        self.x = x
        self.y = y
        self.weights = np.ones((size, 3))
        self.offset = np.array(0.0)
        self.streams = [self.x] * (size)
        self.predicted = np.array(0)
        self.error = 1000.0

    def forward(self):
        """
        Evaluates the output from ALL units for a given set of weights
        and updates the prediction of the whole network
        """
        self.streams = [logistic(params, self.x) for params in self.weights]
        self.predicted = np.sum(np.vstack(self.streams), axis=0) + self.offset

    def forward_unit(self, unit):
        """
        Evaluates the output from ONE unit for a given set of weights
        and updates the prediction of the whole network

        Parameters
        ----------
        unit : the index of the logisitic function
        """
        self.streams[unit] = logistic(self.weights[unit], self.x)
        self.predicted = np.sum(np.vstack(self.streams), axis=0) + self.offset

    def predict(self, x):
        """
        Returns the prediction for the artificial neural network

        Parameters
        ----------
        x : the vector of Adams or chemical potentials
        """
        streams = [logistic(params, x) for params in self.weights]
        return np.sum(np.vstack(streams), axis=0) + self.offset

    def msd(self):
        """
        Calculates the mean-squared error of the model and updates the error
        """
        self.error = np.sum((self.y - self.predicted)**2)
        return self.error

    def huber(self, c=2):
        """
        Calculates the pseudo Huber loss function and
        updates the error of the model

        Parameters
        ----------
        c : number
          the pseudo Huber loss function parameter
        """
        self.error = np.sum(pseudohuber(self.y - self.predicted, c))
        return self.error

    def absolute(self):
        """
        Calculates the absolute loss of the model and updates the error
        """
        self.error = np.sum(np.abs(self.y - self.predicted))
        return self.error

    def randomise(self, pin_max=True):
        """
        Produces initial values of weights by random sampling

        Parameters
        ----------
        pin_max : whether to constrain the random weights so that the sum
          of the units cannot exceed the maximum value of the response variable
        """
        self.offset = self.y.min()
        alpha = np.random.uniform(low=0.0, high=2.0, size=self.size)
        alpha0 = np.sort(-np.random.uniform(
            low=self.x.min(), high=self.x.max(), size=self.size))
        if pin_max:
            # The beta weights are sampled uniformally from a simplex,
            # such that they add up to range of y values.
            beta = simplex_sample(self.size) * (
                self.y.max() - self.offset
            )
        else:
            beta = np.random.uniform(
                low=0, high=1, size=self.size) * self.y.max()
        beta.sort()
        self.weights = np.vstack(
            (alpha0, alpha,
             beta)).T  # Each row is an element in the logistic function.
        self.forward()  # Update the logistic streams

    def randsearch(self, samples):
        """
        Performs a crude initial random search of the space of weights,
        and saves the best set of weights

        Parameters
        ----------
        samples : how many initial guesses of the weights shall be attempted
        """
        error_old = 1000000
        for samps in range(samples):
            self.randomise()
            self.forward()
            error_new = self.msd()
            if error_new < error_old:
                weights = self.weights
                streams = self.streams
                self.forward()
                error_old = self.msd()
        self.weights = weights
        self.streams = streams
        self.forward()
        self.error = error_old

    def fit_unit(self,
                 unit,
                 grad_tol=0.001,
                 pin_max=None,
                 monotonic=True,
                 cost="msd",
                 c=2.0):
        """
        Fits an individual unit of the artificial neural network

        Parameters
        ----------
        unit : int
          the index of the logistic function whose weights
          (parameters) will be optimised
        grad_tol : float
          the tolorance of the gradient that determines when optimisation stops
        pin_max : float
          the value of the maximum value to model can produce.
          If None, it will be fitted
        monotonic : boolean
          whether to insure the artificial neural network is
          monotonically increasing
        cost : string
           the type of cost/loss function that will be used
           can be 'msd', 'abs' or 'huber'
        c : float
          the paramater of the Huber loss function
        """
        if cost == "msd":
            def residuals(params):
                self.weights[unit] = params
                self.forward_unit(unit)
                return (self.msd())

        if cost == "abs":
            def residuals(params):
                self.weights[unit] = params
                self.forward_unit(unit)
                return (self.absolute())

        if cost == "huber":
            def residuals(params):
                self.weights[unit] = params
                self.forward_unit(unit)
                return (self.huber(c))

        if type(pin_max) == float or type(pin_max) == int:
            betas = np.ma.array(self.weights[:, 2], mask=False)
            betas.mask[unit] = True
            unit_max = pin_max - self.offset - np.sum(betas)
        else:
            unit_max = None
        if monotonic:
            fit = optimize.fmin_l_bfgs_b(
                residuals,
                x0=self.weights[unit],
                bounds=[(-10.0, unit_max), (0.0, unit_max), (self.offset,
                                                             unit_max)],
                approx_grad=True,
                pgtol=grad_tol,
                disp=0)
            self.weights[unit] = fit[0]
            self.error = fit[1]
        else:
            fit = optimize.minimize(
                residuals,
                self.weights[unit],
                method='BFGS',
                options={'gtol': grad_tol})
            self.weights[unit] = fit.x
            self.error = fit.fun

    def fit_offset(self, grad_tol=0.001, cost="msd", c=2.0):
        """
        Fits the intercept of the artificial neural network

        Parameters
        ----------
        grad_tol : float
          the tolorance of the gradient that determines when optimisation stops
        cost : string
          the type of cost/loss function that will be used
           can be 'msd', 'abs' or 'huber'
        c : float
          the paramater of the Huber loss function
        """
        diff = self.y - np.sum(np.vstack(self.streams), axis=0)
        if cost == "msd":
            def residuals(offset):
                return np.sum((diff - offset)**2)

        if cost == "abs":
            def residuals(offset):
                return np.sum(np.abs(diff - offset))

        if cost == "huber":
            def residuals(offset):
                return np.sum(pseudohuber(diff - offset, c))

        fit = optimize.fmin_l_bfgs_b(
            residuals,
            x0=np.array(self.offset),
            bounds=[(0.0, None)],
            approx_grad=True,
            pgtol=grad_tol,
            disp=0)
        self.offset = fit[0]
        self.error = fit[1]

    def epoch(self,
              grad_tol=0.001,
              pin_min=None,
              pin_max=None,
              monotonic=True,
              cost="msd",
              c=2.0):
        """
        Fits each unit and intercept of the artificial neural network once

        Parameters
        ----------
        grad_tol : float
          the tolorance of the gradient that determines when optimisation stops
        pin_min : float
          the value of the fixed intercept. If None, it will be fitted.
        pin_max : float
          the value of the maximum value to model can produce.
          If None, it will be fitted
        monotonic : boolean
          whether to insure the artificial neural
          network is monotonically increasing
        cost : string
          the type of cost/loss function that will be used
          can be 'msd', 'abs' or 'huber'
        c : float
          the paramater of the huber loss function
        """
        if pin_min is None:
            self.fit_offset(grad_tol, cost, c)
        if type(pin_min) == float or type(pin_min) == int:
            self.offset = pin_min
        for unit in range(self.size):
            self.fit_unit(unit, grad_tol, pin_max, monotonic, cost, c)

    def train(self,
              iterations=100,
              grad_tol_low=-3,
              grad_tol_high=1,
              pin_min=None,
              pin_max=None,
              monotonic=True,
              cost="msd",
              c=2.0):
        """
        Iterates over a specified amount of training epochs. After each epoch,
        the fitting tolerance will decrease to gradually refine the fits.

        Parameters
        ----------
        iterations : int
          the maximum number of training epochs that will be tested.
        grad_tol_low : int
          ten to the power of this value is the
          initial tolerance of the gradient
        grad_tol_high : int
          ten to the power of this value is the final tolerance of the gradient
        pin_min : float
          the value of the fixed intercept. If left blank, it will be fitted.
        pin_max : float
          the value of the maximum value to model can produce.
          If None, it will be fitted
        monotonic : boolean
          whether to insure the artificial neural network
          is monotonically increasing
        cost : string
          the type of cost/loss function that will be used
          can be 'msd', 'abs' or 'huber'
        c : float
          the paramater of the huber loss function
        """
        tolerances = np.logspace(grad_tol_high, grad_tol_low, iterations)
        for grad_tol in tolerances:
            self.epoch(grad_tol, pin_min, pin_max, monotonic, cost, c)
        self.forward()


def inverse_slp(model, y):
    """
    For a given ANN model, this function returns the value of x that
    corresponds to the value of y. It's basically the inverse function of an
    ANN. Used to determine which B value produces a given average of waters.

    Parameters
    ----------
    model : Slp object
      the ANN (single layer perceptron) to return the inverse function of
    y : float
      the value which you want to find the corresponding value of x

    Return
    ------
    float
      the value of x that produced y in the input ANN
    """
    dev = (model.predicted - y)**2
    x_guess = model.x[np.where(dev == dev.min())][0]

    def errfunc(x):
        return (y - model.predict(x))**2

    solution = optimize.minimize(errfunc, x0=x_guess, method="BFGS")
    return solution.x


def fit_ensemble(x,
                 y,
                 size,
                 repeats=20,
                 randstarts=1000,
                 iterations=100,
                 grad_tol_low=-3,
                 grad_tol_high=1,
                 pin_min=False,
                 pin_max=None,
                 monotonic=True,
                 cost="msd",
                 c=2.0,
                 verbose=True):
    """
    Fits a collection of ANN models to GCMC titration data. Each ANN
    is fitted with different set of initial ANN weights.

    Parameters
    ----------
    x : numpy array
      the B/Adams values or chemical potentials used in a titration
    y : numpy array
      the average value of water at each chemical potential
    repeats : int
      the number of repeats of the fitting. This gives the "ensemble" of models
    others
      As specified in Slp.

    Return
    ------
    list
      a list of Slp objects, each has been fit with different initial weights.
    """
    # A list that contains the fitted single layer perceptrons.
    models = []
    # Initialising the lowest error. The best model that with the lowest error
    error = 1E5
    for i in range(repeats):
        model = Slp(x=x, y=y, size=size)
        model.randsearch(randstarts)
        model.train(
            iterations,
            grad_tol_low=grad_tol_low,
            grad_tol_high=grad_tol_high,
            pin_min=pin_min,
            pin_max=pin_max,
            cost=cost,
            c=c,
            monotonic=True)
        if model.error < error:
            single_model = model
        models.append(model)
        if verbose:
            print "Model %i fitted with error = %f" % (i + 1, model.error)
    return single_model, models


def insertion_pmf(N, gcmc_model, volume, T=298.15):
    """
    Calculates the free energy to insert water from ideal gas to
    the GCMC volume.

    Parameters
    ----------
    N : numpy array
      vector of bound water numbers between which relative
      free energies will be calculated
    gcmc_model : Slp object
      an ANN model (single layer perceptron) that will
      be used to calculate the free energy
    volume : float
      the volume of the GCMC region in Angstroms^3
    T : float
      the temperature of the simulation in Kelvin

    Return
    ------
    numpy array
      relative ideal gas transfer free energies for each water specified in N
    """
    # Boltmann's constant (in kcal/mol/K) multiplied by temperature (in K).
    kT = T * kB

    def integral(gcmc_model, Bi, Bf):  # For numerical integration.
        logies = [
            integrated_logistic(params, Bi, Bf)
            for params in gcmc_model.weights
        ]
        return np.sum(logies)  # Under the convention that 0*log(0) = 0.

    def nonintegral(N, B):
        return N * B

    dG = np.zeros(N.size)
    correction = np.zeros(N.size)
    N_lower = N[0]
    B_lower = inverse_slp(gcmc_model, N_lower)
    for i in range(
            1,
            N.size):  # Starting from the second element, as dG=0 at first N.
        B_upper = inverse_slp(gcmc_model, N[i])
        dG[i] = (nonintegral(N[i], B_upper) - nonintegral(N_lower, B_lower) -
                 integrate.quad(gcmc_model.predict, a=B_lower, b=B_upper)[0]
                 )[0]  # This version numerical integration.
        correction[i] = (N[i] - N_lower) * np.log(volume / standard_volume)
        dG[i] = dG[i] - correction[i]
    return dG * kT


class GCMCMBAR(feb.MBAR):
    result_class = GCMCResult2D

    def __init__(self, lambdas, volume, **kwargs):
        self.lambdas = lambdas
        self.volume = volume
        self.data = []
        self._data_lambdas = []
        self.subdir_glob = 'b_*/'
        self._data_Bs = []
        self.results_name = 'results_inst'

        # deals with a gotcha when user forgets the volume flag
        # when using analysis script interface
        if volume is None:
            raise TypeError(
                'No gcmc volume information passed to GCMCMBAR estimator. '
                'If using a script you may need the --volume flag.')

    def add_data(self, series):
        dat = np.array([series.feenergies[lam]
                        for lam in sorted(series.feenergies)])
        N = dat.shape[-1]
        # add solventson data as the final row of each data entry
        # not an ideal solution but saves the need to rewrite __getitem__
        dat = np.append(dat, series.solventson.reshape((1, N)), axis=0)
        self.data.append(dat)
        self._data_lambdas.append(series.lam[0])
        self._data_Bs.append(series.bvalue)
        return N

    def B_to_chemical_potential(self, B, temperature):
        """Convert a B value to the equivalent chemical potential"""
        beta = 1 / (temperature * kB)
        return (B-np.log(self.volume/thermal_wavelength**3)) / beta

    @property
    def Bs(self):
        return sorted(set(self._data_Bs))

    def calculate(self, temp=300.):
        N = self.data[0].shape[-1]
        N_sims = len(self.data)
        N_Bs = len(self.Bs)
        N_lams = len(self.lambdas)
        beta = 1 / (kB * temp)

        u_kn = np.zeros(shape=(N_sims, N*N_sims))
        ns = np.zeros(N*N_sims)
        for dat, B_sim, lam_sim in zip(
                self.data, self._data_Bs, self._data_lambdas):
            # due to globbing simulation data will not have been
            # read in a logical order, need the below indices to
            # determine the proper position for each set of data
            i_B = self.Bs.index(B_sim)
            i_lam = self.lambdas.index(lam_sim)
            ns = dat[-1]  # water occupancy  data
            es = dat[:-1]  # simulation energies
            for i in xrange(N_lams):
                for j, B in enumerate(self.Bs):
                    mu = self.B_to_chemical_potential(B, temp)
                    u_kn[i*N_Bs+j, (i_lam*N_Bs+i_B)*N: (i_lam*N_Bs+i_B+1)*N] =\
                        beta*(es[i] + ns * mu)
        samples = np.array([N] * N_sims)
        mbar = pymbar.MBAR(u_kn, samples)
        free_energies = mbar.getFreeEnergyDifferences(
            compute_uncertainty=False, return_theta=False)[0] / beta

        # free energy order seems to be reversed with respect to B values
        return GCMCPMF2D(self.lambdas, self.Bs[::-1], self.volume, temp,
                         free_energies[0].reshape([N_lams, N_Bs]))


def plot_surface(lambdas, Bs, Z):
    fig = plt.figure()
    ax = Axes3D(fig)
    lambdas, Bs = np.meshgrid(Bs, lambdas)

    ax.plot_surface(lambdas, Bs, Z,
                    cstride=1, rstride=1)
    return fig
