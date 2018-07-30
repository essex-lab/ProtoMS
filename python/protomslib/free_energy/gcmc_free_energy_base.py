from __future__ import print_function
import matplotlib
import numpy as np
import os
import pymbar
from scipy import optimize
from scipy import integrate
from . import free_energy_base as feb
from .table import Table
from .free_energy_argument_parser import FEArgumentParser
if "DISPLAY" not in os.environ or os.environ["DISPLAY"] == "":
    matplotlib.use('Agg')
import matplotlib.pyplot as plt

thermal_wavelength = 1.00778365325  # of water, in angstroms
kB = 0.0019872  # Boltzmann's constant, in kcal.mol^-1.K^-1
standard_volume = 30.  # Standard volume of a water in bulk water
tip4p_excess = -6.2  # Excess chemical potential of tip4p in bulk tip4p


class GCMCPMF(feb.Series):
    """A potential of mean force object for Grand Canonical calculations.
    Stores occupancy data and PMF for insertion free energies.
    """
    def __init__(self, coordinate, values, volume, temperature, model, pmf):
        """Parameters
        ----------
        coordinate : sequence of numbers
          The values of the controlled coordinate for the PMF
        values : sequence of numbers
          The values of the variable coordinate for the PMF
        volume : number
          The volume of the GCMC region used in a calculation
        temperature : float
          The temperature used in a calculation
        model : Slp object
          A Slp object fitted to the simulation occupancy data
        pmf : PMF object
          A PMF containing insertion free energies
        """
        self.coordinate = coordinate
        self.values = values
        self.volume = volume
        self.temperature = temperature
        self.model = model
        self.pmf = pmf

    @property
    def equilibrium_B(self):
        """Returns the value of B equating to equilibrium with bulk water."""
        beta = 1 / (self.temperature * kB)
        return beta * tip4p_excess + np.log(self.volume / standard_volume)


class GCMCResult(feb.BaseResult):
    def __init__(self, *args):
        if len(args) > 1:
            raise TypeError(
                'GCMCResult objects do not support multiple data instances')
        feb.BaseResult.__init__(self, *args)

    @property
    def B_values(self):
        """Return the simulated B values"""
        return feb.Result.lambdas.fget(self)

    @property
    def equilibrium_B(self):
        """Return the B value representing equilibrium with bulk water"""
        return self.data[0][0].equilibrium_B

    @property
    def occupancies(self):
        """Return a Series object containing occupancy data"""
        return feb.Series(self.B_values, *self.data[0])

    @property
    def insertion_pmf(self):
        """Return a PMF object containing insertion free energies"""
        return feb.PMF(self.data[0][0].pmf.coordinate,
                       *[dat.pmf for dat in self.data[0]])

    @property
    def model(self):
        """Return a Series object containing an average fitted model"""
        model_ys = []
        for rep in self.data[0]:
            rep.model.x = np.linspace(
                min(rep.coordinate), max(rep.coordinate), 100)
            rep.model.forward()
            model_ys.append(rep.model.predicted)
        return feb.Series(rep.model.x, *model_ys)


class GCI(feb.Estimator):
    """Estimate insertion free energies using Grand Canonical Integration."""
    def __init__(self, B_values, volume, results_name='results',
                 nsteps=None, nmin=None, nmax=None, nfits=10, **kwargs):
        """Parameters
        ----------
        B_values: list of floats
          The B values used in a calculation
        volume:
          The volume of the GCMC region used in a calculation
        results_name: string, optional
          Filename of the ProtoMS results file to use
        nsteps: int, optional
          The number of steps to fit in the titration curve.
          By default the number of steps is automatically selected.
        nfits: int, optional
          The number of independent fitting attempts for the neural network
          occupancy model. Increasing the number of fits may help improve
          results for noisy data.
        **kwargs:
          Additional keyword arguments are ignored
        """
        self.B_values = B_values
        self.data = []
        self.subdir_glob = ''
        self.results_name = results_name
        self.volume = volume
        self.nsteps = nsteps
        self.nmin = nmin
        self.nmax = nmax
        self.nfits = nfits

    def add_data(self, series):
        """Save data from a SnapshotResults.series object for later
        calculation. Return the length of the data series added.

        Parameters
        ----------
        series: SnapshotResults.series
            free energy simulation data for a particular B value
        """
        self.data.append(series.solventson)
        return self.data[-1].shape[-1]

    @property
    def N_min(self):
        Ns = np.array(self.data).mean(axis=1)
        return self.nmin if self.nmin is not None else int(round(min(Ns)))

    @property
    def N_max(self):
        Ns = np.array(self.data).mean(axis=1)
        return self.nmax if self.nmax is not None else int(round(max(Ns)))

    def calculate(self, temp):
        """Calculate the free energy difference and return a GCMCPMF object.

        Parameters
        ----------
        temp: float
          Temperature of calculation
        """
        Ns = np.array(self.data).mean(axis=1)

        if self.nsteps is None:
            if self.nmin is not None or self.nmax is not None:
                raise ValueError(
                    'If manually specifying the maximum or minimum number of'
                    ' waters, nsteps must be manually specified as well.')
            else:
                nsteps = self.N_max - self.N_min
        else:
            nsteps = self.nsteps

        model_steps = nsteps if nsteps != 0 else 1
        model = fit_ensemble(
            x=np.array(self.B_values), y=Ns, size=model_steps,
            repeats=self.nfits, verbose=False)[0]
        steps = np.arange(self.N_min, self.N_max+1)
        pmf = feb.PMF(steps, insertion_pmf(steps, model, self.volume))
        return GCMCPMF(self.B_values, Ns, self.volume, temp, model, pmf)


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
            print("Model %i fitted with error = %f" % (i + 1, model.error))
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
    """Estimate free energy differences using the Multistate Bennett Acceptance
    ratio for alchemical simulations performed in the Grand Canonical ensemble
    at multiple B values.
    """
    def __init__(self, lambdas, volume, **kwargs):
        """Parameters
        ----------
        lambdas: list of floats
          lambda values used in a calculation
        volume: float
          the volume of the gcmc region in a calculation
        """
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
        """Save data from a SnapshotResults.series object for later
        calculation. Return the length of the data series added.

        Parameters
        ----------
        series: SnapshotResults.series
            free energy simulation data for a particular simulation
        """
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
        """Convert a B value to the equivalent chemical potential for this
        simulation setup. Depends on GCMC volume of the simulation.

        Parameters
        ----------
        B: float
          the input B value
        temperature: float
          the simulated temperature
        """
        beta = 1 / (temperature * kB)
        return (B-np.log(self.volume/thermal_wavelength**3)) / beta

    @property
    def Bs(self):
        """Return the B value used for this calculation."""
        return sorted(set(self._data_Bs))

    def equilibrium_B(self, temperature):
        """Returns the value of B equating to equilibrium with bulk water."""
        beta = 1 / (temperature * kB)
        return beta * tip4p_excess + np.log(self.volume / standard_volume)

    def calculate(self, temp=300.):
        """Calculate the free energy difference and return a PMF object.

        Parameters
        ----------
        temp: float, optional
              temperature of calculation, default 300 K
        """
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
            for i in range(N_lams):
                for j, B in enumerate(self.Bs):
                    mu = self.B_to_chemical_potential(B, temp)
                    u_kn[i*N_Bs+j, (i_lam*N_Bs+i_B)*N: (i_lam*N_Bs+i_B+1)*N] =\
                        beta*(es[i] + ns * mu)
        samples = np.array([N] * N_sims)
        mbar = pymbar.MBAR(u_kn, samples)
        free_energies = mbar.getFreeEnergyDifferences(
            compute_uncertainty=False, return_theta=False)[0] / beta

        # free energy order seems to be reversed with respect to B values
        closest = np.argmin([abs(B - self.equilibrium_B(temp))
                             for B in self.Bs[::-1]])

        return feb.PMF(
            self.lambdas, free_energies[0].reshape([N_lams, N_Bs])[:, closest])


class TitrationCalculation(feb.FreeEnergyCalculation):
    """Calculate free energy differences using the Grand Canonical Integration
    method. Collates water occupancy data from simulations carried out at
    diferent chemical potentials and gives NVT binding free energies."""
    def __init__(self, root_paths, temperature, volume, nsteps=None,
                 nmin=None, nmax=None, nfits=5, pin_min=None,**kwargs):
        """Parameters
        ---------
        root_paths: a list of strings
          Paths to ProtoMS output directories. Each string specifies an
          the output directory containing a different repeat of a
          calculation. Reported results are averaged over the given repeats.
        temperature: float
          Simulation temperature in degrees Kelvin
        volume: float
          The volume of the simulated GCMC region
        nsteps: integer, optional
          Number of steps to use in fitting the GCI neural network. The
          default value of None attempts to guess the required number of steps.
        nmin: integer, optional
          Lowest number of waters to include when reporting network binding
          free energies
        nmin: integer, optional
          Highest number of waters to include when reporting network binding
          free energies
        nfits: int, optional
          The number of independent fitting attempts for the neural network
          occupancy model. Increasing the number of fits may help improve
          results for noisy data.
        **kwargs:
            Additional keyword arguments are passed to Estimator classes
            during initialisation.
        """
        self.subdir = ''
        feb.FreeEnergyCalculation.__init__(
            self,
            root_paths=[root_paths],
            temperature=temperature,
            estimators=[GCI],
            volume=volume,
            nsteps=nsteps,
            nmin=nmin,
            nmax=nmax,
            nfits=nfits,
	    pin_min=pin_min,
            **kwargs)

    def _path_constructor(self, root_path):
        """Given a root_path (string) construct a full path
        suitable for globbing to find output directories."""
        return os.path.join(root_path, "b_*", self.subdir)

    def _get_lambda(self, path):
        """Given a path (string) extract the contained lambda value"""
        return float(path.split('/')[-2].split('_')[1])

    def calculate(self, subset=(0., 1., 1)):
        """For each estimator return the evaluated potential of mean force.

        Parameters
        ----------
        subset: tuple of two floats and an int
            specify the subset of data to use in the calculation
        """
        results = {}
        for est, legs in self.estimators.items():
            leg_result = GCMCResult()
            for i, leg in enumerate(legs):
                tmp = GCMCResult([
                    rep.subset(*subset).calculate(self.temperature)
                    for rep in leg
                ])
                leg_result += tmp
            results[est] = leg_result
        return results[GCI]

    def _body(self, args):
        """The business logic of the calculation.

        Parameters
        ----------
        args: argparse.Namespace object
            Namespace from argumentparser
        """

        if args.nmin is not None or args.nmax is not None:
            self.header += '\nWARNING! Manually setting the maximum or ' \
                           'minimum number of waters. For free energies to ' \
                           'be meaningful ensure that observed simulation ' \
                           'occupancies cover a sufficient range.\n'
        results = self.calculate(subset=(args.lower_bound, args.upper_bound))

        fig, ax = plt.subplots()
        plot_titration(results, ax)
        self.figures['titration'] = fig

        fig, table = plot_insertion_pmf(results)
        self.tables.append(table)
        self.figures['insertion_pmf'] = fig

        eqb_B = results.equilibrium_B
        eqb_index = np.argmin(
            [abs(B - eqb_B) for B in results.occupancies.coordinate])
        closest_B = results.occupancies.coordinate[eqb_index]

        self.footer += 'The equilibrium B value is %.3f\n' % eqb_B
        self.footer += 'Most similar simulated B value is %.3f\n' % closest_B
        self.footer += 'Occupancy at %.3f is %s\n' % \
                       (closest_B, results.occupancies.values[eqb_index])

        min_index = np.argmin([
            v.value - i*tip4p_excess
            for i, v in enumerate(results.insertion_pmf.values)])
        self.footer += "\nOccupancy at binding PMF minimum is %.3f\n" % \
                       results.insertion_pmf.coordinate[min_index]
        return results


def plot_titration(results, ax, dot_fmt='b'):
    """Convenience function to plot the titration data from repeats."""
    for rep in results.data[0]:
        ax.plot(rep.coordinate, rep.values, 'o')
    results.model.plot(ax, xlabel='B Value', ylabel='Occupancy', color='black')


def plot_insertion_pmf(results, title=''):
    """Convenience function to plot the insertion pmf data from repeats"""
    table = Table(
        title,
        fmts=['%d', '%.3f', '%.3f', '%.3f'],
        headers=['Number of Waters',
                 'Insertion Free Energy',
                 'Network Binding Free Energy',
                 'Water Binding Free Energy'])

    steps = results.data[0][0].pmf.coordinate
    pmf = feb.PMF(steps, *[rep.pmf for rep in results.data[0]])
    prev_fe = 0.0
    for i, (dA, step) in enumerate(zip(pmf, steps)):
        bind_fe = dA - i*tip4p_excess
        table.add_row([step, dA, bind_fe, bind_fe - prev_fe])
        prev_fe = bind_fe.value

    fig, ax = plt.subplots()
    pmf.plot(ax, xlabel="Occupancy")
    return fig, table


def get_gci_arg_parser():
    parser = feb.FEArgumentParser(
        parents=[feb.get_base_arg_parser()],
        conflict_handler='resolve')
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
        help='Override automatic guessing of the number of steps to fit for '
             'titration curve fitting.')
    parser.add_argument(
        '--nmin', type=int,
        help='Override automatic guessing of the minimum number of waters for '
             'tittration curve fitting.')
    parser.add_argument(
        '--nmax', type=int,
        help='Override automatic guessing of maximum number of waters for '
             'titration curve fitting.')
    parser.add_argument(
        '--nfits', type=int, default=10,
        help='The number of independent fitting attempts for the neural '
             'network occupancy model. Increasing the number of fits may '
             'help improve results for noisy data.')
    parser.add_argument(
        '--pin_min', type=float, default=None,
        help='The minimum value when fitting the neural '
             'network occupancy model. Setting this may '
             'help improve models which are poorly fit at low values')
    return parser
