# Author: Gregory Ross

""" 
Routines to calculate free energies with GCMC as outlined in J. Am. Chem. Soc., 2015, 137 (47), pp 14930-14943
    
Can be executed from the command line as a stand-alone program
"""


import sys, os

import numpy as np
from scipy import optimize
from scipy import integrate
from scipy import special

import matplotlib
if not "DISPLAY" in os.environ or os.environ["DISPLAY"] == "" :
  matplotlib.use('Agg') 
import matplotlib.pyplot as plt

import glob
import pickle

import simulationobjects
from simulationobjects import ResultsFile

## Functions
def intersect(a, b):
    """ 
    Finds the intersection of lists. Used for arguement passing from the command line. 

    Parameters
    ----------
    a : list
    b : list
   
    Return
    ------
    list
      List of elements that are found in a and b
    """
    return list(set(a) & set(b))

## Functions specific to grand canonical integration and neural network fitting.
def pseudohuber(r,c):
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
    return (np.sqrt(1+(r/c)**2) -1)*(c**2)

def logistic(params,x):
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
    return params[2]/(1+np.exp(-params[0] - x*params[1])) 


def integrated_logistic(params,Bi,Bf):
    """ 
    The integral of the logistic function, for use in the function "insertion_pmf". The problem is that it is prone to numerical instabilities 
    when the exponents are large and positive. While one could fix the problem, I've opted for numerical integration instead in "insertion_pmf". 
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
    return -(params[2]/params[1])*( np.log(1+np.exp(Bf*params[1]+params[0])) - np.log(1+np.exp(Bi*params[1]+params[0])) )

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
    samp = np.random.uniform(0.0,1.0,size-1)
    samp = np.append([0.0,1.0],samp)
    samp = np.sort(samp)
    return np.diff(samp)

class Slp(object):
    """ 
    Class for fitting and storing a monotonically increasing single layer perceptron with one input and one output.
    
    Fitting is achieved by iteratively optimising one unit at a time, rather than the more typical back-propigation algorithm, to make use of
    SciPy's monotonic optimisation routines. 

    Attributes
    ----------
    size : integer
      the number logsistic terms, equilivalent to the number of 'hidden units' in an artificial neural network
    x : numpy array
      the explanatory variable, which is the Adams value or chemical potential in GCMC
    y : numpy array
      the response variable, which is the average number of inserted grand canonical molecules
    weights : numpy array
      the matrix of coefficients of the logistic terms
    offset : numpy array
      the intercept of the artificial neural network. Grouped together with the weights when fitting
    streams : list
      contains the output of each logistic function or 'unit'
    predicted : numpy array
      the predicted output of the artificial neural network
    error: float
      the error (either mean-squared, Huber, or absolute loss) of the artificial neural network with respect to the response variable.
    """
    def __init__(self,x=np.array(0),y=np.array(0),size=1):
        self.size = size										# The number of logistic functions a.k.a the number of units in the hidden layer.
        self.x = x												# The explanatory variable.
        self.y = y												# The response variabe; what the slp will be fitted to.
        self.weights = np.ones((size,3))						# The parameters of logistic function.
        self.offset = np.array(0.0)								# The intercept term of the slp (grouped together with the weights in most other implementations).
        self.streams = [self.x]*(size)							# The output of each of the logistic function is called a "stream".
        self.predicted = np.array(0)							# The predicted output of the slp.
        self.error = 1000.0										# The error (either mean squared, Huber, or absolute loss) of the slp w.r.t. to the response variable.
    def forward(self):
        """
        Evaluates the output from ALL units for a given set of weights and updates the prediction of the whole network
        """
        self.streams = [logistic(params,self.x) for params in self.weights]
        self.predicted = np.sum(np.vstack(self.streams),axis=0) + self.offset
    def forward_unit(self,unit):
        """
        Evaluates the output from ONE unit for a given set of weights and updates the prediction of the whole network

        Parameters
        ----------
        unit : the index of the logisitic function
        """
        self.streams[unit] = logistic(self.weights[unit],self.x)
        self.predicted = np.sum(np.vstack(self.streams),axis=0) + self.offset
    def predict(self,x):
        """
        Returns the prediction for the artificial neural network

        Parameters
        ----------
        x : the vector of Adams or chemical potentials
        """
        streams = [logistic(params,x) for params in self.weights]
        return np.sum(np.vstack(streams),axis=0) + self.offset
    def msd(self):
        """
        Calculates the mean-squared error of the model and updates the error
        """
        self.error = np.sum((self.y - self.predicted)**2)
        return self.error
    def huber(self,c=2):
        """
        Calculates the pseudo Huber loss function and updates the error of the model
        """
        self.error = np.sum(pseudohuber(self.y - self.predicted,c))
        return self.error
    def absolute(self):
        """
        Calculates the absolute loss of the model and updates the error
        """
        self.error = np.sum(np.abs(self.y - self.predicted))
        return self.error
    def randomise(self,pin_max=True):
        """
        Produces initial values of weights by random sampling

        Parameters
        ----------
        pin_max : whether to constrain the random weights so that the sum of the units cannot exceed the maximum value of the response variable
        """
        self.offset = self.y.min()
        alpha = np.random.uniform(low=0.0,high=2.0,size=self.size)
        alpha0 = np.sort(-np.random.uniform(low=self.x.min(),high=self.x.max(),size=self.size))
        if pin_max==True:
            beta = simplex_sample(self.size)*(self.y.max() - self.offset)			# The beta weights are sampled uniformally from a simplex, such that they add up to range of y values.
        else:
            beta = np.random.uniform(low=0,high=1,size=self.size)*self.y.max()
        beta.sort()												
        self.weights = np.vstack((alpha0,alpha,beta)).T								# Each row is an element in the logistic function. 
        self.forward()																# Update the logistic streams
    def randsearch(self,samples):
        """
        Performs a crude initial random search of the space of weights, and saves the best set of weights

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
    def fit_unit(self, unit, grad_tol=0.001,pin_max=None,monotonic=True, cost="msd", c=2.0):
        """
        Fits an individual unit of the artificial neural network

        Parameters
        ----------
        unit : the index of the logistic function whose weights (parameters) will be optimised
        grad_tol : the tolorance of the gradient that determines when optimisation stops
        pin_max : the value of the maximum value to model can produce. If left blank, it will be fitted
        monotonic : whether to insure the artificial neural network is monotonically increasing
        cost : the type of cost/loss function that will be used
        c : the paramater of the Huber loss function
        """
        if cost=="msd":
            def residuals(params):
                self.weights[unit] = params
                self.forward_unit(unit)
                return(self.msd())
        if cost=="abs":
            def residuals(params):
                self.weights[unit] = params
                self.forward_unit(unit)
                return(self.absolute())
        if cost=="huber":
            def residuals(params):
                self.weights[unit] = params
                self.forward_unit(unit)
                return(self.huber(c))
        if type(pin_max) == float or type(pin_max) == int:
            betas = np.ma.array(self.weights[:,2],mask=False)
            betas.mask[unit]=True
            unit_max = pin_max - self.offset - np.sum(betas)
        else:
            unit_max=None
        if monotonic==True:
            fit = optimize.fmin_l_bfgs_b(residuals,x0=self.weights[unit],bounds=[(-10.0,unit_max),(0.0,unit_max),(self.offset,unit_max)],approx_grad=True,pgtol=grad_tol,disp=0)
            self.weights[unit] = fit[0]
            self.error = fit[1]
        else:
            fit = optimize.minimize(residuals,self.weights[unit],method='BFGS',options={'gtol':grad_tol})
            self.weights[unit] = fit.x
            self.error = fit.fun
    def fit_offset(self,grad_tol=0.001,cost="msd",c=2.0):
        """
        Fits the intercept of the artificial neural network

        Parameters
        ----------
        grad_tol : the tolorance of the gradient that determines when optimisation stops
        cost : the type of cost/loss function that will be used
        c : the paramater of the Huber loss function
        """
        diff = self.y - np.sum(np.vstack(self.streams),axis=0)
        if cost=="msd":
            def residuals(offset):
                return np.sum((diff - offset)**2)
        if cost=="abs":
            def residuals(offset):
                return np.sum(np.abs(diff - offset))
        if cost=="huber":
            def residuals(offset):
                return np.sum(pseudohuber(diff - offset,c))
        fit=optimize.fmin_l_bfgs_b(residuals,x0=np.array(self.offset),bounds=[(0.0,None)],approx_grad=True,pgtol=grad_tol,disp=0)
        self.offset = fit[0]
        self.error = fit[1]
    def epoch(self,grad_tol=0.001,pin_min=None,pin_max=None,monotonic=True,cost="msd",c=2.0):
        """
        Fits each unit and intercept of the artificial neural network once

        Parameters
        ----------
        grad_tol : the tolorance of the gradient that determines when optimisation stops
        pin_min : the value of the fixed intercept. If left blank, it will be fitted.
        pin_max : the value of the maximum value to model can produce. If left blank, it will be fitted
        monotonic : whether to insure the artificial neural network is monotonically increasing
        cost : the type of cost/loss function that will be used
        c : the paramater of the huber loss function
        """
        if pin_min==None:
            self.fit_offset(grad_tol,cost,c)
        if type(pin_min) == float or type(pin_min) == int:
             self.offset = pin_min
        for unit in range(self.size):
            self.fit_unit(unit,grad_tol,pin_max,monotonic,cost,c)
    def train(self,iterations=100,grad_tol_low=-3,grad_tol_high=1,pin_min=False,pin_max=None,monotonic=True,cost="msd",c=2.0):
        """
        Iterates over a specified amount of training epochs. After each epoch, the fitting tolerance will decrease to gradually refine the fits.

        Parameters
        ----------
        iterations : the maximum number of training epochs that will be tested.
        grad_tol_low : ten to the power of this value is the initial tolorance of the gradient
        grad_tol_low : ten to the power of this value is the final tolorance of the gradient
        pin_min : the value of the fixed intercept. If left blank, it will be fitted.
        pin_max : the value of the maximum value to model can produce. If left blank, it will be fitted
        monotonic : whether to insure the artificial neural network is monotonically increasing
        cost : the type of cost/loss function that will be used
        c : the paramater of the huber loss function
        """
        tolerances = np.logspace(grad_tol_high,grad_tol_low,iterations)
        for grad_tol in tolerances:
            self.epoch(grad_tol,pin_min,pin_max,monotonic,cost,c)
        self.forward()


def inverse_slp(model,y):
    """
    For a given ANN model, this function returns the value of x that corresponds to the value of y. It's basically the inverse function of an ANN. Used to determine
    which B value produces a given average of waters.

    Parameters
    ----------
    model : Slp object
      the fitted ANN (single layer perceptron) which you want the inverse function of
    y : float
      the value which you want to find the corresponding value of x

    Return
    ------
    float
      the value of x that produced y in the input ANN
    """
    dev = (model.predicted - y)**2
    x_guess = model.x[np.where(dev == dev.min())][0]
    def errfunc(x): return (y - model.predict(x))**2
    solution = optimize.minimize(errfunc, x0 = x_guess,method="BFGS") 
    return solution.x   


def fit_ensemble(x,y,size,repeats=20,randstarts=1000,iterations=100,grad_tol_low=-3,grad_tol_high=1,pin_min=False,pin_max=None,monotonic=True,cost="msd",c=2.0,verbose=True):
    """
    Fits a collection of ANN models to GCMC titration data. Each ANN is fitted with different set of initial ANN weights.

    Parameters
    ----------
    x : numpy array
      the B/Adams values or chemical potentials used in a titration
    y : numpy array
      the average value of water at each chemical potential
    repeats : int
      the number of repeats of the fitting. This constitutes the "ensemble" of models
    others 
      As specified in Slp.

    Return
    ------
    list
      a list of Slp objects, each has been fit with different initial weights.
    """
    models = []									# A list that contains the fitted single layer perceptrons.
    error = 1E5									# Initialising the lowest error. The best model is the one with the lowest error.
    for i in range(repeats):
        model = Slp(x=x,y=y,size=size)
        model.randsearch(randstarts)
        model.train(iterations,grad_tol_low=grad_tol_low,grad_tol_high=grad_tol_high,pin_min=pin_min,pin_max=pin_max,cost=cost,c=c,monotonic=True)
        if model.error < error:
            single_model = model
        models.append( model )
        if verbose==True:
            print "Model %i fitted with error = %f" % (i+1, model.error) 
    return single_model, models


def fit_boostrap(x,y,size,boot_samps=50,repeats=10,randstarts=1000,iterations=100,grad_tol_low=-3,grad_tol_high=1,pin_min=False,pin_max=None,monotonic=True,cost="msd",c=2.0,verbose=True):
    """
    Fits a collection of ANN models to GCMC bootstrap samples of the titration data.

    Parameters
    ----------
    x : numpy array
      the B/Adams values or chemical potentials used in a titration
    y : numpy array
      the average value of water at each chemical potential
    boot_samps : int
      the number of bootstrap samples of the titration data
    repeats : int
      the number of initial weights that are attempted when fitting each ANN
    others 
      As specified in Slp.

    Return
    ------
    list
      a list of Slp objects, each has been fit to different bootstrap samples of the titration data
    """
    bootstrap_models = []
    indices = range(len(x))
    for boot in range(boot_samps):
        sample_inds = np.random.choice(indices,size=len(x))
        x_sample = x[sample_inds]
        y_sample = y[sample_inds]
        single_model, models = fit_ensemble(x_sample,y_sample,size,repeats,randstarts,iterations,grad_tol_low,grad_tol_high,pin_min,pin_max,monotonic,cost,c,False)
        bootstrap_models.append(single_model)
        if verbose == True: 
            print "...Bootstrap model %i fitted." % (boot+1)
    return bootstrap_models

def insertion_pmf(N,gcmc_model,volume,T=298.15):
    """
    Calculates the free energy to insert water from ideal gas to the GCMC volume.

    Parameters
    ----------
    N : numpy array
      vector of bound water numbers between which relative free energies will be calculated
    gcmc_model : Slp object
      a ANN model (single layer perceptron) that will be used to calculate the free energy
    T : float
      the temperature of the simulation in Kelvin

    Return
    ------
    numpy array
      relative ideal gas transfer free energies for each water specified in N
    """
    kT = T*0.0019872								# Boltmann's constant (in kcal/mol/K) multiplied by temperature (in K).
    Vstandard = 30.0
    if volume is None:
	volume = 30.0
    def integral(gcmc_model,Bi,Bf):					# For numerical integration.
        logies = [integrated_logistic(params,Bi,Bf) for params in gcmc_model.weights]
        return np.sum(logies)								# Under the convention that 0*log(0) = 0. 
    def nonintegral(N,B):
        return N*B # - special.gammaln(N+1)
    dG = np.zeros(N.size)
    correction = np.zeros(N.size)
    N_lower = N[0]
    B_lower = inverse_slp(gcmc_model,N_lower)
    for i in range(1,N.size):                         # Starting from the second element, as dG=0 at first N.
        B_upper = inverse_slp(gcmc_model,N[i])  
        #dG[i] = (nonintegral(N[i],B_upper) - nonintegral(N_lower,B_lower) + integral(gcmc_model,B_lower,B_upper))[0]							# This version uses analytical integration. Prone to instability.
        dG[i] = (nonintegral(N[i],B_upper) - nonintegral(N_lower,B_lower) - integrate.quad(gcmc_model.predict,a=B_lower,b=B_upper)[0])[0]		# This version numerical integration.
        correction[i] = (N[i]-N_lower)*np.log(volume/Vstandard)
        dG[i] = dG[i] - correction[i]
    return dG*kT


def ensemble_pmf(N,gcmc_models,T=298.15):
    energies = np.ones((len(gcmc_models),len(N)))
    for i in range(len(gcmc_models)):
        energies[i] =  insertion_pmf(N,gcmc_models[i],args.volume,T=T)
    return energies.mean(axis=0), energies.std(axis=0)

def ensemble_FreeEnergies(N,gcmc_models,hydration=None,T=None):
    """
    Calculates the free energy to insert water from ideal gas to the GCMC volume from multiple ANNs.

    Parameters
    ----------
    N : numpy array
      vector of bound water numbers between which relative free energies will be calculated.
    gcmc_models : Slp objects
      a list of ANN models (single layer perceptrons) that will be used to calculate the free energy
    hydration : float
      the hydration free energy of water in kcal/mol
    T : float
      the temperature of the simulation in Kelvin

    Return
    ------
    numpy array, numpy array
      matrix of relative ideal gas transfer free energies for each water specified in N, matrix of relative binding free energies for each water specified in N
    """

    if hydration == None: hydration = -6.2
    if T==None: T= 298.15
    dG_samples = np.zeros((len(N),len(gcmc_models)))
    for i in range(len(gcmc_models)):  
        dG_samples[:,i] = insertion_pmf(N,gcmc_models[i],args.volume)								# The free energy to transfer from the gas phase
    dG_binding_samples = dG_samples - hydration*np.tile(N,(dG_samples.shape[1],1)).T	# The free energy to transfer from solution (the binding free energy)
    return dG_samples, dG_binding_samples


# A function to find the percentiles.
def percentile_intervals(yvals,level):
    """
    Calculates the percentiles of a given input. Very much like numpy's np.percentiles. Written before I realised it existed :-/

    Parameters
    ----------
    yvals : numpy array
      vector of values that you want to find the percentiles for
    level : float
      the confidense level you want the percentiles for, eg the top 50% of the data requires level=0.5

    Return
    ------
    numpy array, numpy array, numpy array
      median, lower confidence level, upper confidence level
    """
    if level > 1:
        print "Fatal error: confidense level must be between zero and one."
        return None
    low = (1 - level)/2
    high = 1 - low
    percent_low = int(np.round(low*(yvals.shape[1]-1)))
    percent_high = int(np.round(high*(yvals.shape[1]-1)))
    mid_point = int(np.round(0.5*yvals.shape[1]))
    y_lower = np.zeros(yvals.shape[0])
    y_upper = np.zeros(yvals.shape[0])
    y_median = np.zeros(yvals.shape[0])
    for i in range(yvals.shape[0]):
        y_sorted = np.sort(yvals[i,:])
        y_lower[i] = y_sorted[percent_low]	
        y_upper[i] = y_sorted[percent_high]	
        y_median[i] = y_sorted[mid_point]
    return y_median, y_lower, y_upper


def gaussian_smooth(x,y,sigma=None):
    """
    A Gaussian smoother used in the plotting of the fitted model percentiles.

    Parameters
    ----------
    x : numpy array
      vector of the explanatory variable
    y : numpy array
      vector of the response variable    
    sigma : float
      the standard deviation of over x with which y will be smoothed with

    Return
    ------
    numpy array
      smoothed response variable
    """
    if sigma is None:
        sigma = 0.05
    smooth_y = np.zeros(len(y))
    for i in range(len(x)):
        weights = np.exp(-((x[i] - x)**2)/(2*sigma))
        weights = weights/np.sum(weights)
        smooth_y[i] = np.sum(weights*y)
    return smooth_y


def plot_FitPercentiles(B,N,gcmc_models,resolution=None,level1=None,level1_colour=None,level2=None,level2_colour=None,median_colour=None,smoothness=None):
    """
    Tool to plot a collection of fitted ANNs at two confidence limits.

    Parameters
    ----------
    B : numpy array
      vector of the explanatory variable, the B/Adams value
    N : numpy array
      vector of the response variable, the average number of waters for a given B
    gcmc_models : list of Slp objects 
      the fitted ANNs that will be plotted
    resolution : float
      the number of data points that will be plotted, i.e. resolution of plotted lines
    level1 : float
      the  confidence level (eg top 50% as a fraction) that will be plotted
    level2 : float
      the  confidence level that will be plotted
    level1_color : string
      the color of the level1 confidence region
    level2_color : string
      the color of the level2 confidence region
    median_colour : string
      colour that the median of all the fits will be plotted in
    smoothness : float
      the standard deviation over x that will be used to smooth over y

    Return
    ------
    matplotlib.pyplot object
      plot object
    """
    # Setting the defaults:
    if resolution == None: resolution = 50					# Number of points with which the lines will be drawn.
    if level1 == None: level1=0.50							# Default is also shade the region that contains 50% of the data
    if level2 == None: level2=0.90				        	# Default is NOT to shade a region at a lower percentile level.
    if level1_colour == None: level1_colour="orange"		# Default is also shade the region that contains 50% of the data
    if level2_colour == None: level2_colour="gray"			# Default is to also shade the region that contains 95% of the data
    if median_colour == None: median_colour="red"			# Colour of the median
    if smoothness == None: smoothness = 0.05				# Degree to which the lines will be smoothed (for aesthetic purposes only).
    # Generating the predictions of the different gcmc_models:
    yvals = np.zeros((resolution,len(gcmc_models)))						# Will hold all the predictions from the different models
    x = np.linspace(start=B.min(),stop=B.max(),num=resolution)		# The x values that will be plotted.
    for i in range(len(gcmc_models)):
        gcmc_models[i].x = x
        gcmc_models[i].forward()
        yvals[:,i] = gcmc_models[i].predicted
    # Generating and smoothing the confidense intervals.
    y_median, y_low, y_high = percentile_intervals(yvals,level1)
    y_median, y_lowest, y_highest = percentile_intervals(yvals,level2)
    smooth_median =gaussian_smooth(x,y_median,sigma=smoothness*3)		# The median is made smoother than the error bars.
    smooth_low = gaussian_smooth(x,y_low,sigma=smoothness)
    smooth_high = gaussian_smooth(x,y_high,sigma=smoothness)
    smooth_lowest = gaussian_smooth(x,y_lowest,sigma=smoothness)
    smooth_highest = gaussian_smooth(x,y_highest,sigma=smoothness)
    # Plotting:
    space = 0.07			# The amount of white space to the side of the x and y axis as a fraction of the ranges.
    plt.fill_between(x,smooth_lowest,smooth_highest, facecolor=level2_colour,linewidth=0,alpha=0.3,interpolate=True)
    plt.fill_between(x,smooth_low,smooth_high, facecolor="white",linewidth=0,interpolate=True)
    plt.fill_between(x,smooth_low,smooth_high, facecolor=level1_colour,linewidth=0,alpha=0.4,interpolate=True)
    plt.xlim(B.min()-space*np.ptp(B),B.max()+space*np.ptp(B))
    plt.ylim(N.min()-space*np.ptp(N),N.max()+space*np.ptp(N))
    plt.scatter(B, N,color="black",s=40)
    plt.plot(x, smooth_median,color=median_colour,linewidth=3)
    return plt

def plot_PMFpercentiles(N_range,dG_samples,level1=None,level1_colour=None,level2=None,level2_colour=None,median_colour=None):
    """
    Tool to plot a ensemble of calculated free energies

    Parameters
    ----------
    N_range : numpy array
      vector of the explanatory variable, the B/Adams value
    dG_samples : numpy array
      vector of the calculated free energies, either from repeated fits or bootstrap sampling
    gcmc_models : list of Slp objects 
      the fitted ANNs that will be plotted
    resolution : float
      the number of data points that will be plotted, i.e. resolution of plotted lines
    level1 : float
      the  confidence level (eg top 50% as a fraction) that will be plotted
    level2 : float
      the  confidence level that will be plotted
    level1_color : string
      the color of the level1 confidence region
    level2_color : string
      the color of the level2 confidence region
    median_colour : string
      colour that the median of all the fits will be plotted in
    smoothness : float
      the standard deviation over x that will be used to smooth over y

    Return
    ------
    matplotlib.pyplot object
      plot object
    """
    # Setting the defaults:
    if level1 == None: level1=0.5							# Default is also shade the region that contains 50% of the data
    if level2 == None: level2=0.90				        	# Default is NOT to shade a region at a lower percentile level.
    if level1_colour == None: level1_colour="blue"			# Default is also shade the region that contains 50% of the data
    if level2_colour == None: level2_colour="gray"			# Default is to also shade the region that contains 95% of the data
    if median_colour == None: median_colour="blue"			# Colour of the median
    g_median, g_low, g_high = percentile_intervals(dG_samples,level1)
    g_median, g_lowest, g_highest = percentile_intervals(dG_samples,level2)
    plt.plot(N_range, g_median,color=median_colour,linewidth=1)
    plt.fill_between(N_range,g_lowest,g_highest, facecolor=level2_colour,linewidth=0,alpha=0.3,interpolate=True)
    plt.fill_between(N_range,g_low,g_high, facecolor="white",linewidth=0,interpolate=True)
    plt.fill_between(N_range,g_low,g_high, facecolor=level1_colour,linewidth=0,alpha=0.3,interpolate=True)
    plt.plot(N_range, g_median,color=median_colour,linewidth=3)
    return plt


def minimum_from_free_energy(models,N_range,dG_binding_samples,print_lines=True):
    """
    Calculates the minimum value from free energy profiles. Exact for stastical thermodynamics although requires the computation of the area under the titration curve.

    Parameters
    ----------
    models : list of Slp objects
      multiple fitted ANNs
    N_range : numpy array
      vector of the number of water molecules   
    dG_binding_samples : numpy array
      collection of free energies that have been calculated with models
    print_lines : boolean
      whether to print the minimum free energy state with errors
   """
    min_inds = np.argmin(dG_binding_samples,axis=0)
    min_range = np.round(np.percentile(a=N_range[min_inds],q=[25,75]))
    min_hist,junk = np.histogram(N_range[min_inds],range=(N_range.min()-0.5,N_range.max()+0.5),bins=len(N_range))
    if args.reverse==False:
        Ntemp = N_range
    else:
        Ntemp = N_range[::-1]
    best_Bs = [inverse_slp(model,Ntemp[np.argmax(min_hist)]) for model in models]
    best_Bs = np.array(best_Bs)
    median_B = np.percentile(best_Bs,50)
    upper_B = np.percentile(best_Bs,75)
    lower_B = np.percentile(best_Bs,25)
    if print_lines:
        print "MINIMUM BINDING FREE ENERGY STATE:"
        print "Number of molecules:  Mean   Std. dev   25th Percentile    50th Percentile   75th Percentile"
        print "                      %3.1f %8.1f %14.1f %16.1f %18.1f" %(np.mean(N_range[min_inds]),N_range[min_inds].std(),int(min_range[0]),int(N_range[np.argmax(min_hist)]),int(min_range[1]))
        print "Adams value        :  Mean   Std. dev   25th Percentile    50th Percentile   75th Percentile"
        print "                      %3.1f %8.1f %14.1f %16.1f %18.1f" %(np.mean(best_Bs),best_Bs.std(),lower_B,median_B,upper_B)

# Calculates the minimum free energy state by matching the excess chemical potential to hydration free energy. Valid in the thermodynamic limit, but requires minimal data processing.
def minimum_from_thermo_limit(models,B,dG_hyd,kT=0.592,print_lines=True):
    """
    Calculates the minimum free energy state by matching the excess chemical potential to hydration free energy. Valid in the thermodynamic limit, but requires minimal data processing.

    Parameters
    ----------
    models : list of Slp objects
      multiple fitted ANNs
    B : numpy array
      vector of the Adams values
    dG_hyd : float
      hydration free energy of water in kcal/mol
    kT : float
      temperature multiplied by Boltzmann's constant, ie thermal energy, in kcal/mol
    print_lines : boolean
      whether to print the minimum free energy state with errors
   """
    best_Bs = []
    best_Ns = []
    mu_ex = []
    N_ex = []
    print dG_hyd
    for model in models:
       model.x = np.linspace(start=B.min(),stop=B.max(),num=1000)   
       model.forward()
       N_ex.append(model.predicted)
       excess = kT*(model.x - np.log(model.predicted))
       best_ind = np.argmin((excess- dG_hyd)**2)
       best_Ns.append( model.predicted[best_ind] )
       best_Bs.append( model.x[best_ind] )
       mu_ex.append(excess)
    best_Ns = np.array(best_Ns)
    best_Bs = np.array(best_Bs)
    if print_lines:
        print "THERMODYNAMIC EQUILIBRIUM STATE:"
        print "Number of molecules:  Mean   Std. dev   25th Percentile    50th Percentile   75th Percentile"
        print "                      %3.1f %8.1f %14.1f %16.1f %18.1f" %(np.mean(best_Ns),best_Ns.std(),np.percentile(best_Ns,25),np.percentile(best_Ns,50),np.percentile(best_Ns,75))
        print "Adams value        :  Mean   Std. dev   25th Percentile    50th Percentile   75th Percentile"
        print "                      %3.1f %8.1f %14.1f %16.1f %18.1f" %(np.mean(best_Bs),best_Bs.std(),np.percentile(best_Bs,25),np.percentile(best_Bs,50),np.percentile(best_Bs,75))
    return mu_ex, N_ex, best_Ns

def minimum_from_bequil(models,B,N,dG_hyd,kT=0.592,print_lines=True):
    """
    Calculates N* from Bequil.

    Parameters
    ----------
    models : list of Slp objects
      multiple fitted ANNs
    B : numpy array
      vector of the Adams values
    dG_hyd : float
      hydration free energy of water in kcal/mol
    kT : float
      temperature multiplied by Boltzmann's constant, ie thermal energy, in kcal/mol
    print_lines : boolean
      whether to print the minimum free energy state with errors
   """
    if args.volume is None:
	print 'B EQUILIBRIUM CONDITION CANNOT BE CALCULATED WITHOUT THE GCMC VOLUME'
	print 'please use -v flag to calculate the optimum water occupancy'
	return None, None, None
    else:
	Bequil = dG_hyd*(1.0/0.592) + np.log(args.volume/30.0)
	closestBfit = (np.abs(single_model.x - Bequil)).argmin() 
	Nstarfit = []
	Nstarcalc = []
	
	if Bequil < B.min() or Bequil > B.max():
	    	print 'Equilibrium B value is not within range sampled in simulation, so an estimation of N* cannot be made.'
		return Bequil, Nstarcalc, Nstarfit
	else:
	    for model in models:
	        model.x = np.linspace(start=B.min(),stop=B.max(),num=len(model.predicted))
	        Nstarfit.append(model.predicted[closestBfit])
	    Nstarfit = np.array(Nstarfit)
	    if print_lines:
	        print "B EQUILIBRIUM CONDITION FROM FITTED MODELS:"
	        print "Bequil:",np.round(Bequil,1)
	        print "Number of molecules:  Mean   Std. dev   25th Percentile    50th Percentile   75th Percentile"
	        print "                      %3.1f %8.1f %14.1f %16.1f %18.1f" %(np.mean(Nstarfit),Nstarfit.std(),np.percentile(Nstarfit,25),np.percentile(Nstarfit,50),np.percentile(Nstarfit,75))
	    
	    Bclose = []
	    nearBvalues = np.abs(B - Bequil)
	    for index,each in enumerate(nearBvalues):
	        if each <=nearBvalues.min():                    # this works if Bequil is exactly between two B values, it will find an average of both
	            Bclose.append(B[index])
	            Nstarcalc.append(N[index])
	    Bclose = np.array(Bclose)
	    Nstarcalc = np.array(Nstarcalc)
	    print "B EQUILIBRIUM CONDITION FROM SIMULATED B VALUE:"
	    print "Bequil:",np.round(Bequil,1)
	    print "Similar simulated B values:", np.unique(Bclose) 
	    print "Average N* at simululated at B values:",np.unique(Bclose),"is", np.round(Nstarcalc.mean(),1)
	
	return Bequil, Nstarcalc, Nstarfit

def get_arg_parser():
  import argparse

  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Program to analyse and plot free energies from GCMC simulations")
  parser.add_argument('-d','--directories',nargs="+",help="the directories containing the GCMC simulation data, default=None",default=None)
  parser.add_argument('-f','--file',help="the name of the file to analyse. Default is results.",default="results")
  parser.add_argument('-p','--plot',nargs="+",choices=["titration","fit","percentiles","pmf","all"],help="whether to plot the GCMC simulation data or analysis, default=None",default=None)
  parser.add_argument('-c','--calc',nargs="+",choices=["fit","pmf","minimum","excess","all"], help="fit an artificial neural network to the data, the potential of mean force to bind a specified number of waters, locate the minimum of the pmf",default=None)
  parser.add_argument('-s','--skip',help="the number of initial snapshots that will be discarded when fitting",type=int,default=0)
  parser.add_argument('-b','--bootstraps',help="the number of bootstrap samples to perform.",type=int,default=None)
  parser.add_argument('-i','--input',help="the pickled neural networks that will be read, bypassing the fitting procedure, default=ANNs.pickle",default=None)
  parser.add_argument('-o','--out',help="the pickled file name where the fitted neural networks will be output, default=ANNs.pickle",default=None)
  parser.add_argument('--steps', help="the number of steps observed in the titration data that will be fitted to",type=int,default=1)
  parser.add_argument('--range',nargs="+",help="range of water molecules to calculate free energy",default=None)
  parser.add_argument('--reverse',action='store_true',help="reverse the direction of integration, starting from the highest number of waters to the lowest number, default=False",default=False)
  parser.add_argument('--fit_options',help="additional options to be passed to the artificial neural network",default=None)
  parser.add_argument('-v','--volume', type=float,help="volume of the GCMC insertion region",default=None)
  return parser
  
if __name__ == '__main__' :

  args = get_arg_parser().parse_args()

  dG_hyd = -6.2
  print "CALCULATING VOLUME CORRECTION:"
  if args.volume is None:
      print "No GCMC volume given. No volume correction will be applied.\n A correction of kBT ln(Vsystem/Vstandard) is required for GCMC calculations"
  else:
      print "GCMC volume:", args.volume
      print "Standard volume: 30.0"
      print "Volume correction of kBT ln (Vsystem/Vstandard) will be applied" 

  if args.directories is None:
      print "\nError. Please list the directories containing the GCMC simulation data. Exiting program.\n"
      sys.exit()

  if args.plot==None and args.calc == None:
      print "\nNo plots or calculations specified. Exiting program.\n"
      sys.exit()

  if args.plot is not None:
      if len(intersect(args.plot,["fit","percentiles", "excess","pmf","all"])) >= 1 and args.calc == None and args.input == None:
          print "\nError. You are asking to plot things without specifying the corresponding calculation or have not input an ANN. For instance, try typing '--plot titration' or, add '--calc all' when calling this script. Exiting program.\n"
          sys.exit()

  ANNs=None
  if args.input is not None:
    try:
        ANNs = pickle.load( open( args.input, "rb" ) )
        single_model = ANNs["single"]
        models = ANNs["collection"]
    except IOError:
        print "\n Warning. No artificial neural network called '%s' detected. Exiting program" % args.input
        sys.exit()

  # Read in the GCMC results data from multiple GCMC output folders and calculating the mean number of "on" waters, after skipping a specified number of frames.
  directories = args.directories
  N = []
  B = []
  print "\nREADING GCMC DATA:"
  for dirs in directories:
    folders =  glob.glob(dirs)
    if len(folders)==0:
        print "\nError. No folder(s) matching '%s'. Exiting program.\n" % args.directories
        sys.exit()
    for folder in folders:
        results = ResultsFile()
        fpath = os.path.join ( folder, args.file )
        resultsfiles = [ fpath ] if os.path.exists ( fpath ) else []
        resultsfiles += glob.glob(fpath+"-*")
        
        if len(resultsfiles) > 1: # It is assumed the results are from ProtoMS 2.
            results.read([folder,args.file])
        elif len(resultsfiles)==1:# It is assumed the results are from ProtoMS 3.
            results.read(folder+ "/"+args.file)
        else:
           print "\nError. No results file matching %s. Exiting program\n" % folder+"/"+args.file 
           sys.exit()
        solventson = np.array([snap.solventson for snap in results.snapshots])			# Counting how many molecules have been inserted.
        mean_on = np.mean(solventson[args.skip:solventson.size])
        N.append(mean_on)
        adams = np.array([snap.bvalue for snap in results.snapshots])					# Recording the B-value for each window (should be constant, but averaging just in case).
        mean_adams = np.mean(adams[(args.skip-1):adams.size])
        B.append(mean_adams)

  # Checking to make sure the data makes sense.
  N = np.array(N)
  B = np.array(B)
  if B.size != N.size:
      print "\nFatal error: Length of array of Adams values must equal the length of the array of the number of inserted waters.\n"
      sys.exit()
  else:
      print "......GCMC data has been read.\n"   

  
  # Specifying the default options for artificial neural network, and reading in values from the command line.
  fit_dict = {"monotonic":True, "repeats":10, "randstarts":1000, "iterations":100, "grad_tol_low":-3, "grad_tol_high":1, "pin_min":None, "pin_max":None, "cost":"msd", "c":2.0, "verbose":False}
  if args.fit_options is not None:
      options = args.fit_options.split()
      for i in range(len(options)):
          if options[i] in fit_dict:
              if type(fit_dict[options[i]]) == str:
                  fit_dict[options[i]]= str(options[i+1])
              if type(fit_dict[options[i]]) == float or fit_dict[options[i]] == None:
                  fit_dict[options[i]]= float(options[i+1])
              if type(fit_dict[options[i]]) == int:
                  fit_dict[options[i]]= int(options[i+1])
              if type(fit_dict[options[i]]) == bool:
                  fit_dict[options[i]]= bool(options[i+1])

  # Fitting the artificial neural network, which is crucial for all other analysis.
  if args.calc is not None and ANNs is None:
      print "FITTING TO TITRATION DATA:"
      if args.bootstraps is None:
          single_model, models = fit_ensemble(x=B,y=N,size=args.steps,verbose=fit_dict["verbose"],pin_min=fit_dict["pin_min"],pin_max=fit_dict["pin_max"],cost=fit_dict["cost"],c=fit_dict["c"],randstarts=fit_dict["randstarts"],repeats=fit_dict["repeats"],iterations=fit_dict["iterations"])
          print "......Neural network fitting complete."      
      else:
          samples = args.bootstraps
          models = fit_boostrap(x=B,y=N,size=args.steps,boot_samps=samples,pin_min=fit_dict["pin_min"],cost=fit_dict["cost"],c=fit_dict["c"],randstarts=fit_dict["randstarts"],repeats=fit_dict["repeats"],iterations=fit_dict["iterations"])              
          single_model, rubbish = fit_ensemble(x=B,y=N,size=args.steps,verbose=fit_dict["verbose"],pin_min=fit_dict["pin_min"],cost=fit_dict["cost"],c=fit_dict["c"],randstarts=fit_dict["randstarts"],repeats=fit_dict["repeats"],iterations=fit_dict["iterations"])
          print "......Neural network bootstrap fitting complete." 
      if args.out is not None:    
          ANNs = {"single":single_model,"collection":models}
          pickle.dump( ANNs, open( args.out, "wb" ) )

  # Defining the number of water molecules that will be used for integration and plotting.
  if args.range == None:										# The default is to integrate between the integers present in N.
      N_range = np.arange(np.round(N.min()),np.round(N.max())+1)
  else:
      b1 = float(args.range[0])
      b2 = float(args.range[1])
      if b2-b1 > 0 : 
          N_range = np.arange(b1,b2+1)
      else: 
          N_range = np.arange(b1,b2-1,-1)

  # Calculating free energies using integration and the fitted neural network.
  if args.calc is not None:
    if len(intersect(args.calc,["pmf","all"])) > 0:
      if args.reverse==False:
          dG_single = insertion_pmf(N_range,single_model,args.volume)
          dG_samples, dG_binding_samples = ensemble_FreeEnergies(N_range,models)
      else: 
          N_range = N_range[::-1]
          dG_single = insertion_pmf(N_range,single_model,args.volume)
          dG_samples, dG_binding_samples = ensemble_FreeEnergies(N_range,models)

    if len(intersect(args.calc,["pmf","all"])) > 0:
      results = np.vstack((N_range,dG_samples.mean(axis=1), dG_samples.std(axis=1),np.percentile(dG_samples,25,axis=1),np.percentile(dG_samples,50,axis=1),np.percentile(dG_samples,75,axis=1) )).T
#      if args.bootstraps is None:
#          error_type = "'Std. dev. of fit'"
#      else:
#          error_type = "'Std. dev. of bootstrap'"
      print "\nFREE ENERGIES:"
      if args.bootstraps is None and args.input is None:
          print "  Quoted errors are from multiple repeats of the fitting."
      elif args.bootstraps is not None and args.input is None:
          print "  Quoted errors are from %i bootstrap samples." % args.bootstraps
      elif  args.input is not None:
          print "  Quoted errors are from the input model(s)." 
      print "          |----------------------IDEAL GAS TRANSFER FREE ENERGIES--------------------|   |-BINDING FREE ENERGIES-|"
      print "'# Waters' 'Mean'  'Std. dev.'  '25th Percentile'       'Median'      '75th Percentile'    'Mean'        'Median'"
      for row in results:
          print " %5.2f %9.2f %9.2f %15.2f %19.2f %18.2f %15.2f %14.2f" % (row[0], row[1], row[2],row[3],row[4],row[5],row[1] - dG_hyd*(row[0]-N_range[0]),row[4] - dG_hyd*(row[0]-N_range[0]) )
      print "Binding free energies are transfer free energies minus the hydration free energy of water (%2.1f kcal/mol) multiplied by the number of waters.\n" % dG_hyd

    # Calculating the equilibrium number of bound waters.
    if len(intersect(args.calc,["minimum","pmf","all"])) > 0:
      if len(intersect(args.calc,["pmf","all"])) > 0:
        print "\n"
        minimum_from_free_energy(models,N_range,dG_binding_samples) 
        print "\n"
        Bequil, Nstarcalc, Nstarfit = minimum_from_bequil(models,B,N,dG_hyd,kT=0.592)
      else:
        print "\n"
        Bequil, Nstarcalc, Nstarfit = minimum_from_bequil(models,B,N,dG_hyd,kT=0.592)
 

  # Plotting the requested results.
  FigNum = 0
  if args.plot is not None:
    if len(intersect(args.plot,["titration","fit","all"])) > 0:
      FigNum += 1
      plt.figure("GCMC Titration")
      currfig = plt
      currfig.scatter(B, N,color="black")
      if len(intersect(args.plot,["fit", "all"])) > 0:
        single_model.x = np.linspace(start=B.min(),stop=B.max(),num=100)
        single_model.forward()
        currfig.plot(single_model.x,single_model.predicted,color="red",linewidth=3)
      currfig.xlabel("Adams parameter (B)",fontsize=15)
      currfig.ylabel("Average number of waters",fontsize=15)
      if args.calc is not None:
        currfig.suptitle("GCMC titration data and fitted model",fontweight="bold")
      else:
        currfig.suptitle("GCMC titration data",fontweight="bold")   
      currfig.savefig("Titration.png")
      currfig.show(block=False)

    if len(intersect(args.plot,["percentiles", "all"])) > 0:
      FigNum += 1
      plt.figure("Fitted ANN Percentiles")
      currfig = plot_FitPercentiles(B,N,models,level1=0.50,level2=0.90)
      currfig.xlabel("Adams parameter (B)",fontsize=15)
      currfig.ylabel("Average number of waters",fontsize=15)
      currfig.suptitle("Median and percentiles of fitted models",fontweight="bold")
      currfig.savefig("Fit_percentiles.png")
      currfig.show(block=False)

    if len(intersect(args.plot,["pmf", "all"])) > 0 and len(intersect(args.calc,["pmf", "all"])) > 0 :
      FigNum += 1 
      plt.figure("Binding Free Energy")
      currfig = plot_PMFpercentiles(N_range,dG_binding_samples)
      currfig.xlabel("Number of inserted waters",fontsize=15)
      currfig.ylabel("Binding free energy (kcal/mol)",fontsize=15)
      currfig.suptitle("Binding free energy profile",fontweight="bold")
      currfig.savefig("Binding_Free_Energy.png")
      currfig.show(block=False)

    # if len(intersect(args.plot,["excess", "all"])) > 0 and len(intersect(args.calc,["minimum","excess","all",])) > 0:
    #   mu_ex, N_ex, best_Ns = minimum_from_bequil(models,B,dG_hyd,kT=0.592,print_lines=False)
    #   FigNum += 1 
    #   plt.figure("Excess Chemical Potential")
    #   for u in range(len(models)):
    #     plt.plot(N_ex[u],mu_ex[u],alpha=0.1,color="blue",linewidth=3)
    #   plt.axhline(y=dG_hyd,color="grey",linewidth=2)
    #   plt.axvline(x=np.percentile(best_Ns,50),color="grey",linewidth=2)
    #   plt.xlabel("Number of inserted waters",fontsize=15)
    #   plt.savefig("Excess_Chem_Potential.png")
    #   plt.ylabel("Excess chemical potential (kcal/mol)",fontsize=15)
    #   plt.show(block=False)
  
  
  
  if args.plot is not None:      
      print "\nType enter to quit\n>"
      raw_input()

