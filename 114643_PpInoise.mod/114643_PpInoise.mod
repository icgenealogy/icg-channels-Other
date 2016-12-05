://=============================================================================
://  pInoise - fluctuating current model for synaptic
://            background activity (static i + fluctuating i)
://=============================================================================
://
://  THEORY:
://    Synaptic background activity is represented by a fluctuating
://    current I_syn(t) described by an one-dimensional Ornstein-Uhlenbeck
://    stochastic process:
://
://      d I_syn / dt = -(I_syn - I_0) / tau_I + sqrt(D_I) * Ft
://
://    where I_0 denotes the average current, tau_I the noise
://    time constant, D_I the noise diffusion coefficient and Ft
://    Gaussian white noise of unit standard deviation.
://
://    If tau_I is zero, I_syn describes white noise. D_I and the
://    standard deviation sigma_I of the noise are related by:
://
://      D_I = 2 * sigma_I^2 / tau_I
://
://  NUMERICAL RESOLUTION:
://    The numerical scheme for integration of the OU process takes advantage
://    of the fact that the process is Gaussian, which led to an exact
://    update rule independent of the time step dt (see Gillespie DT,
://    Am J Phys 64: 225, 1996):
://
://      x(t+dt) = x(t) * exp(-dt/tau) + A * N(0,1)
://
://    where A = sqrt( D*tau/2 * (1-exp(-2*dt/tau)) ) and N(0,1) is a normal
://    random number (avg=0, sigma=1)
://
://  REMARKS:
://	- mechanism uses virtual ion z to inject current
://	- current inverted for use with DSP board
://
://  PARAMETERS:
://    The mechanism takes the following standard parameters:
://
://    I_0 = -0.257 (nA)        : average current (-65mV for E_pas=-80mV)
://    std_I = 0.25 (nA)        : standard deviation (sigma_V around 4mV)
://    tau_I = 2.0 (ms)         : time constant of noise
://
://  IMPLEMENTATION: Michael Rudolph, UNIC/CNRS Paris, March 2005
://
://=============================================================================

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
  POINT_PROCESS pInoise
  RANGE active
  RANGE I_0, std_I, tau_I, D_I
  RANGE del, dur, seed
  USEION z READ ez WRITE iz VALENCE 1
}

UNITS {
  (nA) = (nanoamp)
}

PARAMETER {
  active = 0		:// 1 - process active, 0 - process not active
  
  dt              (ms)

  del = 0	  (ms)
  dur = 5000 	  (ms)
  seed = 0
  
  I_0    = 0.1    (nA)  :// mean current
  std_I  = 0.1    (nA)  :// standard deviation of noisy current
  tau_I  = 2.0    (ms)  :// time constant of noise
}

ASSIGNED {
  iz    (nA)        	:// fluctuating current
  ez	(mV)
  iz1   (nA)            :// fluctuating current
  D_I   (nA nA /ms)     :// diffusion coefficient
  exp_I
  amp_I (nA)
}

INITIAL {
  set_seed(seed)
  
  iz1 = 0
  
  if(tau_I != 0) {
    D_I = 2 * std_I * std_I / tau_I
    exp_I = exp(-dt/tau_I)
    amp_I = std_I * sqrt( (1-exp(-2*dt/tau_I)) )
  }
}

BREAKPOINT {
  if (active)  {
    if ((t >= del) && (t < dur)) {
      SOLVE oup
      iz = - I_0 - iz1
    } else { iz = 0 }
  }
  else { iz = 0 }   
}

PROCEDURE oup() {       :// use Scop function normrand(mean, std_dev)
  if (tau_I!=0) { iz1 = exp_I * iz1 + amp_I * normrand(0,1) }
  if (tau_I==0) { iz1 = std_I * normrand(0,1) }
}

PROCEDURE Update() {
  set_seed(seed)
  
  if(tau_I != 0) {
    D_I = 2 * std_I * std_I / tau_I
    exp_I = exp(-dt/tau_I)
    amp_I = std_I * sqrt( (1-exp(-2*dt/tau_I)) )
  }
}

://=============================================================================
