://=============================================================================
://  pGnoise  -  fluctuating conductance model for synaptic bombardment
://=============================================================================
://
://  THEORY:
://    Synaptic bombardment is represented by a stochastic model containing
://    two fluctuating conductances g_e(t) and g_i(t) descibed by:
://
://      Isyn = g_e(t) * [V - E_e] + g_i(t) * [V - E_i]
://      d g_e / dt = -(g_e - g_e0) / tau_e + sqrt(D_e) * Ft
://      d g_i / dt = -(g_i - g_i0) / tau_i + sqrt(D_i) * Ft
://
://    where E_e, E_i are the reversal potentials, g_e0, g_i0 are the average
://    conductances, tau_e, tau_i are time constants, D_e, D_i are noise
://    diffusion coefficients and Ft is a gaussian white noise of unit
://    standard deviation.
://
://    g_e and g_i are described by an Ornstein-Uhlenbeck (OU) stochastic
://    process where tau_e and tau_i represent the "correlation" (if tau_e
://    and tau_i are zero, g_e and g_i are white noise). The estimation
://    of OU parameters can be made from the power spectrum:
://
://      S(w) =  2 * D * tau^2 / (1 + w^2 * tau^2)
://
://    and the diffusion coeffient D is estimated from the variance:
://
://      D = 2 * sigma^2 / tau
://
://  NUMERICAL RESOLUTION:
://    The numerical scheme for integration of OU processes takes advantage
://    of the fact that these processes are gaussian, which led to an exact
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
://     - conductance cannot be negative!
://
://  PARAMETERS:
://    The mechanism takes the following parameters:
://
://    E_e = 0  (mV)            : reversal potential of excitatory conductance
://    E_i = -75 (mV)           : reversal potential of inhibitory conductance
://
://    g_e0 = 0.0121 (umho)     : average excitatory conductance
://    g_i0 = 0.0573 (umho)     : average inhibitory conductance
://
://    std_e = 0.0030 (umho)    : standard dev of excitatory conductance
://    std_i = 0.0066 (umho)    : standard dev of inhibitory conductance
://
://    tau_e = 2.728 (ms)       : time constant of excitatory conductance
://    tau_i = 10.49 (ms)       : time constant of inhibitory conductance
://
://  IMPLEMENTATION: A. Destexhe, Laval University, 1999
://
://  IMPLEMENTATION (for use in experiments): Michael Rudolph, 
://					      UNIC/CNRS Paris, March 2005
://
://=============================================================================

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
  POINT_PROCESS pGnoise
  POINTER vcell
  RANGE active
  RANGE g_e, g_i, E_e, E_i, g_e0, g_i0, g_e1, g_i1
  RANGE std_e, std_i, tau_e, tau_i, D_e, D_i
  RANGE del, dur, seed
  RANGE g_eout, g_iout
  USEION z READ ez WRITE iz VALENCE 1
}

UNITS {
  (nA) = (nanoamp)
  (mV) = (millivolt)
  (umho) = (micromho)
}

PARAMETER {
  active = 0		:// 1 - process active, 0 - process not active
  
  dt (ms)

  del = 0	  (ms)
  dur = 5000 	  (ms)
  seed = 0
  
  E_e= 0 	  (mV)  :// reversal potential of excitatory conductance
  E_i= -75 	  (mV)  :// reversal potential of inhibitory conductance

  g_e0= 0.0121    (umho):// average excitatory conductance
  g_i0= 0.0573    (umho):// average inhibitory conductance

  std_e= 0.0030   (umho):// standard dev of excitatory conductance
  std_i= 0.0066   (umho):// standard dev of inhibitory conductance

  tau_e= 2.728    (ms)  :// time constant of excitatory conductance
  tau_i= 10.49    (ms)  :// time constant of inhibitory conductance
}

ASSIGNED {
  vcell	(mV)		:// membrane voltage
  iz    (nA)        	:// fluctuating current
  ez	(mV)
  g_e   (umho)          :// total excitatory conductance
  g_i   (umho)          :// total inhibitory conductance
  g_e1  (umho)          :// fluctuating excitatory conductance
  g_i1  (umho)          :// fluctuating inhibitory conductance
  D_e   (umho umho /ms) :// excitatory diffusion coefficient
  D_i   (umho umho /ms) :// inhibitory diffusion coefficient
  exp_e
  exp_i
  amp_e (umho)
  amp_i (umho)

  g_eout (umho)         :// 50000 total excitatory conductance
  g_iout (umho)      	:// 50000 total inhibitory conductance
}

INITIAL {
  set_seed(seed)
  
  exp_e = 0
  exp_i = 0
  amp_e = 0
  exp_i = 0

  g_e1 = 0
  g_i1 = 0
  
  if(tau_e != 0) {
    D_e = 2 * std_e * std_e / tau_e
    exp_e = exp(-dt/tau_e)
    amp_e = std_e * sqrt( (1-exp(-2*dt/tau_e)) )
  }
  
  if(tau_i != 0) {
    D_i = 2 * std_i * std_i / tau_i
    exp_i = exp(-dt/tau_i)
    amp_i = std_i * sqrt( (1-exp(-2*dt/tau_i)) )
  }
}

BREAKPOINT {
  if (active)  {
    if ((t >= del) && (t < dur)) {
      SOLVE oup
      
      if(tau_e==0) { 
        g_e1 = std_e * normrand(0,1) 
      }
      
      if(tau_i==0) { 
        g_i1 = std_i * normrand(0,1) 
      }
      
      g_e = g_e0 + g_e1
      if(g_e < 0) { g_e = 0 }
      
      g_i = g_i0 + g_i1
      if(g_i < 0) { g_i = 0 }
      
      iz = g_e * (vcell - E_e) + g_i * (vcell - E_i)

      g_eout = 50000.0 * g_e
      g_iout = 50000.0 * g_i
      
    } else { 
      iz = 0 
      g_eout = 0
      g_iout = 0
    }
  }
  else { 
    iz = 0 
    g_eout = 0
    g_iout = 0
  }
}

PROCEDURE oup() {       :// use Scop function normrand(mean, std_dev)
  if(tau_e!=0) { 
    g_e1 =  exp_e * g_e1 + amp_e * normrand(0,1) 
  }
  
  if(tau_i!=0) { 
    g_i1 =  exp_i * g_i1 + amp_i * normrand(0,1) 
  }
}

PROCEDURE Update() {
  set_seed(seed)
  
  if(tau_e != 0) {
    D_e = 2 * std_e * std_e / tau_e
    exp_e = exp(-dt/tau_e)
    amp_e = std_e * sqrt( (1-exp(-2*dt/tau_e)) )
  }
  
  if(tau_i != 0) {
    D_i = 2 * std_i * std_i / tau_i
    exp_i = exp(-dt/tau_i)
    amp_i = std_i * sqrt( (1-exp(-2*dt/tau_i)) )
  }
}

://=============================================================================
