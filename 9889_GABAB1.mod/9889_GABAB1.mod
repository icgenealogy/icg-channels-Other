NEURON {  POINT_PROCESS GABAB1 }
:   Simple synaptic mechanism derived for first order kinetics of
:   binding of transmitter to postsynaptic receptors.
:
:     GABA SYNAPSE (GABA-B receptors)
:
:     Parameters estimated from patch clamp recordings of GABAB PSP's in
:     rat hippocampal slices (Otis et al, J. Physiol. 463: 391-407, 1993)
:     and sharp electrode recordings in hippocampal pyramidal cells 
:     (Solis & Nicoll, J. Neurosci. 12: 3466-3472, 1992).  
:
:     To account for the time course of the GABAB IPSP with a first order
:     kinetics, it is assumed that the transmitter stays for 150 ms in the 
:     synaptic cleft.  This is very unrealistic, but it provides reasonable
:     summation of PSP's.
:
:   A. Destexhe , The Salk Institute, May 1993.
PARAMETER {
  Cmax    = 1     (mM)            : max transmitter concentration
  Cdur    = 150   (ms)            : transmitter duration (rising phase)
  Alpha   = 0.01  (/ms mM)        : forward (binding) rate
  Beta    = 0.005 (/ms)           : backward (unbinding) rate
  Erev    = -95   (mV)            : reversal potential (potassium)
  Prethresh = 0                   : voltage level nec for release
  Deadtime = 1    (ms)            : mimimum time between release events
  gmax            (umho)          : maximum conductance
}
INCLUDE "synq.inc"
:________________________________________________________________

