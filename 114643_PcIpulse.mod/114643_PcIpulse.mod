://=============================================================================
://  cIpulse  -  generator for a constant current pulse
://=============================================================================
://
://  REMARKS: 
://	- mechanism uses virtual ion z to inject current
://	- current inverted for use with DSP board
://
://  IMPLEMENTATION: Michael Rudolph, UNIC/CNRS Paris, March 2005
://	             Michael.Rudolph@iaf.cnrs-gif.fr
://
://=============================================================================

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
  POINT_PROCESS cIpulse
  RANGE active
  RANGE del, dur, amp 
  USEION z READ ez WRITE iz VALENCE 1
}

UNITS {
  (nA) = (nanoamp)
  (mV) = (millivolt)
}

PARAMETER {
  active = 0		:// 1 - process active, 0 - process not active
  
  del = 0	(ms)	:// start of pulse generation
  dur = 1e+09	(ms)	:// total duration of pulse
  amp = 0	(nA)	:// amplitude of pulse
}

ASSIGNED {
  iz	(nA)
  ez	(mV)
}

INITIAL {
  iz = 0
}

BREAKPOINT {
  if (active)  {
    if ((t >= del) && (t < dur)) { iz = -amp }
    else { iz = 0 }
  }
  else { iz = 0 }
}

PROCEDURE Update() {
  :// nothing to do
}

://=============================================================================
