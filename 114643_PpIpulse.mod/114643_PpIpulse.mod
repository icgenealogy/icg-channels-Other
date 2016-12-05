://=============================================================================
://  pIpulse  -  generator for periodic current pulse
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
  POINT_PROCESS pIpulse
  RANGE active
  RANGE del, dur, i0, amp, pfreq, pdur
  USEION z READ ez WRITE iz VALENCE 1
}

UNITS {
  (nA) = (nanoamp)
  (mV) = (millivolt)
}

PARAMETER {
  active = 0		:// 1 - process active, 0 - process not active
  
  del = 0	(ms)	:// start of pulse generation
  dur = 10000	(ms)	:// total duration of pulse
  i0 = 0	(nA)	:// offset current
  amp = 0	(nA)	:// amplitude of pulse
  pfreq = 1	(Hz)	:// frequency of pulses
  pdur = 500	(ms)	:// length of each pulse
}

ASSIGNED {
  iz		(nA)
  ez		(mV)
  Npulses		:// counter of pulses
  pstart	(ms)	:// start of next pulse
  pstop		(ms)	:// stop of current pulse
}

INITIAL {
  iz = 0
  Npulses = 0
  pstart = del
  pstop = pstart + pdur

	VERBATIM
 	printf("%f \n", pstart);
	ENDVERBATIM
}

BREAKPOINT {
  if (active) {
    if ((t >= pstart) && (t < dur)) { 

	VERBATIM
 	printf("%f \n", active);
	ENDVERBATIM

      iz = -(amp + i0)
      Npulses = Npulses + 1
      pstart = del + Npulses*1000/pfreq
      if (1000/pfreq <= pdur) { pstop = del + dur }
    } 
    else if (t >= pstop) { 
      iz = -i0 
      pstop = pstart + pdur
    }
  }
  else { iz = 0 }
}

PROCEDURE Update() {
  VERBATIM
    Npulses = (int) ((t-del)*pfreq/1000);
  ENDVERBATIM
  
  pstart = del + Npulses*1000/pfreq
  
  if (1000/pfreq <= pdur) { pstop = del + dur }
  else { pstop = pstart + pdur }
}

://=============================================================================
