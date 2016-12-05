://=============================================================================
://  pIpulse  -  generator for periodic conductance pulses with chosen Erev
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
  POINT_PROCESS pGpulse
  RANGE active
  RANGE g_e,g_i,g_eout,g_iout,del, dur, i0, amp_ge,amp_gi, pfreq, pdur,E_e, E_i
  USEION z READ ez WRITE iz VALENCE 1
  POINTER vcell
}

UNITS {
  (nA) = (nanoamp)
  (mV) = (millivolt)
(umho) = (micromho)
}

PARAMETER {
  active = 0		:// 1 - process active, 0 - process not active
  
  del = 0	(ms)	:// start of pulse generation
  dur = 10000	(ms)	:// total duration of pulse
  ://i0 = 0	(nA)	:// offset current
  E_e= 0 	  (mV)  :// reversal potential of excitatory conductance
  E_i= -75 	  (mV)  :// reversal potential of inhibitory conductance

  amp_ge = 0	(umho)	:// amplitude of excitation pulse
  amp_gi = 0	(umho)	:// amplitude of inhibition pulse
  pfreq = 10	(Hz)	:// frequency of pulses
  pdur = 500	(ms)	:// length of each pulse
}

ASSIGNED {
  vcell	(mV)		:// membrane voltage
  iz		(nA)
  ez		(mV)
  Npulses		:// counter of pulses
  pstart	(ms)	:// start of next pulse
  pstop		(ms)	:// stop of current pulse
  g_e   (umho)          :// total excitatory conductance
  g_i   (umho)          :// total inhibitory conductance
  g_eout (umho)         :// 50000 total excitatory conductance
  g_iout (umho)      	:// 50000 total inhibitory conductance
}

INITIAL {
  iz = 0
  g_e = 0
  g_i = 0
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

      g_e = amp_ge
      ://g_i = 0
      Npulses = Npulses + 1
      pstart = del + Npulses*1000/pfreq
      if (1000/pfreq <= pdur) { pstop = del + dur }

          } 
    else if (t >= pstop) { 
      g_e = 0
      ://g_i = amp_gi

      pstop = pstart + pdur
    }
  }
  else { g_e = 0 
         ://g_i = 0
}

iz = g_e * (vcell - E_e) 

      g_eout = 50000.0 * g_e
      ://g_iout = 50000.0 * g_i



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
