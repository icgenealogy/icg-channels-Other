://=============================================================================
://  eCompRB  -  full electrode compensation
://		 (procedure by Romain Brette, 2 Jan 2005)
://=============================================================================
://
:// REMARKS:
://	Compensates the membrane potential recording in a
://	single-electrode experiment in real time,
://	by substracting from the recorded potential the convolution
://	of the injected current with the electrode impulse
://	response (kernel computed elsewhere).
://
://	Input:
://		vcell		raw recorded potential	(pointer)
://		i_out		total injected current	(pointer)
://	Output:
://		vr		potential after correction
://		vr10		vr * 10
://		vr50		vr * 50
://		ue		voltage drop across the electrode
://		ue10		ue * 10
://		ue50		ue * 50
://	Parameters:
://		active		switch (1 = On, 0 = Off ; when off, vr = vcell)
://		t_start		starting time of the compensation
://	
://	The kernel is loaded in kernel[] at initialization from
://	file 'kernel_filename'.
://	
://=============================================================================
://
:// Romain Brette, Jan 2005
:// brette@ccr.jussieu.fr
://
://=============================================================================

:// Set the kernel filename here

VERBATIM
#define kernel_filename	"kernel.txt"
ENDVERBATIM

TITLE Full Electrode Compensation

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS eCompRB
	POINTER vcell,i_out
	RANGE active,active_pulses
	RANGE vr, vr10, vr50, ue, ue10, ue50
	RANGE kernel
	RANGE ksize
	RANGE imax
	RANGE t_start
://variables for Gpulses
RANGE g_e,g_i,g_eout,g_iout,del, dur, i0, amp_ge,amp_gi, pfreq, pdur,E_e, E_i

	USEION z READ ez WRITE iz VALENCE 1
}

UNITS {
	  (nA) = (nanoamp)
  (mV) = (millivolt)
(umho) = (micromho)
}

PARAMETER {
	active = 0		:// compensation switch * DEFAULT ON *
	active_pulses = 0	
	t_start = 0	(ms)	:// compensation start time
	ksize = 20		:// kernel size
	imax = 2	(nA)

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
	vcell	(mV)
	vr	(mV)
	vr10	(mV)
	vr50	(mV)
	ue	(mV)
	ue10	(mV)
	ue50	(mV)
	i_out	(nA)
	kernel[500]		:// the electrode filter
	previous_i[500]	(nA)	:// history of injected current
	iz	(nA)
	ez	(mV)

       Npulses		:// counter of pulses
  pstart	(ms)	:// start of next pulse
  pstop		(ms)	:// stop of current pulse
  g_e   (umho)          :// total excitatory conductance
  g_i   (umho)          :// total inhibitory conductance
  g_eout (umho)         :// 50000 total excitatory conductance
  g_iout (umho)      	:// 50000 total inhibitory conductance
}

INITIAL {
	FROM k = 0 TO 499 {
		previous_i[k]=0
	}
	vr = vcell
	vr10 = 10*vr
      vr50 = 50.4*vr
  g_e = 0
  g_i = 0
  Npulses = 0
  pstart = del
  pstop = pstart + pdur

	iz = 0
}

BREAKPOINT {
	SOLVE compensate
}

:// Main compensation procedure
PROCEDURE compensate() {
	if (active && (i_out > imax || i_out < -imax)) {
		i_out = 0
	}

		vr = vcell
		if (active == 1 && t > t_start) {
			:// Shift the history of current and compute the convolution
			:// (optimized fast C version)

			VERBATIM
			{
				int k,kmax;
				double *pi;
				double *kern;
				kmax = (int) ksize;
				pi = &(previous_i[kmax-1]);
				kern = &(kernel[kmax-1]);
				for(k = 1;k<kmax;k++) {
					vr-=(*pi=*(pi-1))*(*kern);
					pi--;
					kern--;
				}
				vr+=i_out*(*kern);
			}
			*previous_i = -i_out;
			ENDVERBATIM
		}
		vr10 = vr * 10
		vr50 = vr * 50.4
		ue = vcell-vr
		ue10 = ue * 10
	      ue50 = ue * 50

  if (active_pulses) {
if ((t >= pstart) && (t < dur)) { 



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

iz = g_e * (vr - E_e) 

      g_eout = 50000.0 * g_e
      ://g_iout = 50000.0 * g_i




}

:// ***************************************************************************
:// compensate() NMODL version (replaced by an optimized C version in the code)
:// ***************************************************************************
://		FROM k = 1 TO ksize-1 {
://			previous_i[ksize-k] = previous_i[ksize - k -1 ]
://		}
://		previous_i[0] = -iz	:// Beware: the injected current is inverted (ION convention)
://
://		:// Compensate the voltage recording
://		vr = vcell
://
://		FROM k = 0 TO ksize-1 {
://			vr = vr - previous_i [k] * kernel [k]
://		}
://
://		VERBATIM
://		{
://			int k;
://			for(k=0;k<(int)ksize;k++)
://				vr = vr - previous_i[k]*kernel[k];
://		}
://		ENDVERBATIM
://
://		vr10 = vr * 10
://		vr50 = vr * 50

://=============================================================================
