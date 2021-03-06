TITLE AMPA receptors with 3-state DCO kinetic scheme
COMMENT
Author: Elena Saftenku, 2003
ENDCOMMENT
NEURON {
	POINT_PROCESS AMPA_D2
      POINTER pglu      
	NONSPECIFIC_CURRENT iampa
	RANGE  Erev, RelProb 
	RANGE gampa,gbarampa,ko,kd,kr,kb,kc,Ro
	}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
	(pS) = (picosiemens)
}
PARAMETER {	
	gbarampa	= 1200   (pS)					 
	ko	= 25.39 (/ms) 				
	kc	= 4	(/ms)				
	kd	= 5.11 (/ms) 		
	kr	= 0.065(/ms)			
      kb    =0.44 (mM) 
	Erev	= 0	(mV)
      RelProb=0.46
}
ASSIGNED {
	v (mV)		
	iampa (nA)	:current	
	gampa  (pS)	: conductance	
      pglu(mM): glutamate concentration
Ro	: receptor occupancy
celsius  (degC)
}
STATE {	
	C
	O
	D
}
INITIAL {
	C=1
	O=0
	D=0	
}

BREAKPOINT {
	SOLVE states METHOD derivimplicit
	Ro=(pglu^2+ (kb^2+2*kb*pglu)*(O+D))/(pglu+kb)^2
       gampa = gbarampa *  O * RelProb
	iampa = (1e-6) * gampa * (v - Erev)
}

DERIVATIVE states {
	LOCAL kon,kdn,Q10
Q10 = 2^((celsius-37(degC))/10(degC))
kdn = kd*pglu^2/(pglu+kb)^2
kon = ko*pglu^2/(pglu+kb)^2	
C'=Q10*(-(kdn+kon)*C+kr*D+kc*O)
D'=Q10*(-kr*D+kdn*C)
O'=Q10*(-kc*O+kon*C)	
}

