TITLE transmitter release

COMMENT
-----------------------------------------------------------------------------

 Simple (minimal?) model of transmitter release

 - single compartment, need calcium influx and efflux

 - Ca++ binds to a "fusion factor" protein F leading to an activated form FA.
   Assuming a cooperativity factor of 4 (see Augustine & charlton, 
   J Physiol. 381: 619-640, 1986), one obtains:

	F + 4 Cai <-> FA	(kb,ku)

 - FA binds to presynaptic vesicles and activates them according to:

	FA + V <-> VA		(k1,k2)

   VA represents the "activated vesicle" which is able to bind to the
   membrane and release transmitter.  Presynaptic vesicles (V) are 
   considered in excess.

 - VA releases nt transmitter molecules in the synaptic cleft

	VA  ->  nt T		(k3)

   This reaction is the slowest and a constant number of transmitter per 
   vesicule is considered (nt).  

 - Finally, T is hydrolyzed according to a first-order reaction

	T  ->  ...		(kh)


   References:

   Destexhe, A., Mainen, Z.F. and Sejnowski, T.J. Synthesis of models for
   excitable membranes, synaptic transmission and neuromodulation using a 
   common kinetic formalism, Journal of Computational Neuroscience 1: 
   195-230, 1994.

   Destexhe, A., Mainen, Z.F. and Sejnowski, T.J.  Kinetic models of 
   synaptic transmission.  In: Methods in Neuronal Modeling (2nd edition; 
   edited by Koch, C. and Segev, I.), MIT press, Cambridge, 1998, pp 1-25.

  (electronic copy available at http://cns.iaf.cnrs-gif.fr)

   For a more realistic model, see Yamada, WM & Zucker, RS. Time course
   of transmitter release calculated from simulations of a calcium
   diffusion model. Biophys. J. 61: 671-5682, 1992.


  Written by A. Destexhe, Salk Institute, December 1993; modified 1996

-----------------------------------------------------------------------------
ENDCOMMENT


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX rel
	USEION ca READ cai WRITE cai
	RANGE T,FA,CA,Fmax,Ves,b,u,k1,k2,k3,nt,kh
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mM) = (milli/liter)
}

PARAMETER {

	Ves = 0.1 	(mM)		: conc of vesicles
	Fmax = 0.001	(mM)		: conc of fusion factor F
	b = 1e16 	(/mM4-ms)	: ca binding to F
	u = 0.1  	(/ms)		: ca unbinding 
	k1 = 1000   	(/mM-ms)	: F binding to vesicle
	k2 = 0.1	(/ms)		: F unbinding to vesicle
	k3 = 4   	(/ms)		: exocytosis of T
	nt = 10000			: nb of molec of T per vesicle
	kh = 10  	(/ms)		: cst for hydolysis of T
}

ASSIGNED {
}

STATE {
	FA	(mM)
	VA	(mM)
	T	(mM)
	cai	(mM) 
}

INITIAL {
	FA = 0
	VA = 0
	T = 0
	cai = 1e-8
}

BREAKPOINT {
	SOLVE state METHOD derivimplicit : see http://www.neuron.yale.edu/phpBB/viewtopic.php?f=28&t=592
}

LOCAL bfc , kfv

DERIVATIVE state {

	bfc = b * (Fmax-FA-VA) * cai^4
	kfv = k1 * FA * Ves

	cai'	= - bfc + 4 * u * FA
	FA'	= bfc - u * FA - kfv + k2 * VA
	VA'	= kfv - (k2+k3) * VA
	T'	= nt * k3 * VA - kh * T
}	

