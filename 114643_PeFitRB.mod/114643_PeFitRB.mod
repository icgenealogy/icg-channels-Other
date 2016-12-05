://=============================================================================
://  eFitRB  -  electrode fitting via white noise injection
://		(procedure by Romain Brette)
://=============================================================================
://
:// REMARKS:
://	- automatic full electrode compensation
:// 	- noise injection is used to determine the linear response 
://	  of the electrode
://
://	- Finds the kernel by least square minimization 
://	  (Wiener-Hopf / Levinson algorithm)
://
:// In this file:
://	- we use uniform random numbers instead of normal ones
://	- we use the full Wiener-Hopf linear optimization method
://	- we find tau by least square fitting (exponential)
://	- we use golden section search to find R and Re
://
://=============================================================================
://
:// Romain Brette - Last update: 12 Dec 2005
:// brette@di.ens.fr
://
://=============================================================================

:// Set the kernel filename here

VERBATIM
#define kernel_filename		"kernel.txt"
#define rawkernel_filename	"rawkernel.txt"
ENDVERBATIM

TITLE Automatic full electrode compensation

NEURON {
	POINT_PROCESS eFitRB
	POINTER vcell		:// membrane potential
	RANGE active
	RANGE amp
	RANGE ksize, tail
	RANGE t_start, t_stop
	RANGE v0, RplusRe
	RANGE R, tau, Re
	RANGE n
	USEION z READ ez WRITE iz VALENCE 1
}

UNITS {
	(MOhm) = (megaohm)
	(mV) = (millivolt)
}

PARAMETER {
	dt		(ms)
	amp = 0.5	(nA)	:// maximum value of the noise
	t_start = 300	(ms)	:// starting time for the compensation (injection starts at time 0)
	t_stop = 19000	(ms)	:// stopping time for the compensation procedure

	ksize = 100		:// size of the kernel (in time steps)
	tail = 50		:// tail starts at this timestep

	active = 0		:// On/Off switch
}

ASSIGNED {
	iz		(nA)
	ez		(mV)
	n
	vcell		(mV)
	v0		(mV)
	i0		(nA)

	Re		(MOhm)	:// electrode resistance
	R		(MOhm)	:// neuron input resistance
	tau		(ms)	:// membrane time constant

	VI[500]
	II[500]
	previous_i[500]	(nA)
	tkern[500]
	yk[500]

	kernel[500]
	r[500]
	g[500]
	h[500]

	saved

	invVar			:// inverse of noise variance

	invN
}

INITIAL {
	iz = 0
	n = 0
	v0 = vcell
	FROM k = 0 TO 499 {
		VI[k] = 0
		II[k] = 0
		previous_i[k] = 0
	}

	saved = 0

	invVar = (1 / (amp*amp))* 3
	i0=0
	invN=0

	R=0
	Re=0
	tau=0
}

BREAKPOINT {
	SOLVE injection
}

PROCEDURE injection() {
	if (active == 1) {	:// recalculate the current every 'stepsize' time steps
        
		:// Shift the history of injected current
		VERBATIM
		{
			int k, kmax;
			double *pi;

			kmax = (int) ksize;
			pi = &(previous_i[kmax-1]);

			for(k = 1;k<kmax;k++) {
				*pi=*(pi-1);
				pi--;
			}
			*previous_i = -iz;
		}
		ENDVERBATIM
		://NMODL EQUIVALENT
		://----------------
		://FROM k = 1 TO ksize-1 {
		://	previous_i[ksize-k] = previous_i[ksize - k -1 ]
		://}
		://previous_i[0] = -iz

		if (t>t_start && t <t_stop) {	:// start the procedure at tstart
			n = n + 1
			invN=1/n		:// multiplications are faster than divisions
		
			:// Calculate v0 = mean v (online averaging)
			v0 = v0 +(vcell - v0)*invN

			:// Calculate i0 = mean i (online averaging)
			i0 = i0 +(previous_i[0] - i0)*invN

			:// Calculate the V-I correlations and I autocorrelation
			VERBATIM
			{
				int k, kmax;
				double *pi;
				double *VIloc;
				double *IIloc;

				kmax = (int) ksize;
				pi = previous_i;
				VIloc = VI;
				IIloc = II;

				for(k = 0;k<kmax;k++) {
					*VIloc+= (vcell*(*pi)*invVar - (*VIloc))*invN;
					*IIloc+= ((*previous_i)*(*pi)*invVar-(*IIloc))*invN;

					pi++;
					VIloc++;
					IIloc++;
				}
			}
			ENDVERBATIM

			:// NMODL EQUIVALENT
			:// ----------------
			://FROM k = 0 TO ksize-1 {
			://	kernel[k] = kernel[k] + ((vcell-v0)*previous_i[k]/(sigma*sigma)-kernel[k])/n
			://	Re = Re + kernel[k] - pass_kernel[k]
			://}
			:// vr = vcell
			://
			://FROM k = 0 TO ksize-1 {
			://	true_kernel[k] = (kernel[k] - pass_kernel[k])*Re/kern_int
			://	vr = vr - true_kernel[k]*previous_i[k]
			://}
		}

		if (t<t_stop) {
			:// Inject white noise
			iz = -(2*scop_random()-1)*amp
		} else {
			iz = 0
			:// Save the kernel at the end of the simulation
			if (saved == 0) {
				saved = 1
				save_kernel()
			}
		}
	}
}

:// Do the calculations and
:// save the kernel to a file
PROCEDURE save_kernel() {
	:// Adjust VI
	VERBATIM
	{
		int k;

		for(k=0;k<(int)ksize;k++) {
			VI[k]-= v0*i0*invVar;
		}
	}
	ENDVERBATIM

	VERBATIM
	{

		FILE *f_write;
		int k,nel;
		double alpha,lambda,tot;

		if (f_write=fopen("VI.txt","w")) {
			for(k=0;k<(int)ksize;k++) {
				fprintf(f_write,"%lf\n",VI[k]);
			}
			fclose(f_write);
		}

		if (f_write=fopen("II.txt","w")) {
			for(k=0;k<(int)ksize;k++) {
				fprintf(f_write,"%lf\n",II[k]);
			}
			fclose(f_write);
		}

		/* The autocorrelation matrix MI is formed from II
			with MI(i,j)=II(|i-j|)		(C convention for arrays)
		   MI is a symetric Toeplitz matrix
		   We use Levinso's method (O(N^2))
		 	to solve MI*K=VI */
		for(k=0;k<(int)ksize;k++) {
			r[(int)ksize-1+k]=II[k];
			r[(int)ksize-1-k]=II[k];
		}
		nel=(int)ksize;

		/* Levinsons method */
		{
			int j,k,m,m1,m2;
			double pp,pt1,pt2,qq,qt1,qt2,sd,sgd,sgn,shn,sxn;

			for(k=0;k<(int)ksize;k++) {
				g[k]=0;
				h[k]=0;
				kernel[k]=0;
			}

			kernel[0]=VI[0]/r[nel-1];
			g[0]=r[nel-2]/r[nel-1];
			h[0]=r[nel]/r[nel-1];
			for (m=0;m<nel;m++) {
				m1=m+1;
				sxn = -VI[m1];
				sd = -r[nel-1];
				for (j=0;j<m;j++) {
					sxn += r[nel-1+m1-j]*kernel[j];
					sd += r[nel-1+m1-j]*g[m-j];
				}
				kernel[m1]=sxn/sd;
				for (j=0;j<m;j++) kernel[j] -= kernel[m1]*g[m-j];
				if (m1 == nel-1) break;
				sgn = -r[nel-1-m1];
				shn = -r[nel-1+m1];
				sgd = -r[nel-1];
				for (j=0;j<m;j++) {
					sgn += r[nel-1+j-m1]*g[j];
					shn += r[nel-1+m1-j]*h[j];
					sgd += r[nel-1+j-m1]*h[m-j];
				}
				g[m1]=sgn/sgd;
				h[m1]=shn/sd;
				k=m;
				m2=(m+1) >> 1;
				pp=g[m1];
				qq=h[m1];
				for (j=0;j<m2;j++) {
					pt1=g[j];
					pt2=g[k];
					qt1=h[j];
					qt2=h[k];
					g[j]=pt1-pp*qt2;
					g[k]=pt2-pp*qt1;
					h[j]=qt1-qq*pt2;
					h[k--]=qt2-qq*pt1;
				}
			}
		}

		kernel[0]=0.;
		kernel[1]=0.;

		if (f_write=fopen("rawkernel.txt","w")) {
			for(k=0;k<(int)ksize;k++) {
				fprintf(f_write,"%lf\n",kernel[k]);
			}
			fclose(f_write);
		}

		/* First guess for R, Re and tau (fit the tail of the kernel) */
		/* exponential fit */
		{
			double x,y, sx2y, sylny, sxy, sxylny, sy;
			int k;

			sylny=0;
			sx2y=0;
			sxy=0;
			sxylny=0;
			sy=0;
			for(k=tail;k<(int)ksize;k++) {
				y = kernel[k];
				x = k*dt;
				if (y > 0.) {
					sylny = sylny + y * log (y);
					sx2y = sx2y + x * x * y;
					sxy = sxy + x * y;
					sxylny = sxylny + x * y * log(y);
					sy = sy + y;
				}
			}

			R = exp((sx2y * sylny - sxy * sxylny) / (sy * sx2y - sxy*sxy));
			tau = - (sy * sx2y - sxy*sxy) / (sy * sxylny - sxy * sylny);
		}

		R=R*tau/dt;
		/* Find Re */
		Re=0;
		for(k=0;k<(int)ksize;k++) {
			Re=Re+kernel[k]-(R*dt/tau)*exp(-k*dt/tau);
		}

		/* fprintf(stderr,"First estimation:\n");
		   fprintf(stderr,"Re = %f\n",Re);
		   fprintf(stderr,"R = %f\n",R);
		   fprintf(stderr,"tau = %f\n",tau); */

		/* Remove the membrane kernel */
		/* Find R and Re */
		/* Golden section search */
		{
			static double f(double x) {
				/* This function removes the membrane kernel 
				   and estimates the error on the tail */
				double tailerr=0.;

				alpha = x*dt/tau;
				lambda = exp(-dt/tau);
				tot=0.;
				yk[0] = alpha/(alpha+1)*kernel[0];
				tkern[0] = kernel[0]-yk[0];
				tot+=tkern[0];
				for(k=1;k<(int)ksize;k++) {
					yk[k] = (alpha*kernel[k]+lambda*yk[k-1])/(alpha+1.);
					tkern[k]=kernel[k]-yk[k];
					tot+=tkern[k];
				}

				/* Minimize the tail */
				for(k=tail;k<(int)ksize;k++)
					tailerr+=tkern[k]*tkern[k];

				return tailerr;
			}

			double taugold=0.5*(sqrt(5.)-1.);
			double a=.25*R/Re;
			double b=4.*R/Re;
			double x1=a+(1.-taugold)*(b-a);
			double x2=a+taugold*(b-a);
			double f1=f(x1);
			double f2=f(x2);

			while((b-a)>0.0001) {
				if (f1>f2) {
					a=x1;
					x1=x2;
					f1=f2;
					x2=a+taugold*(b-a);
					f2=f(x2);
				} else {
					b=x2;
					x2=x1;
					f2=f1;
					x1=a+(1.-taugold)*(b-a);
					f1=f(x1);
				}
			}

			Re=tot;
			R=a*Re;	
		}
		/* fprintf(stderr,"R=%f MOhm\tRe=%f MOhm\n",R,Re); */

		/* Write the file (up to the tail) */
		if (f_write=fopen(kernel_filename,"w")) {
			for(k=0;k<(int)tail;k++)
				fprintf(f_write,"%lf\n",tkern[k]);
			fclose(f_write);
		}

	}
	ENDVERBATIM
}

://=============================================================================
