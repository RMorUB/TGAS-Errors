c ---------------------------------------------------------------------
c                                	FOR TGAS-HIPPARCOS SUBSET
C**********************************************************************************************
c Subroutine that computes and assign TGAS errors for TYCHO-2 stars only to observables
c Errors are computed following the paper Michalik et al. 2015
c 
c Note: No correlation is assumed between astrometric errors. No dependency with scaning law 
c is assumed. The error for mu_alpha, mu_delta is imposed to be the same as mu as an 
c approximation.
c
c The code allows to compute the astrometric errors.
c
c Universitat de Barcelona
c Contact: Roger Mor, Merce Romero-Gomez
c email: rmor@am.ub.es
c
c This is the First Version and work is still in progress:
c Work is in progress to include the correlations
c Work is in progress to add errors as a function of G magnitude
c Future Work: Include dependencies on the scanning law
c
c ---------------------------------------------------------------------
c
c Updated: September 2015
c
c Input values:
c -------------
c     a(i): Astrometric true values for the star (not affected by errors)
c         a(1): Equatorial heliocentric right ascension (units: radians)
c         a(2): Equatorial heliocentric declination: delta (units: radians) 
c         a(3): Parallax: pi (units: mas)
c         a(4): Equatorial proper motions in right ascension in true arcs on the sky: mua_{*}=mua*cos(delta) (units: mas/yr)
c         a(5): Equatorial proper motions: mud (units: mas/yr)
c
c     V_T: V tycho magnitude. (TGAS is not going to contain new magnitudes for the stars)
c
c
c
c Output values: 
c -------------
c Parameters of the star expected in TGAS (affected by errors)
c TGAS errors assigned to each parameter:  
c         TGASER(1): Standard deviation in Equatorial heliocentric right ascension in true arcs on the sky:: alpha_{*}=alpha*cos(delta)  (units: mas)
c         TGASER(2): Standard deviation in Equatorial heliocentric declination: delta (units: mas) 
c         TGASER(3): Standard deviation in Parallax: pi (units: mas)
c         TGASER(4): Standard deviation in Equatorial proper motions in right ascension in true arcs on the sky: mua_{*}=mua*cos(delta) (units: mas/yr)
c         TGASER (5): Standard deviation in Equatorial proper motions: mud (units: mas/yr)
c   Astrometric data affected by errors  
c     ao(i): Astrometric values for the star affected by Gaia errors
c         TGASOB(1): Observed Equatorial heliocentric right ascension (units: radians)
c         TGASOB(2): Observed Equatorial heliocentric declination: delta (units: radians) 
c         TGASOB(3): Observed Parallax: pi (units: mas)
c         TGASOB(4): Observed Equatorial proper motions in right ascension in true arcs on the sky: mua_{*}=mua*cos(delta) (units: mas/yr)
c         TGASOB(5): Observed Equatorial proper motions: mud (units: mas/yr)

      Subroutine TGASHIPErrors(V_T,a,TGASOB,TGASER)
  
      Double Precision V_T, gasdev
      Double Precision a(5),TGASER(5), TGASOB(5)
      Integer Dseed

	INCLUDE 'const_math.h'
        INCLUDE 'const_ast.h'

      Dseed=time() ! number of seconds from 1970
      
c      Dseed=214749.d0 # Previous line can be substituted by this one if desired.







c       Computation of the observed astrometric quantities 

c       The error in right ascension TGASER(1) denotes true arc on the sky
c       so the right ascension shall be converted to that before the 
c       random error is assigned 
c       alpha_{*}=alpha*cos(delta) 
 


        a(1)=a(1)*dcos(a(2))


c  Conversion of (alpha*,delta) from radians to mas
 
	a(1)=a(1)/mas
	a(2)=a(2)/mas

	
c Here we assign the values for the standard deviation depending on magnitude.

      if(V_T.lt.7) then
		TGASER(1)=0.116
		TGASER(2)=0.116
		TGASER(3)= 0.180
		TGASER(4)= 0.017
		TGASER(5)= 0.017

		
      elseif(V_t.ge.7.and.V_T.lt.8) then
		TGASER(1)=0.120
		TGASER(2)=0.120
		TGASER(3)= 0.192
		TGASER(4)= 0.021
		TGASER(5)= 0.021
		 

      elseif(V_t.ge.8.and.V_T.lt.9) then 
		TGASER(1)=0.125
		TGASER(2)=0.125
		TGASER(3)= 0.198
		TGASER(4)= 0.029
		TGASER(5)= 0.029
      elseif(V_t.ge.9.and.V_T.lt.10) then
		TGASER(1)=0.133
		TGASER(2)=0.133
		TGASER(3)= 0.217
		TGASER(4)= 0.039
		TGASER(5)= 0.039
      elseif(V_t.ge.10.and.V_T.lt.11) then 
		TGASER(1)=0.154
		TGASER(2)=0.154
		TGASER(3)= 0.253
		TGASER(4)= 0.058
		TGASER(5)= 0.058

      elseif(V_t.ge.11.and.V_T.lt.12) then 
		TGASER(1)=0.128
		TGASER(2)=0.128
		TGASER(3)= 0.211
		TGASER(4)= 0.087
		TGASER(5)=0.087
      elseif(V_t.ge.12) then
		TGASER(1)=0.151
		TGASER(2)=0.151
		TGASER(3)= 0.248
		TGASER(4)= 0.135
		TGASER(5)= 0.135

      endif  

	do i=1,5
	
	TGASOB(i)=gasdev(Dseed, a(i),TGASER(i))

	enddo
	
      


c       The true and observed (alpha*,delta) are converted from mas to radiants 
	a(1)=a(1)*mas
	a(2)=a(2)*mas
	TGASOB(1)=TGASOB(1)*mas
	TGASOB(2)=TGASOB(2)*mas


c       and the alpha_{*} is converted to alpha
        TGASOB(1)=TGASOB(1)/dcos(a(2))
c       we also convert back the input value
        a(1)=a(1)/dcos(a(2))


        





      end
