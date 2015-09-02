c ---------------------------------------------------------------------
c
c Subroutine that computes and assign TGAS errors for TYCHO-2 only subset and HIPPARCOS subset
c to observables. Errors are computed following the paper Michalik et al. 2015
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
c*********************************************************************
c     subset: INTEGER Variable to select the TGAS subset             *
c--------------------------------------------------------------------*
c                                                                    *     
c            subset=1 : For Tycho Only subset of TGAS                *
c            subset=2 : For HIPPARCOS subset of TGAS                 * 
c*********************************************************************
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

      Subroutine TGASErrors(V_T,a,TGASOB,TGASER, subset)
  
      Double Precision V_T
      Double Precision a(5),TGASER(5), TGASOB(5)
      Integer Dseed, subset

      INCLUDE 'const_math.h'
      INCLUDE 'const_ast.h'

      Dseed=time() ! number of seconds from 1970 (If you prefere a fixed seed use the line below)
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

	
c Here we assign the values for the standard deviation depending on magnitude and on subset.
      if(subset.eq.1) then
      	if(V_T.lt.7) then
		TGASER(1)=0.244d0
		TGASER(2)=0.244d0
		TGASER(3)= 0.399d0
		TGASER(4)= 0.198d0
		TGASER(5)= 0.198d0

		
      	elseif(V_t.ge.7.and.V_T.lt.8) then
		TGASER(1)=0.198d0
		TGASER(2)=0.198d0
		TGASER(3)= 0.348d0
		TGASER(4)= 0.264d0
		TGASER(5)= 0.264d0
		 

      	elseif(V_t.ge.8.and.V_T.lt.9) then 
		TGASER(1)=0.191d0
		TGASER(2)=0.191d0
		TGASER(3)= 0.327d0
		TGASER(4)= 0.403d0
		TGASER(5)= 0.403d0
      	elseif(V_t.ge.9.and.V_T.lt.10) then
		TGASER(1)=0.230d0
		TGASER(2)=0.230d0
		TGASER(3)= 0.407d0
		TGASER(4)= 0.680d0
		TGASER(5)= 0.680d0
      	elseif(V_t.ge.10.and.V_T.lt.11) then 
		TGASER(1)=0.329d0
		TGASER(2)=0.329d0
		TGASER(3)= 0.601d0
		TGASER(4)= 1.145d0
		TGASER(5)= 1.145d0

      	elseif(V_t.ge.11.and.V_T.lt.12) then 
		TGASER(1)=0.379d0
		TGASER(2)=0.379d0
		TGASER(3)= 0.722d0
		TGASER(4)=1.522d0
		TGASER(5)=1.522d0
      	elseif(V_t.ge.12) then
		TGASER(1)=0.349d0
		TGASER(2)=0.349d0
		TGASER(3)= 0.702d0
		TGASER(4)= 1.615d0
		TGASER(5)= 1.615d0

      	endif  

      elseif(subset.eq.2) then

      	if(V_T.lt.7) then
		TGASER(1)=0.116d0
		TGASER(2)=0.116d0
		TGASER(3)= 0.180d0
		TGASER(4)= 0.017d0
		TGASER(5)= 0.017d0

		
      	elseif(V_t.ge.7.and.V_T.lt.8) then
		TGASER(1)=0.120d0
		TGASER(2)=0.120d0
		TGASER(3)= 0.192d0
		TGASER(4)= 0.021d0
		TGASER(5)= 0.021d0
		 

      	elseif(V_t.ge.8.and.V_T.lt.9) then 
		TGASER(1)=0.125d0
		TGASER(2)=0.125d0
		TGASER(3)= 0.198d0
		TGASER(4)= 0.029d0
		TGASER(5)= 0.029d0
      	elseif(V_t.ge.9.and.V_T.lt.10) then
		TGASER(1)=0.133d0
		TGASER(2)=0.133d0
		TGASER(3)= 0.217d0
		TGASER(4)= 0.039d0
		TGASER(5)= 0.039d0
      	elseif(V_t.ge.10.and.V_T.lt.11) then 
		TGASER(1)=0.154d0
		TGASER(2)=0.154d0
		TGASER(3)= 0.253d0
		TGASER(4)= 0.058d0
		TGASER(5)= 0.058d0

      	elseif(V_t.ge.11.and.V_T.lt.12) then 
		TGASER(1)=0.128d0
		TGASER(2)=0.128d0
		TGASER(3)= 0.211d0
		TGASER(4)= 0.087d0
		TGASER(5)=0.087d0
      	elseif(V_t.ge.12) then
		TGASER(1)=0.151d0
		TGASER(2)=0.151d0
		TGASER(3)= 0.248d0
		TGASER(4)= 0.135d0
		TGASER(5)= 0.135d0

      	endif  

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
