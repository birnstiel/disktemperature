subroutine OPACITY(model,top,shape,ross,rho,nt,t0,t1,output)
	IMPLICIT NONE
	INTEGER IT,NT      
	REAL*8 eD, eG, rho, T, T0, T1, dT, aKext
	DIMENSION eD(5,6) 
	DIMENSION eG(71,71)
	CHARACTER*3 model  
	CHARACTER*1 top    
	CHARACTER*1 shape  
	LOGICAL ross

	REAL*8,intent(out) :: OUTPUT(2,NT)


	! Initialization of all necessary data for a chosen dust model:

    call init_d(ross,model,top,shape,eD)
    
    ! Initialization of data for the gas model:
    
    call init_g(ross,eG)
    
	!  
    !  Start loop by temperatures:
    !
	T = T0
	dT = DEXP(DLOG(T1/T0)/(NT-1))
	
	DO IT = 1, NT
		!-----------------------------------------------------------------------------
		! Calculation of Rosseland or Planck mean opacities using a chosen density, 
		! temperature, silicate dust model, shape of the grains and their topology:
		!-----------------------------------------------------------------------------
		CALL cop(eD,eG,rho,T,aKext)
		! Write results to the output: 
		OUTPUT(1,IT)=T
		OUTPUT(2,IT)=aKext
		!         WRITE(99,'(2X,F10.4,1X,E15.6)') T, aKext
		! Jump to the next temperature:
		T = T * dT         
      END DO

	! Exit:        
END