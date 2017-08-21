!==========================================================================
! THE OMEXDIA model, implemented in FORTRAN
! Karline Soetaert, nioz-yerseke
!==========================================================================

!==========================================================================
! initialise the common block with parameter values, 
! followed by thicknesses, porosities, bioturbation values
!==========================================================================

      SUBROUTINE initomexdia (steadyparms)
      IMPLICIT NONE
      EXTERNAL steadyparms

      INTEGER,PARAMETER :: N=100
      INTEGER,PARAMETER :: nc = 2*N + 3*(N+1) + 26

      DOUBLE PRECISION parms(nc)
      COMMON /myparms/parms

       CALL steadyparms(nc, parms)
       
      RETURN
      END SUBROUTINE

!==========================================================================
! Initialise the forcing function common block (Carbon flux)
!==========================================================================

      SUBROUTINE initforc (steadyforcs)
      IMPLICIT NONE
      EXTERNAL steadyforcs

      INTEGER,PARAMETER :: N=1 

      DOUBLE PRECISION forcs(N)
      COMMON /myforcs/forcs

       CALL steadyforcs(N, forcs)
       
      RETURN
      END SUBROUTINE

!==========================================================================
!==========================================================================
! subroutine calculating the rate of change of
! the omexdia model - here the quantities in the common blocks are named
!==========================================================================
!==========================================================================

      SUBROUTINE omexdiamod (neq, t, Conc, dConc, yout, ip)
      IMPLICIT NONE

!......................... declaration section.............................
      INTEGER           :: neq, ip(*), i
      INTEGER,PARAMETER :: N=100         

      DOUBLE PRECISION  :: t, Conc(6*N), dConc(6*N), yout(*)
      DOUBLE PRECISION  :: por(N),intpor(N+1),Db(N+1),dx(N),dxInt(N+1)
      DOUBLE PRECISION  :: Fdet(N),Sdet(N),O2(N),NO3(N),NH3(N),ODU(N)
      DOUBLE PRECISION  :: dFdet(N),dSdet(N),dO2(N),dNO3(N),                    &
     &                     dNH3(N),dODU(N)

      DOUBLE PRECISION  :: cflux,Cprod(N),Nprod(N),Rescale(N),TOC(N)
      DOUBLE PRECISION  :: Oxicminlim(N),Denitrilim(N),Anoxiclim(N) 
      DOUBLE PRECISION  :: Oxicmin(N),Denitrific(N),anoxicmin(N)
      DOUBLE PRECISION  :: nitri(N),oduox(N),odudepo(N)
      DOUBLE PRECISION  :: pdepo, Sum
      DOUBLE PRECISION  :: Flux(N+1), DS(N+1)
       
      DOUBLE PRECISION  :: MeanFlux,rFast,rSlow ,pFast,w,NCrFdet,               & 
     &  NCrSdet,bwO2,bwNO3,bwNH3,bwODU,NH3Ads,rnit,ksO2nitri,                   &
     &  rODUox,ksO2oduox,ksO2oxic,ksNO3denit,kinO2denit,kinNO3anox,             &
     &  kinO2anox,DispO2,DispNO3,DispNH3,DispODU,TOC0   

      COMMON /myparms     /MeanFlux,rFast,rSlow ,pFast,w,NCrFdet,               & 
     &  NCrSdet,bwO2,bwNO3,bwNH3,bwODU,NH3Ads,rnit,ksO2nitri,                   &
     &  rODUox,ksO2oduox,ksO2oxic,ksNO3denit,kinO2denit,kinNO3anox,             &
     &  kinO2anox,DispO2,DispNO3,DispNH3,DispODU,TOC0,dx,dxint,                 &
     &  por,intpor,Db   

      DOUBLE PRECISION CarbonFlux
      COMMON /myforcs/ CarbonFlux

! output variables
      DOUBLE PRECISION  :: O2flux, O2deepflux, NO3flux, NO3deepflux
      DOUBLE PRECISION  :: NH3flux,NH3deepflux,ODUflux,ODUdeepflux
      DOUBLE PRECISION  :: partDenit, partAnoxic,partOxic
      COMMON /myout       /O2flux, O2deepflux, NO3flux, NO3deepflux,            &
     &  NH3flux,NH3deepflux,ODUflux, ODUdeepflux, Cflux, partDenit,             &
     &  partAnoxic,partOxic, Cprod,Nprod,                                       &         
     &  TOC, Oxicmin,Denitrific,anoxicmin,nitri,oduox,odudepo

      CHARACTER(len=80) msg
!............................ statements ..................................

!     check memory allocated to output variables
      IF (ip(1) < 912)  CALL rexit("nout should be at least 912") 

! from Conc to fdet, sdet, o2,...
      DO I = 1, N
        Fdet(I) = Conc(I)
        Sdet(I) = Conc(N+I)
        O2(I)   = Conc(2*N+I)
        NO3(I)  = Conc(3*N+I)
        NH3(I)  = Conc(4*N+I)
        ODU(I)  = Conc(5*N+I)
      ENDDO
      
      TOC = (FDET + SDET)*1200d0*1e-9/2.5 + TOC0

! --------------------------------------------------------------------------
! Rate of change due to transport 
! --------------------------------------------------------------------------
       
! Carbon deposition flux
      CFlux  =  CarbonFlux
       
      CALL diff1D (N, FDET, 0.d0, 0.d0, cFlux*pFast, 0.d0, 0.d0, 0.d0,          &
     &             1, 3, Db, 1.d0-intpor, dx, dxint, Flux, dFDET)

      CALL diff1D (N, SDET, 0.d0,0.d0,cFlux*(1.d0-pFast),0.d0,0.d0,0.d0,        &
     &             1, 3, Db,1.d0-intpor, dx, dxint, Flux, dSDET)


      Ds = dispO2 * intpor  ! effective diffusion coefficient 
      CALL diff1D (N, O2 ,bwO2 ,0.d0, 0.d0, 0.d0, 0.d0, 0.d0,                   &
     &             2,  3, Ds ,intpor, dx, dxint, Flux, dO2)     
      O2flux      = Flux(1)
      O2deepflux  = Flux(N+1)
                        
      Ds  = dispNO3*intpor 
      CALL diff1D (N, NO3 ,bwNO3 ,0.d0, 0.d0, 0.d0, 0.d0, 0.d0,                 &
     &             2,  3, Ds ,intpor, dx, dxint, Flux, dNO3)     
      NO3flux     = Flux(1)
      NO3deepflux = Flux(N+1)

      Ds = dispNH3/(1+NH3Ads)*intpor
      CALL diff1D (N, NH3 ,bwNH3 ,0.d0, 0.d0, 0.d0, 0.d0, 0.d0,                 &
     &             2,  3, Ds ,intpor, dx, dxint, Flux, dNH3)     
      NH3flux     = Flux(1)
      NH3deepflux = Flux(N+1)

      Ds = dispODU*intpor
      CALL diff1D (N, ODU ,bwODU ,0.d0, 0.d0, 0.d0, 0.d0, 0.d0,                 &
     &             2,  3, Ds ,intpor, dx,dxint, Flux, dODU)     
      ODUflux     = Flux(1)
      ODUdeepflux = Flux(N+1)

! --------------------------------------------------------------------------
! Rate of change due to biogeochemistry 
! --------------------------------------------------------------------------

! Production of DIC and DIN, expressed per cm3 LIQUID/day

      Cprod= (rFast*FDET        +rSlow*SDET        ) * (1.d0-por)/por
      Nprod= (rFast*FDET*NCrFdet+rSlow*SDET*NCrSdet) * (1.d0-por)/por

! Oxic mineralisation, denitrification, anoxic mineralisation

! first the limitation terms
      Oxicminlim = O2/(O2+ksO2oxic)                ! limitation terms
      Denitrilim = (1.d0-O2/(O2+kinO2denit)) * NO3/(NO3+ksNO3denit)
      Anoxiclim  = (1.d0-O2/(O2+kinO2anox)) * (1-NO3/(NO3+kinNO3anox))
      Rescale    = 1.d0/(Oxicminlim+Denitrilim+Anoxiclim)

! then the mineralisation rates
      OxicMin    = Cprod*Oxicminlim*Rescale        ! oxic mineralisation
      Denitrific = Cprod*Denitrilim*Rescale        ! Denitrification
      AnoxicMin  = Cprod*Anoxiclim *Rescale        ! anoxic mineralisation

! reoxidation and ODU deposition
      Nitri      = rnit  *NH3*O2/(O2+ksO2nitri)
      OduOx      = rODUox*ODU*O2/(O2+ksO2oduox)

      pDepo      = MIN(1.d0,0.233*(w*365)**0.336 )
      OduDepo    = AnoxicMin*pDepo 

! --------------------------------------------------------------------------
! Update the rate of change with rates due to biogeochemical processes
! --------------------------------------------------------------------------
      dFDET = dFDET  - rFast*FDET
      dSDET = dSDET  - rSlow*SDET
      dO2   = dO2    - OxicMin          -2.d0* Nitri -      OduOx
      dNH3  = dNH3   + (Nprod                - Nitri) / (1.d0+NH3Ads)
      dNO3  = dNO3   - 0.8d0*Denitrific      + Nitri 
      dODU  = dODU   + AnoxicMin                  - OduOx - OduDepo

! from dfdet, dsdet, do2,... to dconc, and calculate integrated rates
      partDenit = 0.D0
      partOxic = 0.D0
      partAnoxic = 0.D0
      DO I = 1, N
         partDenit  = partDenit + Denitrific(I) * por(I)*dx(I)
         partAnoxic = partAnoxic + AnoxicMin(I) * por(I)*dx(I)
         partOxic   = partOxic + OxicMin(I) * por(I)*dx(I)
         dConc(I)      =  dFdet(I)
         dConc(N+I)    =  dSdet(I) 
         dConc(2*N+I)  =  dO2(I)  
         dConc(3*N+I)  =  dNO3(I) 
         dConc(4*N+I)  =  dNH3(I) 
         dConc(5*N+I)  =  dODU(I) 
      ENDDO
      Sum = partDenit+ partAnoxic + partOxic
      partDenit = partDenit / Sum
      partOxic = partOxic / Sum
      partAnoxic = partAnoxic / Sum
      CALL getout(yout)
      
      RETURN
      END SUBROUTINE

!==========================================================================
! put output variables in one vector
!==========================================================================

      SUBROUTINE getout(yout)
      DOUBLE PRECISION :: yout(*), out(912)
      INTEGER :: i

      COMMON /myout    /out
      DO i = 1, 912
       yout(i) = out (i)
      ENDDO       
 
      END SUBROUTINE getout
       
 
!==============================================================================
! Diffusion in a 1-dimensional finite difference grid 
! all inputs are vectors
! subroutine from ReacTran in isnt\doc\fortran directory
!==============================================================================

      SUBROUTINE diff1d (N, C, Cup, Cdown, fluxup, fluxdown, aup, adown,        &
     &                   BcUp, BcDown, D, VF, dx, dxaux,                        &
     &                   Flux, dC) 
      IMPLICIT NONE
      INTEGER N                  ! length of C
C input
      DOUBLE PRECISION C(N)

C Boundary concentrations (used if Bc..=2,4), fluxes (used if Bc= 1) 
C and convection coeff (used if Bc=4)
      DOUBLE PRECISION Cup, Cdown, fluxup, fluxdown, aup, adown

C Diffusion, volume fraction
      DOUBLE PRECISION D(N+1), VF(N+1)

C grid size, distance from mid to mid
      DOUBLE PRECISION dx(N), dxaux(N+1)

C boundary concitions (1= flux, 2=conc, 3 = 0-grad, 4=convect)
      INTEGER BcUp, BcDown   

C output: fluxes and rate of change
      DOUBLE PRECISION Flux(N+1), dC(N)

C locals 
      INTEGER I
      DOUBLE PRECISION AVF

C -------------------------------------------------------------------------------

C Flux - first internal cells

      DO I = 2,N
        Flux(I) = -VF(I)*D(I) * (C(I)-C(I-1)) /dxaux(I)
      ENDDO

C Then the outer cells 
C upstream boundary
      IF (BcUp .EQ. 1) THEN
        Flux(1) = fluxup

      ELSE IF (BcUp .EQ. 2) THEN
        Flux(1) = -VF(1)*D(1) * (C(1)-Cup) /dxaux(1)

      ELSE IF (BcUp .EQ. 3) THEN
        Flux(1) = 0.D0

      ELSE 
        Flux(1) = aup * (Cup - C(1))

      ENDIF

C downstream boundary
      IF (BcDown .EQ. 1) THEN
        Flux(N+1) = fluxdown

      ELSE IF (BcDown .EQ. 2) THEN
        Flux(N+1) = -VF(N+1)*D(N+1) * (Cdown-C(N)) /dxaux(N+1)

      ELSE IF (BcDown .EQ. 3) THEN
        Flux(N+1) =0.D0

      ELSE
        Flux(N+1) = -adown * (Cdown-C(N))

      ENDIF


C Rate of change = negative flux gradient
      DO I = 1,N
        AVF   = 0.5 * (VF(I)+VF(I+1))
        dC(I) = -(Flux(I+1) - Flux(I)) / AVF / dx(I)
      ENDDO
    
      RETURN
      END SUBROUTINE diff1D
