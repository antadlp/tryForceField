

* ********************************************************
*           Interactive Molecular Dynamics                *
*               Module: pressure cooker                   *
*                         by                              *
*                Dr. Ruben Santamaria                     *
*                 rso@fisica.unam.mx                      *
*            Dept. of Theoretical Physics                 *
*                Institute of Physics                     *
*               Univ. of Mexico, UNAM                     *
* ********************************************************
      program pressure_cooker
      implicit double precision (a-h,o-z)
      include 'parametros'

      character atsymbol(nmax)*5
      dimension atomass(nmax)

* positions

      dimension x(nmax),    y(nmax),    z(nmax)
      dimension xnew(nmax), ynew(nmax), znew(nmax)
      dimension Xn(nmax),   Yn(nmax),   Zn(nmax)
      dimension x0(nmax),   y0(nmax),   z0(nmax)

* speeds

      dimension vx(nmax),    vy(nmax),    vz(nmax)
      dimension vxnew(nmax), vynew(nmax), vznew(nmax)
      dimension Vxn(nmax),   Vyn(nmax),   Vzn(nmax)

* ab-initio gradients

      dimension Gradx(nmax),    Grady(nmax),    Gradz(nmax)
      dimension Gradxnew(nmax), Gradynew(nmax), Gradznew(nmax)

* kinetic & potential energies

      dimension enekin(-1:maxt), enepot(-1:maxt)

* temperatures, pressure, force constants, etc

      dimension tCage(maxt), tConf(maxt), press(maxt), sK(nmax)
      dimension ftest(nmax), sum(nmax), knt(m1)
      dimension zK(nmax)

      integer seed
      real*8 zbqlu01

* ------------------------------------------
* identify initial cpu processes & directories
* ------------------------------------------

      call infoUser

      print*, "* ---------------------------------------------- *"
      print*, "          INTERACTIVE MOLECULAR DYNAMICS          "
      print*, "             Module:  Pressure Cooker"
      print*
      print*, "  CITATION:"
      print*, "  IMD program by Ruben Santamaria, IFUNAM, 2014"
      print*, "* ---------------------------------------------- *"

      print*, "internal modules (alphabetical order)"
      print*, "analytic"
      print*, "classify"
      print*, "gauss"
      print*, "main-olla"
      print*, "positions"
      print*, "powerSer"
      print*, "printGauss"
      print*, "speeds"
      print*
      print*, "external modules"
      print*, "abEne"
      print*, "averageT"
      print*, "breakBond"
      print*, "cartesian"
      print*, "citation"
      print*, "distances"
      print*, "eneSpring"
      print*, "estimateCPU"
      print*, "figures"
      print*, "groupCharge"
      print*, "infoUser"
      print*, "initConditOlla"
      print*, "inptTera"
      print*, "inptNwchem"
      print*, "kinEne"
      print*, "kparse"
      print*, "largest"
      print*, "makeMovie"
      print*, "massCenter"
      print*, "molInfo"
      print*, "periodicTable"
      print*, "perturb"
      print*, "polar"
      print*, "restartOlla"
      print*, "setWFs"
      print*, "scriptAB"
      print*, "shrink"
      print*, "slowSpeeds"
      print*, "taylor"
      print*, "units"
      print*, "workDir"
      print*, "writeMol"
      print*, "zRndGenerator"

* -------------------------
* collect molecular info
* -------------------------

      call molInfo (npar,ncage, atsymbol,atomass, x,y,z, dtIni,dt,Time,
     &     rfac,vfac, tempEquil, sK,xi,radiusf,ifreqShr,ifreqStor,
     &     irestDyn,iwavFunct,nsteps,ibegin,trest,ngrp1,ngrp2,ngrp3,
     &     itest,imethod)

     
      !open (unit=32, file= 'verERRORES.out',status= 'new')

* npar = numero total de atomos

* print input to the output file & create
* the script to run the ab-initio calculation

      !MODIFICACION RESORTES: COMENTAREAR

      !if (imethod.eq.1) call inptTera (3)
      !if (imethod.eq.2) call inptNwchem (3,npar,atsymbol,x,y,z)
      !call scriptAB (imethod)

      !MODIFICACION RESORTES: COMENTAREAR END

      call zbqlini (seed)


      !MODIFICACION RESORTES: COMENTAREAR
* restart the dynamics


!      if (irestDyn.eq.2) then ! <----
!                                     !
!      call restartOlla (2,ibegin,dt,trest,npar,ncage,x,y,z,vx,vy,vz,
!     & Gradx,Grady,Gradz, x0,y0,z0,sK,refPE,ibTime1,ibTime2,ibTime3)
!                                     !
!      write(*,*) "restarting the dynamics from:"
!      write(*,"(a,i6,1x,f9.2)") "step , time", ibegin, trest
!                                     !
!                                     ! reuse the wave function
!      if (imethod.eq.1) then         !
!      call inptTera (2)              !
!      call writeMol (npar,atsymbol,x,y,z)
!      end if                         !
!                                     !
!      if (imethod.eq.2) then         !
!      call inptNwchem (2,npar,atsymbol,x,y,z)
!      end if                         !
!                                     !
!      goto 120                       !
!                                     !
!      end if ! <---------------------
!
* -------------------------------
* initial conditions on positions & speeds
* -------------------------------


      !MODIFICACION RESORTES: COMENTAREAR END
      
      imzero= 0

      call massCenter (ncage,imzero,atomass,x,y,z,cmx,cmy,cmz)

      call initConditOlla (npar,ncage, atomass, x,y,z, vx,vy,vz,
     &    x0,y0,z0,sK, cmx,cmy,cmz, tempEquil,
     &    nsteps,ibTime1,ibTime2,ibTime3)

* ----------------------------------------------
* write the molec, submit a job & get new forces
* ----------------------------------------------

      call perturb (npar,ncage, x,y,z, vx,vy,vz,
     &    xnew,ynew,znew, vxnew,vynew,vznew, rfac,vfac)

      tt= -dt  ! present time

      write (*,*) "--------------------"
      write (*,"(1x,a,i6,1x,f9.2)") "step & time", -1, tt

      !MODIFICACION RESORTES: COMENTAREAR
!* create the external input file
!
!      if (imethod.eq.1) then
!      call inptTera (iwavFunct)
!      call writeMol (npar,atsymbol,xnew,ynew,znew)
!      end if
!
!      if (imethod.eq.2) then
!      call inptNwchem (iwavFunct,npar,atsymbol,xnew,ynew,znew)
!      end if
!
!* submit a job, give seconds to close files & get forces
!
!      if (itest.eq.2) then ! <----------
!                                        !
!      if (imethod.eq.1) call system ("./script.tera")
!      if (imethod.eq.2) call system ("./script.nwchem")
!                                        !
!      call sleep (02)                   !

!
!      call abEne (imethod,dtIni,dt,tt,npar,
!     &  Gradxnew,Gradynew,Gradznew,Epot,refPE)
!                                        !
!      if (imethod.eq.1) call groupCharge (tt,ngrp1,ngrp2,ngrp3)
!                                        !
!      end if ! <------------------------


!MODIFICACION RESORTES: COMENTAREAR END

!MODIFICACION RESORTES: CODIGO NUEVO      
!      call abEne (imethod,dtIni,dt,tt,npar,
!     &  Gradxnew,Gradynew,Gradznew,Epot,refPE)


      call abEneR (imethod,dtIni,dt,tt,npar,
     &  Gradxnew,Gradynew,Gradznew,Epot,refPE,
     & xnew, ynew, znew)

!MODIFICACION RESORTES: CODIGO NUEVO END
      

      call largest (tt,npar,ncage, xnew,ynew,znew, radConf,radCage,Vol)
      !call eneSpring (tt,npar,ncage,xnew,ynew,znew,x0,y0,z0,sK,Espring)
      call kinEne (npar,ncage, atomass, vxnew,vynew,vznew,
     &     EkinConf,EkinCage, tempConf,tempCage,Vol,dynPress)

* 1 nano-metre= 10 Angst

      Volnano= Vol/1000.
      !Etot= EkinConf + EkinCage + Epot + Espring
      Etot= EkinConf + EkinCage + Epot 

      write(*,11) "  energies [Hartree] t=  ", tt
      write(*,12) "kinetic energy confined  ", EkinConf
      write(*,12) "kinetic energy cage      ", EkinCage
      write(*,12) "potential ab energy      ", Epot
      !write(*,12) "spring energy            ", Espring
      write(*,12) "total energy             ", Etot
      write(*,13) "temperature cage [K]     ", tempCage
      write(*,13) "temperature conf [K]     ", tempConf
      write(*,13) "dynamic pressure [atm]   ", dynPress
      write(*,13) "confinement rad [Ang]    ", radConf
      write(*,13) "cage radius [Ang]        ", radCage
      write(*,13) "confinement vol [nanom^3]", Volnano
      print*

      do i= 1, npar ! <--
                         !
         x(i)= xnew(i)   ! shuffle variables
         y(i)= ynew(i)   !
         z(i)= znew(i)   !
                         !
        vx(i)= vxnew(i)  !
        vy(i)= vynew(i)  !
        vz(i)= vznew(i)  !
                         !
        Gradx(i)= Gradxnew(i)
        Grady(i)= Gradynew(i)
        Gradz(i)= Gradznew(i)
                         !
      end do ! <---------

      cpu1= 0
      if (imethod.eq.1) call estimateCPU (cpu1)

* ------------------------------------------
* evolve positions & speeds with Taylor from t=-dt to t= 0
* ------------------------------------------

      call taylor (npar, ncage, atomass,
     &   x,y,z, vx,vy,vz, Gradx,Grady,Gradz,
     &   xnew,ynew,znew, vxnew,vynew,vznew, dt)

      tt= 0  ! present time

      write (*,*) "--------------------"
      write (*,"(1x,a,i6,1x,f9.2)") "step & time", 0, tt

! MODIFICACION RESORTES: COMENTAREAR
* reuse the wave function
!
!      if (imethod.eq.1) then
!      call inptTera (2)
!      call writeMol (npar,atsymbol,xnew,ynew,znew)
!      end if
!
!      if (imethod.eq.2) then
!      call inptNwchem (2,npar,atsymbol,xnew,ynew,znew)
!      end if
!
!* submit a job, give seconds to close files & get forces
!
!      if (itest.eq.2) then ! <----------
!                                        !
!      if (imethod.eq.1) call system ("./script.tera")
!      if (imethod.eq.2) call system ("./script.nwchem")
!                                        !
!      call sleep (02)                   !
!      call abEne (imethod,dtIni,dt,tt,npar,
!     &  Gradxnew,Gradynew,Gradznew,Epot,refPE)
!                                        !
!      if (imethod.eq.1) call groupCharge (tt,ngrp1,ngrp2,ngrp3)
!                                        !
!      end if ! <------------------------
!
! MODIFICACION RESORTES: COMENTAREAR END

      call largest (tt,npar,ncage, xnew,ynew,znew, radConf,radCage,Vol)
      call eneSpring (tt,npar,ncage,xnew,ynew,znew,x0,y0,z0,sK,Espring)
      call kinEne (npar,ncage, atomass, vxnew,vynew,vznew,
     &     EkinConf,EkinCage, tempConf,tempCage,Vol,dynPress)

      Volnano= Vol/1000.
      !Etot= EkinConf + EkinCage + Epot + Espring
      Etot= EkinConf + EkinCage + Epot 

      write(*,11) "  energies [Hartree] t=  ", tt
      write(*,12) "kinetic energy confined  ", EkinConf
      write(*,12) "kinetic energy cage      ", EkinCage
      write(*,12) "potential ab energy      ", Epot
      !write(*,12) "spring energy            ", Espring
      write(*,12) "total energy             ", Etot
      write(*,13) "temperature cage [K]     ", tempCage
      write(*,13) "temperature conf [K]     ", tempConf
      write(*,13) "dynamic pressure [atm]   ", dynPress
      write(*,13) "confinement rad [Ang]    ", radConf
      write(*,13) "cage radius [Ang]        ", radCage
      write(*,13) "confinement vol [nanom^3]", Volnano
      print*

      do i= 1, npar ! <--
                         !
        x(i)= xnew(i)    ! shuffle variables
        y(i)= ynew(i)    !
        z(i)= znew(i)    !
                         !
        vx(i)= vxnew(i)  !
        vy(i)= vynew(i)  !
        vz(i)= vznew(i)  !
                         !
        Gradx(i)= Gradxnew(i)
        Grady(i)= Gradynew(i)
        Gradz(i)= Gradznew(i)
                         !
      end do ! <---------

* ---------------
* estimate the cpu
* ---------------

      cpu2= 0
      if (imethod.eq.1) call estimateCPU (cpu2)

* cpu time per step

      cputime= (cpu1 + cpu2)/ 2.

* ab-initio computations= nsteps + 2
* sleeping time= nsteps*2 + 30
* radiation + charge uptake + sleeping time computations= 6*cputime+6*2

      cpu= cputime*(nsteps+2.) + nsteps*2+30 + 6*cputime+6*2  ! in sec

      ndays=   int(cpu/ (3600.*24.))
      nhours=  int(cpu/3600. -ndays*24.)
      minutes= int(cpu/60 -ndays*24.*60.-nhours*60.)
      nsec=    int(mod(cpu,60.))

      write(*,14) "number of time steps to evaluate", nsteps
      write(*,20) ndays, nhours, minutes, nsec
      print*

* --------------------------
* initialization of parameters
* --------------------------

  120 sEkin=  0.d0   ! sum of kinetic energies
      sEpot=  0.d0   ! sum of pot energies
      sEtot=  0.d0   ! sum of total energies
      sEtot2= 0.d0   ! sum of squares of kin energies
      stCage= 0.d0   ! sum of cage temperatures
      stConf= 0.d0   ! sum of confined atoms temperatures
      sPres=  0.d0   ! sum of pressures
      sVol=   0.d0   ! sum of volumes

      do k= 1, nsteps
       enekin(k)= 0.d0
       enepot(k)= 0.d0
       tConf(k)=  0.d0
       tCage(k)=  0.d0
       press(k)=  0.d0
      end do

      do i= 1, ncage
       sum(i)= 0.d0  ! summations for every atom speed
      end do

      icount= 0  ! counter of points for the histogram
      kj= 0      ! counter of frames of the simulation
      seed= 15

c     xl=  0.01       ! scan speeds in = [0,xl]
      xl=  0.03       ! an ordinary Maxwell distrib. funct.
      wid= xl/1000.   ! sub-intervals to classify speeds

      nintervals= int(xl/wid) ! total num of sub-intervals
      kdim= 0

      do i= 1, nintervals  ! initialize counters
       kdim= kdim +1
       knt(i)= 0
      end do

      write (*,18) "interval to scan atom speeds [0,xl]", 0.0, xl
      write (*,16) "resolution to classify speeds", wid
      write (*,14) "number of intervals in the histogram", kdim

* open the file for the animation frames

      open (unit=12, file= 'movie.xyz',status= 'unknown')

      tempConf= 0.d0
      tt= trest  ! starting time
      iAft= 10   ! time after breaking bond

* avoid radiation effects when the dynamics just starts

c     ibTime1= 50123
c     ibTime2= 50123
c     ibTime3= 50123

* ----------------------
      print*, "evolving coords and speeds with velVerlet"
* ----------------------

      do k= ibegin+1, nsteps, 1 ! <------
                                         !
      tt= tt + dt                        !
      kj= kj + 1                         !
                                         !
      write (*,*) "--------------------" !
      write (*,"(1x,a,i6,1x,f9.2)") "step & time", k, tt
                                         !
                                         !
* parameters of the distrib funct (depend of the time step)
      call powerSer (kj,dt,xi,c0,c1,c2,c3,aa,bb,cc)
                                         !
                                ! -----------------
                                ! get new positions
                                         !
      call positions (tt,dt,npar,ncage,atomass,
     &  x,y,z, vx,vy,vz, Gradx,Grady,Gradz,
     &  xnew,ynew,znew, Xn,Yn,Zn, x0,y0,z0,sK,
     &  tempEquil,tempConf,c1,c2,cc)     !
                                         !
                                         !
!MODIFICACION RESORTES: NUEVO CODIGO 
                                 ! get new forces
                                         !
!                                         ! wavefunction type
!      call setWFs (k,imethod,ibTime1,ibTime2,ibTime3,iAft,
!     &   npar,atsymbol,xnew,ynew,znew)   !
!                                         !
!      if (itest.eq.2) then               !---
!                                             !
!      if (imethod.eq.1) call system ("./script.tera")
!      if (imethod.eq.2) call system ("./script.nwchem")
!                                             !
!      call sleep (01)                        !
!      call abEne (imethod,dtIni,dt,tt,npar,  !
!     &  Gradxnew,Gradynew,Gradznew,Epot,refPE)
!                                             !
!      if (imethod.eq.1) call groupCharge (tt,ngrp1,ngrp2,ngrp3)
!                                             !
!      end if                             !---
                                         !
      call abEneR (imethod,dtIni,dt,tt,npar,  !
     &  Gradxnew,Gradynew,Gradznew,Epot,refPE,
     &  xnew, ynew, znew)
                                        !
!MODIFICACION RESORTES: NUEVO CODIGO END 

        
                                         ! --------------
                                 ! get new speeds
                                         !
      call speeds (tt,dt,npar,ncage,atomass,
     & x,y,z, vx,vy,vz, Gradx,Grady,Gradz, xnew,ynew,znew,
     & vxnew,vynew,vznew, Gradxnew,Gradynew,Gradznew,
     & Xn,Yn,Zn, Vxn,Vyn,Vzn, x0,y0,z0,sK,
     & tempEquil,tempConf,c0,c1,c2, aa,bb,cc, ftest)
                                         !
                                         !
      do i= 1, ncage                     ! count speeds for histogram
       icount= icount + 1                !
       sum(i)= sum(i) + ftest(i)         ! every cage atom has its
       point= ftest(i)                   ! own summation
       call classify (nintervals,point,knt,xl,wid)
      end do                             !
                                         !
                                         !
                                  ! ------------
                                  ! get energies
                                         !
      call largest (tt,npar,ncage, xnew,ynew,znew, radConf,radCage,Vol)
      !call eneSpring (tt,npar,ncage,xnew,ynew,znew,x0,y0,z0,sK,Espring)
      call kinEne (npar,ncage, atomass, vxnew,vynew,vznew,
     &     EkinConf,EkinCage, tempConf,tempCage,Vol,dynPress)
                                         !
                                         !
      !Etot= EkinConf + EkinCage + Epot + Espring
      Etot= EkinConf + EkinCage + Epot 
                                         !
                                         !
                               ! ------------------
                               ! store dynamic data
                                         !
      enekin(kj)= EkinCage + EkinConf    !
      enepot(kj)= Epot
      !enepot(kj)= Epot + Espring         !
      tCage(kj)= tempCage                !
      tConf(kj)= tempConf                !
      press(kj)= dynPress                !
                                         !
      Volnano= Vol/1000.                 ! vol in nano^3
                                         !
      write(*,11) "  energies [Hartree] t=  ", tt
      write(*,12) "kinetic energy confined  ", EkinConf
      write(*,12) "kinetic energy cage      ", EkinCage
      write(*,12) "potential ab energy      ", Epot
      !write(*,12) "spring energy            ", Espring
      write(*,12) "total energy             ", Etot
      write(*,13) "temperature cage [K]     ", tempCage
      write(*,13) "temperature conf [K]     ", tempConf
      write(*,13) "dynamic pressure [atm]   ", dynPress
      write(*,13) "confinement rad [Ang]    ", radConf
      write(*,13) "cage radius [Ang]        ", radCage
      write(*,13) "confinement vol [nanom^3]", Volnano
                                         !
                                         !
      sEkin=  sEkin + EkinCage +EkinConf ! sum kin energies
      !sEpot=  sEpot + Epot + Espring     ! sum pot energies
      sEpot=  sEpot + Epot     ! sum pot energies
      sEtot=  sEtot + Etot               ! sum tot energies
      sEtot2= sEtot2 + Etot**2           ! sum (tot energies)^2
      stCage= stCage + tempCage          ! sum cage temperatures
      stConf= stConf + tempConf          ! sum confined atoms temps
      sPres=  sPres + dynPress           ! sum pressures
      sVol=   sVol + Vol                 ! sum volumes
                                         !
                                         !
                               ! -----------------
                               ! shuffle variables
      do i= 1, npar                      !
         x(i)= xnew(i)                   !
         y(i)= ynew(i)                   !
         z(i)= znew(i)                   !
                                         !
        vx(i)= vxnew(i)                  !
        vy(i)= vynew(i)                  !
        vz(i)= vznew(i)                  !
                                         !
        Gradx(i)= Gradxnew(i)            !
        Grady(i)= Gradynew(i)            !
        Gradz(i)= Gradznew(i)            !
      end do                             !
                                         !
      call makeMovie (k,npar, atsymbol,x,y,z)
                                         !
                                         !
                                 ! --------------
                                 ! cage shrinking
                                         !
c     call distances (tt,npar, x,y,z)    !
                                         !
      if (mod(k,ifreqShr).eq.0) then     !
      print*, "shrinking the cage"       !
      call shrink (npar, ncage, x0,y0,z0, x,y,z, radiusf)
      end if                             !
                                         !
!MODIFICACION CODIGO: COMENTAREAR        ! ------------------
!                               ! implicit radiation
!                                         !
!                                         !
!      if ((k.eq.ibTime1).or.(k.eq.ibTime2).or.(k.eq.ibTime3)) then
!      write(*,"(/,a)") "ionizing radiation: bond stochastically broken"
!                                         !
!      call breakBond (npar,ncage,atsymbol,atomass,x,y,z)
!      call inptTera (4)                  ! write new input file
!      call writeMol (npar,atsymbol,x,y,z)! and coords
!                                         !
!      ene1= Epot                         ! pot ene from previous step
!                                         !
!                                         ! get new forces
!      if (itest.eq.2) then               !---
!                                             !
!      if (imethod.eq.1) call system ("./script.tera")
!      if (imethod.eq.2) call system ("./script.nwchem")
!                                             !
!      call sleep (01)                        !
!      call abEne (imethod,dtIni,dt,tt,npar,  !
!     &  Gradx,Grady,Gradz,Epot,refPE)        !
!                                             !
!      if (imethod.eq.1) call groupCharge (tt,ngrp1,ngrp2,ngrp3)
!                                             !
!      end if                             !---
!                                         !
!      ene2= Epot                         ! in hartree
!      radEne= (ene2-ene1)* 27.2107       ! in eV
!      wLength= (455.619348)/(ene2-ene1)  ! in Angst
!                                         !
!      write(*,"(/,a,f9.5)") "ionization energy [eV]", radEne
!      write(*,"(a,e14.5)")  "energy wavelength [Angst]", wLength
!      call makeMovie (k,npar, atsymbol,x,y,z)
!      end if                             !
!                                         !
!                                 ! -------------
!                                 ! charge uptake
!                                         !
!                                         !
!      if ((k.eq.ibTime1+iAft).or.(k.eq.ibTime2+iAft).or.
!     &    (k.eq.ibTime3+iAft)) then      !
!      write(*,"(/,a)") "charge uptake from environment"
!                                         !
!      call inptTera (1)                  ! write initial input file
!      call writeMol (npar,atsymbol,x,y,z)! and coords
!                                         !
!      if (itest.eq.2) then               !---
!                                             !
!      if (imethod.eq.1) call system ("./script.tera")
!      if (imethod.eq.2) call system ("./script.nwchem")
!                                             !
!      call sleep (01)                        !
!      call abEne (imethod,dtIni,dt,tt,npar,  !
!     &  Gradx,Grady,Gradz,Epot,refPE)        !
!                                             !
!      if (imethod.eq.1) call groupCharge (tt,ngrp1,ngrp2,ngrp3)
!                                             !
!      end if                             !---
!      end if                             !
!                                         !
!MODIFICACION CODIGO: COMENTAREAR END

                              ! -------------------
                              ! store relevant data
                                         !
      if (mod(k,ifreqStor).eq.0) then    !
      write(*,"(/,a)") "storing data for restart"
                                         !
      call restartOlla (1,k,dt,tt,npar,ncage,x,y,z,vx,vy,vz,
     &   Gradx,Grady,Gradz, x0,y0,z0,sK,refPE,ibTime1,ibTime2,ibTime3)
                                         !
      print*                             ! average temps
      call averageT (k,tCage,ifreqStor,averTemp)
      write(*,"(1x,a,f9.2)") "average cage temp [K]", averTemp
                                         !
      call averageT (k,tConf,ifreqStor,averTemp)
      write(*,"(1x,a,f9.2)") "average confined temp [K]", averTemp
      end if                             !
                                         !
      end do ! <-------------------------

      close (12) ! close movie file

* -------------------------
* averages of main energies
* -------------------------

      Ekin=  sEkin  /dble(kj)    ! kin energy
      Epot=  sEpot  /dble(kj)    ! pot energy
      Etot=  sEtot  /dble(kj)    ! tot energy
      Etot2= sEtot2 /dble(kj)    ! (tot energy)^2
      TempCage= stCage /dble(kj) ! cage temp
      TempConf= stConf /dble(kj) ! confined atoms temp
      Pres=  sPres  /dble(kj)    ! pressure
      Vol=   sVol   /dble(kj)    ! volume

      Volnano= Vol/1000.

* energy fluctuation

      Efluct= dsqrt(dabs(Etot2 - Etot**2))/ Etot

* average speed of cage particles per frame and per atom
* (~Maxwell average speed for the given temperature)

      sum2= 0.d0

      do i= 1, ncage
       sum2= sum2 + sum(i)
      end do

      avSpd= sum2/ (dble(kj)*dble(ncage))

* info for the Vprom speed:
* http://hyperphysics.phy-astr.gsu.edu/hbase/kinetic/kintem.html

* [sqrt(RT/m)]= sqrt [J/(mol*Kelvin) * Kelvin/amu]= fac *Angst/fs

      fac= 1.d-5* sqrt(1000.)

      Vprom= fac* sqrt(2.*8.3144621 *tempEquil/atomass(1)) ! in Angst/fs

* create the speeds histogram and the analytic Maxwell speed distrib. funct.

      call printGauss (nintervals,icount,knt,xl,wid)
      call analytic (atomass,tempEquil,xl,wid)

      write(*,*)  "--------------------"
      write(*,14) "number of intervals in the histogram", kdim
      write(*,14) "number of classified speeds", icount
      write(*,16) "average rnd speed of cage particles [Angst/fs]",avSpd
      write(*,16) "maxwellian speed [Angst/fs]", Vprom
      write(*,11) "maxwellian speed [m/s]", Vprom *1.d+5

      print*
      print*, "the speed histogram of the cage particles was created"
      print*, "the analytic Maxwell speed distrib. funct. was created"

* ------------------------------
* printing of physical variables
* ------------------------------

* print energies in terms of t

      print*, "--------------------------------------------"
      print*, "time   Ekin [Hartree] Epot [Hartree] Tcage [K]  Tconf [K]
     &  Pdyn [GPa]"

      tt= 0.d0
      sum1= 0.
      sum2= 0.

      do kk= 1, nsteps, 1 ! <---
                                !
       tt= tt + dt              !
       sum1= sum1 + tCage(kk)   !
       sum2= sum2 + tCage(kk)**2!
                                !
       write(*,10) tt, enekin(kk),enepot(kk),
     &      tCage(kk),tConf(kk),press(kk)
                                !
      end do ! <----------------

      tempfluct= dsqrt(dabs(sum2 - sum1**2))/ sum1

      print*
      write(*,*) "number of computed time steps= ", kj
      write(*,*)

      write(*,*) "  average values [Hartree]"
      write(*,12) "kinetic ene             ", Ekin
      write(*,12) "potential ene           ", Epot
      write(*,12) "total ene               ", Etot
      write(*,12) "(total ene)^2           ", Etot2
      write(*,12) "ene fluctuat            ", Efluct
      write(*,13) "cage temperature [K]    ", TempCage
      write(*,13) "confined atoms temp [K] ", TempConf
      write(*,13) "cage temp fluctuat      ", tempfluct
      write(*,13) "dynamic pressure [atm]  ", Pres
      write(*,13) "volume of conf [nanom^3]", Volnano

      print*, "--------------------------------------------"
      write(*,*) "final parameters of the dynamics"
      write(*,16) "time step [fs] ", dt
      write(*,16) "total simulation time [fs] ", Time
      write(*,16) "rfac      ", rfac
      write(*,16) "vfac      ", vfac
      write(*,14) "npar      ", npar
      write(*,14) "ncage     ", ncage
      write(*,22) "charge grps", ngrp1,ngrp2,ngrp3
      write(*,16) "tempEquil ", tempEquil
      write(*,16) "spring K  ", sK(1)  ! all are the same
      write(*,16) "friction  ", xi
      write(*,16) "radiusf   ", radiusf
      write(*,14) "ifreqShr  ", ifreqShr
      write(*,14) "ifreqStor ", ifreqStor
      write(*,14) "irestDyn  ", irestDyn
      write(*,14) "iwavFunct ", iwavFunct
      write(*,*)

      print*, "created files: movie.xyz, analytic, rndgauss"

c     call units
  100 call citation (imethod)
      call figures
c     call system ("rm mol.xyz a.out script.tera terachem.inp")

   10 format (1x,f8.2,2(1x,d12.5), 3(1x,f10.3))
   11 format (1x,a,2x,f10.3)
   12 format (1x,a,2x,d16.8)
   13 format (1x,a,2x,f14.4)
   14 format (1x,a,2x,i10)
   16 format (1x,a,2x,f16.8)
   18 format (1x,a,2f9.5)
   20 format (1x,"shortest estimated cpu time for the dynamics",/,
     &     i3," days:", i3," hrs:", i3," min:", i3," sec")
   22 format (1x,a,3i6)

      !close(32)

      stop
      end
* *****************************
* create an analytic gaussian function

* info on the normal distribution function:
* http://hyperphysics.phy-astr.gsu.edu/hbase/kinetic/maxspe.html#c3
* http://en.wikipedia.org/wiki/Normal_distribution
* variance= sigma^2
* standard deviation= sigma
* *****************************
      subroutine analytic (atomass,tempEquil,xl,wid)
      implicit double precision (a-h,o-z)
      include 'parametros'

      dimension atomass(nmax)
      real*8 kB

      open (16,file="analytic",status="unknown")

      pi= 3.1415926535898

* Boltzmann constant

      kB= 1.3806504d-23    ! in Joule/Kelvin

* facKT *(Angst/fs)^2= (Joule/Kelvin)*(Kelvin/amu)= J/amu

      facKT= 1.0d17/ 1.6605402

* normalization factor

      v1= atomass(1)/ (2.*kB*tempEquil) *1./facKT ! in (fs/Angst)^2
      fac= 4.*pi*(v1/pi)**(3./2.) ! in (fs/Angst)^3

* define intervals

      nn= 125
      dy= 0.0002
      y= -dy
      ik= 0

* loop over y-values

      do ki= 1, nn ! <---
                         !
      ik= ik + 1         !
      y= y + dy          !
                         ! in fs/Angst
      distrib= fac*y**2 *exp(-v1*y**2)
                         !
      write (16,"(2f9.4)") y, distrib
                         !
      end do ! <---------

      close (16)

      return
      end
* *****************************
* classify speeds in sub-intervals to build a histogram
* *****************************
      subroutine classify (nintervals,point,knt,xl,wid)
      implicit double precision (a-h,o-z)
      include 'parametros'

      dimension knt(m1)  ! m1= ncage*Time~ 180*1000= 180 mil

      j= 0
      sum= 0.d0  ! values of point >= 0

      do ki= 1, nintervals ! <--
                                ! creat the j^th
      j=j+1                     ! sub-interval with
                                ! limits [sum,v1]
      sum= sum + wid            !
      v1=  sum + wid            !
                                !
c     write(*,"(i5,2f10.5)") j, sum, v1
                                !
                                ! is the point in the
                                ! interval [sum,v1)?
                                !
      if ((point.ge.sum).and.(point.lt.v1)) knt(j)= knt(j)+1
                                !
      end do ! <----------------

      return
      end
* *******************************
* evolve positions using vVerlet integrators
* *******************************
      subroutine positions (time,dt,npar,ncage,atomass,
     &  x,y,z, vx,vy,vz, Gradx,Grady,Gradz,
     &  xnew,ynew,znew, Xn,Yn,Zn, x0,y0,z0,sK,
     &  tempEquil,tempConf,c1,c2,cc)

      implicit double precision (a-h,o-z)
      include 'parametros'

      external gauss

      dimension atomass(nmax)

      dimension  x(nmax),  y(nmax),  z(nmax)
      dimension vx(nmax), vy(nmax), vz(nmax)
      dimension Gradx(nmax), Grady(nmax), Gradz(nmax)

      dimension xnew(nmax), ynew(nmax), znew(nmax)
      dimension Xn(nmax), Yn(nmax), Zn(nmax)

      dimension x0(nmax), y0(nmax), z0(nmax), sK(nmax)

      integer ilabel(nmax)

      real*8 kB

* Boltzmann constant

      kB= 1.3806504d-23    ! in Joule/Kelvin= Kg (m/s)^2 /Kelvin

* faccl * Angst= (1/amu)*(Hartree/Bohr)*(fs^2)

      faccl= 0.49616184    ! factors to change units

* fpring *(amu/fs^2) = Newton/cm

      fspring= 1.d0/ 16.605402

* ------------------
* evolve positions of cage atoms
* ------------------

      do i= 1, ncage ! <------------
                                    !
      amas= atomass(i)              !
                          ! -----------------
                          ! systematic forces
!MODIFICACION RESORTES: NUEVO CODIGO!
!      a1x= -faccl *Gradx(i) /amas   !
!      a1y= -faccl *Grady(i) /amas   ! ab-initio forces
!      a1z= -faccl *Gradz(i) /amas   !
!                                    !
      a1x= fspring *Gradx(i) /amas   !
      a1y= fspring *Grady(i) /amas   ! ab-initio forces
      a1z= fspring *Gradz(i) /amas   !

!MODIFICACION RESORTES: NUEVO CODIGO! END


                                    ! spring forces F= -k*dx
!      a2x= -fspring *sK(i) *(x(i)-x0(i)) /amas
!      a2y= -fspring *sK(i) *(y(i)-y0(i)) /amas ! in Angst/fs^2
!      a2z= -fspring *sK(i) *(z(i)-z0(i)) /amas
!      
      a2x= 0.0
      a2y= 0.0 
      a2z= 0.0
                                    !
      ax= a1x + a2x                 ! total
      ay= a1y + a2y                 !
      az= a1z + a2z                 !
                                    !
                          ! -----------------
                          ! stochastic forces
                                    !
      pkT= kB*tempEquil/amas *cc    ! in Angst^2
                                    !
      etax= gauss()                 ! random numbers
      etay= gauss()                 !
      etaz= gauss()                 ! stochastic positions
                                    ! Eqs 1.13.15 of chapter Langevin
      Xn(i)= etax * sqrt(pkT)       !
      Yn(i)= etay * sqrt(pkT)       ! in Angst
      Zn(i)= etaz * sqrt(pkT)       !
                                    !
                                    !
                          ! ------------------
                          ! position variables
                                    !
* x(t)= x(t-dt) + c1*v(t-dt)*dt + c2*F(t-dt)/m *dt^2 + Xn(dt)
                                    !
                                    ! in Angst
      xnew(i)= x(i) + c1*vx(i)*dt + c2*ax * dt**2 + Xn(i)
      ynew(i)= y(i) + c1*vy(i)*dt + c2*ay * dt**2 + Yn(i)
      znew(i)= z(i) + c1*vz(i)*dt + c2*az * dt**2 + Zn(i)
                                    !
                                    !
                          ! -------------------
                          ! debugging positions
                                    !
c     if ((time.le.1.0).and.(i.le.5)) then
c     if (i.le.3) then              !
c                                   !
c     write(*,15) i, sqrt(pkT)      !
c     write(*,51) "c1,c2, dt", c1,c2,dt
c     write(*,51) "aki posicion", x(i),c1*vx(i)*dt,c2*ax*dt**2,Xn(i)
c     write(*,51) "aki posicion", vx(i),vy(i),vz(i)
c                                   !
c  51 format (a,4f10.4)             !
c     end if                        !
                                    !
      end do   ! <------------------

* ------------------------
* get atoms close to the walls
* ------------------------

* threshold dist to the cage

      thrDist= 3.50
      icont= 0

c     do i= ncage+1, npar ! <-------  atoms in the cage
c     ilabel(i)= 0                  !
c     do j= 1, ncage                ! atoms of the cage
c                                   !
c     dx= (x(i)-x(j))**2            !
c     dy= (y(i)-y(j))**2            !
c     dz= (z(i)-z(j))**2            !
c     dist= sqrt(dx+dy+dz)          ! find atoms close
c                                   ! to the cage
c     if (dist.le.thrDist) then     !
c     icont= icont +1               !
c     ilabel(i)= 1                  ! tag this atom
c     goto 20                       !
c     end if                        !
c                                   !
c     end do                        !
c  20 continue                      !
c     end do ! <--------------------

c     write(*,"(a,1x,f6.3)") "distance from cage for cooling down",
c    &               thrDist
c     write(*,"(a,1x,i4)") "num of atoms close to walls", icont

* ------------------
* evolve positions of confined atoms
* ------------------

* threshold temp

      thrTemp= tempEquil + 1200.

      do i= ncage+1, npar ! <-------
                                    !
      amas= atomass(i)              !
                                    !
                           ! -----------------
                           ! systematic forces
                                    !
!MODIFICACION RESORTES: NUEVO CODIGO
!                              
!      ax= -faccl *Gradx(i) /amas    ! ab-initio forces
!      ay= -faccl *Grady(i) /amas    !
!      az= -faccl *Gradz(i) /amas    !
! 
      a1x= fspring *Gradx(i) /amas   !
      a1y= fspring *Grady(i) /amas   ! ab-initio forces
      a1z= fspring *Gradz(i) /amas   !

                                   !--
!MODIFICACION RESORTES: NUEVO CODIGO
                                       !
c     if ((ilabel(i).eq.1).and.(tempConf.gt.thrTemp)) then
c     if (ilabel(i).eq.1) then         ! cool down this atom
c                                      !
c                                      !
c                             ! -----------------
c                             ! stochastic forces
c                                      !
c     pkT= kB*tempEquil/amas *cc       ! in Angst^2
c                                      !
c     etax= gauss()                    ! random numbers
c     etay= gauss()                    !
c     etaz= gauss()                    ! stochastic positions
c                                      ! Eqs 1.13.15 of chapter Langevin
c     Xn(i)= etax * sqrt(pkT)          !
c     Yn(i)= etay * sqrt(pkT)          ! in Angst
c     Zn(i)= etaz * sqrt(pkT)          !
c                                      !
* x(t)= x(t-dt) + c1*v(t-dt)*dt + c2*F(t-dt)/m *dt^2 + Xn(dt)
c                                      !
c                                      ! in Angst
c     xnew(i)= x(i) + c1*vx(i)*dt + c2*ax * dt**2 + Xn(i)
c     ynew(i)= y(i) + c1*vy(i)*dt + c2*ay * dt**2 + Yn(i)
c     znew(i)= z(i) + c1*vz(i)*dt + c2*az * dt**2 + Zn(i)
c     end if                           !
                                    !--
                                    !
c     if (ilabel(i).eq.0) then      ! normal atom
* x(t)= x(t-dt) + v(t-dt)*dt + F(t-dt)/m *dt^2/2
                                    ! in Angst
      xnew(i)= x(i) + vx(i)*dt + ax/2.d0 *dt**2
      ynew(i)= y(i) + vy(i)*dt + ay/2.d0 *dt**2
      znew(i)= z(i) + vz(i)*dt + az/2.d0 *dt**2
c     end if                        !
                                    !
      end do   ! <------------------

      return
      end
* ************************************
* power series of coefficients
* this module is called once
* ************************************
      subroutine powerSer (kj,dt,xi,c0,c1,c2,c3,aa,bb,cc)
      implicit double precision (a-h,o-z)
      include 'parametros'

* xi in 1/ps,  1 ps= 1000 fs

      hxi= xi/1000.0   ! in 1/fs

      z= abs(hxi)*dt   ! dimensionless

* expansions from Paterlini's article (refer to the Appendix,
* some of the numbers there are mistaken)
* Chem. Phys. vol. 236, 243-252, 1998

* refer to expressions \label{eqmopafl.9} of my
* Langevin chapter for the parameters, c0, c1, etc

* -----------------------
* c0= e^{-z}  ;  c1= (1-c0)/z  ;  c2= (1-c1)/z  ;  c3= (1/2-c2)/z
* -----------------------

* the coefficients are dimensionless

      c0= 1.    -z     +z**2/2.  -z**3/6.  +z**4/24. -z**5/120.
      c1= 1.    -z/2.  +z**2/6.  -z**3/24. +z**4/120.
      c2= 1./2. -z/6.  +z**2/24. -z**3/120.
      c3= 1./6. -z/24. +z**2/120.

      if (kj.eq.1) then ! <------
                                 !
      print*                     !
      write(*,*) "dynamics & distrib function factors"
      write(*,*) "(for internal use)"
      write(*,13) "c0=", c0      !
      write(*,13) "c1=", c1      !
      write(*,13) "c2=", c2      !
      write(*,13) "c3=", c3      !
      print*                     !
                                 !
      end if ! <-----------------

* facKT *(Angst/fs)^2= (Joule/Kelvin)*(Kelvin/amu)

      facKT= 1.0d17/1.6605402

* -----------------------
* refer to expressions \ref{intteqom.3} of my
* Langevin chapter for the parameters aa, bb, cc
* -----------------------

* we include the factor (kB*Tref/mass)
* after the call to this module

      aa= 1. -z +2.*z**2/3. -z**3/3. +2.*z**4/15. -2.*z**5/45.
c     aa= 2.*kB*T/mass *z *aa *facKT   <------------ in (Angst/fs)^2
      aa= 2.*z*aa *facKT

      bb= 1.-z + 7.*z**2/12. -z**3/4. +31.*z**4/360. -z**5/40.
c     bb= kB*T/mass *z*dt *bb *facKT  <-------- in Angst^2/fs
      bb= z*dt*bb *facKT

      cc= 1.-3.*z/4.+7.*z**2/20.-z**3/8.+31.*z**4/840.-3.*z**5/320.
c     cc= 2.*kB*T/(3.*mass) *z*dt**2 *cc *facKT  <-------- in Angst^2
      cc= 2./3. *z*dt**2 *cc *facKT

* print without the factor facKT

      if (kj.eq.1) then ! <------
                                 !
      write(*,13) "a=", aa/facKT !
      write(*,13) "b=", bb/facKT !
      write(*,13) "c=", cc/facKT !
      print*                     !
                                 !
      write(*,14) "ac=", aa*cc/facKT**2
      write(*,14) "ac-b^2=", (aa*cc-bb**2)/facKT**2
      write(*,14) "sqrt((ac-b^2)/c)=", sqrt((aa*cc-bb**2)/(cc*facKT))
      write(*,14) "b/c=", bb/cc  !
                                 !
      end if ! <-----------------

   13 format (1x,a,1x,f10.4)
   14 format (1x,a,d13.4)

      return
      end
* *****************************
* print the gaussian curve for the random speeds
* *****************************
      subroutine printGauss (nintervals,npoints,knt,xl,wid)
      implicit double precision (a-h,o-z)
      include 'parametros'

      dimension knt(m1)

      open (16,file="rndgauss", status="unknown")

* numeric integration of the function: base * height
* vnorm= wid*n1 + wid*n2 + wid*n3= dx*(n1+n2+...)= wid*npoints

      vnorm= wid*npoints  ! normalization factor
      j= 0
      sum= 0.d0

* print number of points at each sub-interval

      do ki= 1, nintervals ! <--
                                !
      j= j+1                    ! j^th interval
                                !
      sum= sum + wid            !
      v1=  sum + wid            ! get central point
      vc= (sum + v1)/ 2.        ! of the interval
                                !
      write (16,"(2f12.3)") vc, dble(knt(j))/vnorm
                                !
      end do ! <----------------

      close (16)

      return
      end
* *******************************
* evolve speeds using vVerlet integrators
* *******************************
      subroutine speeds (time,dt,npar,ncage,atomass,
     & x,y,z, vx,vy,vz, Gradx,Grady,Gradz,
     & xnew,ynew,znew, vxnew,vynew,vznew, Gradxnew,Gradynew,Gradznew,
     & Xn,Yn,Zn, Vxn,Vyn,Vzn, x0,y0,z0,sK,
     & tempEquil,tempConf,c0,c1,c2, aa,bb,cc, ftest)

      implicit double precision (a-h,o-z)
      include 'parametros'

      external gauss

      dimension atomass(nmax), ftest(nmax)

      dimension  x(nmax),  y(nmax),  z(nmax)
      dimension vx(nmax), vy(nmax), vz(nmax)
      dimension Gradx(nmax), Grady(nmax), Gradz(nmax)

      dimension  xnew(nmax),  ynew(nmax),  znew(nmax)
      dimension vxnew(nmax), vynew(nmax), vznew(nmax)
      dimension Gradxnew(nmax), Gradynew(nmax), Gradznew(nmax)

      dimension  Xn(nmax),  Yn(nmax),  Zn(nmax)
      dimension Vxn(nmax), Vyn(nmax), Vzn(nmax)

      dimension x0(nmax), y0(nmax), z0(nmax), sK(nmax)

      integer ilabel(nmax)

      real*8 kB

* Boltzmann constant

      kB= 1.3806504d-23    ! in Joule/Kelvin= Kg (m/s)^2 /Kelvin

* fvel *(Angst/fs)= (1/amu) *(Hartree/Bohr) *fs

      fvel= 0.49616184     ! factors to change units

* fspring *(amu/fs^2) = Newton/cm

      fspring= 1.d0/ 16.605402

* facKT *(Angst/fs)^2= (Joule/Kelvin)*(Kelvin/amu)

      facKT= 1.0d17/1.6605402

* -----------------------
* evolve speeds of cage atoms
* -----------------------

      do i= 1, ncage ! <------------
                                    !
      amas= atomass(i)              !
                                    !
                        ! ---------------------
                        ! old systematic forces
                                    !
!MODIFICACION RESORTES: NUEVO CODIGO
!      a1ax= -fvel *Gradx(i) /amas   ! ab-initio
!      a1ay= -fvel *Grady(i) /amas   ! in Angst/fs^2
!      a1az= -fvel *Gradz(i) /amas   !
!                                    !
      a1ax= -fspring *Gradx(i) /amas   ! ab-initio
      a1ay= -fspring *Grady(i) /amas   ! in Angst/fs^2
      a1az= -fspring *Gradz(i) /amas   !
                                    ! F= -k*dX
!MODIFICACION RESORTES: NUEVO CODIGO

!      a1bx= -fspring* sK(i)* (x(i)-x0(i)) /amas
!      a1by= -fspring* sK(i)* (y(i)-y0(i)) /amas ! in Angst/fs^2
!      a1bz= -fspring* sK(i)* (z(i)-z0(i)) /amas
! 
      a1bx= 0.0
      a1by= 0.0
      a1bz= 0.0
!
      a1x= a1ax + a1bx              ! total
      a1y= a1ay + a1by              !
      a1z= a1az + a1bz              !
                                    !
                        ! ---------------------
                        ! new systematic forces
                                    !
!MODIFICACION RESORTES: NUEVO CODIGO

!      a2ax= -fvel *Gradxnew(i) /amas! ab-initio
!      a2ay= -fvel *Gradynew(i) /amas! in Angst/fs^2
!      a2az= -fvel *Gradznew(i) /amas!
! 
      a2ax= fspring *Gradxnew(i) /amas! ab-initio
      a2ay= fspring *Gradynew(i) /amas! in Angst/fs^2
      a2az= fspring *Gradznew(i) /amas!
                                    !
!MODIFICACION RESORTES: NUEVO CODIGO
                                    ! F= -k*dX
!      a2bx= -fspring* sK(i)* (xnew(i)-x0(i)) /amas
!      a2by= -fspring* sK(i)* (ynew(i)-y0(i)) /amas ! in Angst/fs^2
!      a2bz= -fspring* sK(i)* (znew(i)-z0(i)) /amas
!                                    !
      a2bx= 0.0
      a2by= 0.0
      a2bz= 0.0
      
!
      a2x= a2ax + a2bx              ! total
      a2y= a2ay + a2by              !
      a2z= a2az + a2bz              !
                                    !
                          ! -----------------
                          ! stochastic speeds
                                    !
      pKT= kB*tempEquil/amas        !
      pa= pKT *aa                   ! (Angst/fs)^2
      pb= pKT *bb                   !  Angst^2/fs
      pc= pKT *cc                   !  Angst^2
                                    !
      etax= gauss()                 ! random numbers
      etay= gauss()                 ! Eqs 1.13.15 of chapter Langevin
      etaz= gauss()                 !
                                    ! in Angst/fs
      Vxn(i)= etax *sqrt(pa-pb**2/pc) +pb/pc *Xn(i)
      Vyn(i)= etay *sqrt(pa-pb**2/pc) +pb/pc *Yn(i)
      Vzn(i)= etaz *sqrt(pa-pb**2/pc) +pb/pc *Zn(i)
                                    !
                                    !
                                    !
                 ! -------------------------------------
                 ! reproduce a Maxwellian distrib funct.
                                    !
* http://hyperphysics.phy-astr.gsu.edu/hbase/kinetic/maxspe.html#c3
* run in test form, use the interval below for speeds and plot the
* analytic and rndgauss data        !
                                    !
*     xl=  0.03     ! interval to scan speeds = [-xl,xl]
*     wid= xl/1000. ! resolution of sub-intervals to scan speeds
                                    !
c     Vxn(i)= etax *sqrt(pKT*facKT) ! to reproduce
c     Vyn(i)= etay *sqrt(pKT*facKT) ! the Maxwell
c     Vzn(i)= etaz *sqrt(pKT*facKT) ! distrib. funct.
                                    !
                                    !
                                    ! test function just for rnd speeds
      ftest(i)= sqrt (Vxn(i)**2 + Vyn(i)**2 + Vzn(i)**2)
                                    !
                                    !
                            ! ---------------
                            ! speed variables
                                    !
* v(t)= c0*v(t-dt) + (c1-c2)*a(t-dt)*dt + c2*a(t)*dt + Vn(dt)
                                    !
                                    ! in Angst/fs
      vxnew(i)= c0*vx(i) + (c1-c2)*a1x*dt + c2*a2x*dt +Vxn(i)
      vynew(i)= c0*vy(i) + (c1-c2)*a1y*dt + c2*a2y*dt +Vyn(i)
      vznew(i)= c0*vz(i) + (c1-c2)*a1z*dt + c2*a2z*dt +Vzn(i)
                                    !
                                    !
                           ! ----------------
                           ! debugging speeds
                                    !
c     if ((time.lt.9.65).and.(i.eq.1)) then
c     if (i.le.3) then              !
c                                   !
c     write(*,51) "aki vx   ",vx(i),vy(i),vz(i)
c     write(*,51) "aki a1ax ",a1ax,a1ay,a1az
c     write(*,51) "aki a1bx ",a1bx,a1by,a1bz
c     write(*,51) "aki a2ax ",a2ax,a2ay,a2az
c     write(*,51) "aki a2bx ",a2bx,a2by,a2bz
c     write(*,51) "aki Vxn  ",Vxn(i),Vyn(i),Vzn(i)
c     write(*,51) "aki vxnew",vxnew(i),vynew(i),vznew(i)
c     write(*,51) "aki ", sqrt(pa-pb**2/pc) ,pb* Xn(i)/ pc
c     print*                        !
c                                   !
c  51 format (a,3f10.4)             !
c     end if                        !
                                    !
      end do   ! <------------------

* ------------------------
* get atoms close to the walls
* ------------------------

* threshold dist to the cage

      thrDist= 3.50
      icont= 0

c     do i= ncage+1, npar ! <-------  atoms in the cage
c     ilabel(i)= 0                  !
c     do j= 1, ncage                ! atoms of the cage
c                                   !
c     dx= (x(i)-x(j))**2            !
c     dy= (y(i)-y(j))**2            !
c     dz= (z(i)-z(j))**2            !
c     dist= sqrt(dx+dy+dz)          ! find atoms close
c                                   ! to the cage
c     if (dist.le.thrDist) then     !
c     icont= icont +1               !
c     ilabel(i)= 1                  ! tag this atom
c     goto 20                       !
c     end if                        !
c                                   !
c     end do                        !
c  20 continue                      !
c     end do ! <--------------------

c     write(*,"(a,1x,f6.3)") "distance from cage for cooling down",
c    &               thrDist
c     write(*,"(a,1x,i4)") "num of atoms close to walls", icont

* ------------------
* evolve speeds of confined atoms
* ------------------

* threshold temperature

      thrTemp= tempEquil + 1200.

      do i= ncage+1, npar ! <----
                                 !
      amas= atomass(i)           !
                        ! ----------------
                        ! systematic forces

!MODIFICACION RESORTE: NUEVO CODIGO      !

!      a1x= -fvel *Gradx(i) /amas ! ab-initio
!      a1y= -fvel *Grady(i) /amas !
!      a1z= -fvel *Gradz(i) /amas !
!                                 !
!      a2x= -fvel *Gradxnew(i) /amas ! ab-initio
!      a2y= -fvel *Gradynew(i) /amas
!      a2z= -fvel *Gradznew(i) /amas
                                 !
      a1x= fspring *Gradx(i) /amas ! ab-initio
      a1y= fspring *Grady(i) /amas !
      a1z= fspring *Gradz(i) /amas !
                                 !
      a2x= fspring *Gradxnew(i) /amas ! ab-initio
      a2y= fspring *Gradynew(i) /amas
      a2z= fspring *Gradznew(i) /amas
                                 !---
!MODIFICACION RESORTE: NUEVO CODIGO END
                                     !
c     if ((ilabel(i).eq.1).and.(tempConf.gt.thrTemp)) then
c     if (ilabel(i).eq.1) then       ! cooldown this atom
c                                    !
c                          ! -----------------
c                          ! stochastic speeds
c                                    !
c     pKT= kB*tempEquil/amas         !
c     pa= pKT *aa                    ! (Angst/fs)^2
c     pb= pKT *bb                    !  Angst^2/fs
c     pc= pKT *cc                    !  Angst^2
c                                    !
c     etax= gauss()                  ! random numbers
c     etay= gauss()                  ! Eqs 1.13.15 of chapter Langevin
c     etaz= gauss()                  !
c                                    ! in Angst/fs
c     Vxn(i)= etax *sqrt(pa-pb**2/pc)+ pb/pc *Xn(i)
c     Vyn(i)= etay *sqrt(pa-pb**2/pc)+ pb/pc *Yn(i)
c     Vzn(i)= etaz *sqrt(pa-pb**2/pc)+ pb/pc *Zn(i)
c                                    !
* v(t)= c0*v(t-dt) + (c1-c2)*a(t-dt)*dt + c2*a(t)*dt + Vn(dt)
c                                    !
c                                    ! in Angst/fs
c     vxnew(i)= c0*vx(i) + (c1-c2)*a1x*dt + c2*a2x*dt +Vxn(i)
c     vynew(i)= c0*vy(i) + (c1-c2)*a1y*dt + c2*a2y*dt +Vyn(i)
c     vznew(i)= c0*vz(i) + (c1-c2)*a1z*dt + c2*a2z*dt +Vzn(i)
c     end if                         !
                                     !
                                 !---
                                 !
c     if (ilabel(i).eq.0) then   ! normal atom
* v(t)= v(t-dt) + [F(t-dt)/m + F(t)/m]/2 * dt
                                 ! in Angst/fs
      vxnew(i)= vx(i)+ (a1x + a2x)/ 2.d0*dt
      vynew(i)= vy(i)+ (a1y + a2y)/ 2.d0*dt
      vznew(i)= vz(i)+ (a1z + a2z)/ 2.d0*dt
c     end if                     !
                                 !
      end do   ! <---------------

      return
      end
* *******************************
* make the distribution of rnd numbers
* to satisfy a gaussian function with
* variance 1 and zero mean
* *******************************
      function gauss()
      implicit double precision (a-h,o-z)

      real*8 zbqlu01

* zbqlu01 is in [0,1]

   30 v1= 2.*zbqlu01() -1.
      v2= 2.*zbqlu01() -1.

      r= v1*v1 + v2*v2

      if (r.ge.1) goto 30

      fac= sqrt (-2. *log(r)/r)
      gauss= v2*fac

      return
      end
* *****************************
* ***********************************
* get ab-initio energy & gradients
      subroutine abEneR (imethod,dtIni,dt,time,npar,
     &      Grdx,Grdy,Grdz,Epot,refPE, X, Y, Z)

      implicit double precision (a-h,o-z)
      include 'parametros'

      character qline*(maxchr)
      dimension itok(maxtok), jtok(maxtok)

      dimension xcrds(10)
      dimension Grdx(nmax), Grdy(nmax), Grdz(nmax)
      dimension X(nmax), Y(nmax), Z(nmax), Id(nmax)
      dimension GX(nmax), GY(nmax), GZ(nmax), zK(nmax)
      dimension VTk(nmax)


      external kparse

* initialize variables

      xenergy= 0.d0
      Epot=  0.d0

      do i= 1, 10
      xcrds(i)= 0.d0
      end do

      do k= 1, npar
       Grdx(k)= 0.d0
       Grdy(k)= 0.d0
       Grdz(k)= 0.d0
      end do


      
* -----------------------------------
* ---------------------------------
* obtain the largest forces in x,y,z
* ---------------------------------

      big=  -100.d0
      xbig= -100.d0
      ybig= -100.d0
      zbig= -100.d0

      !Vector de Relaciones entre de la particula k
      !con las 3 que interactua


      open (unit=31, file = 'settingGeometry.dat', status ='unknown')
      
      open (unit=33, file = 'settingConstantes.dat', status ='unknown')
     

      do k = 1, npar

        j = 2*(k-1)
        read(31,*) Id(j+1), Id(j+2)

      end do

 
      do k = 1, npar

        j = 2*(k-1)
        read(33,*) zK(j+1), zK(j+2)
      end do

      close (33)
      close (31)

      do k = 1, npar
        
        GX(k) = 0.d0
        GY(k) = 0.d0
        GZ(k) = 0.d0

        VTk(k) = 0.d0
        VTotal = 0.d0
        
        do i = 1, 2
         
          j = 2*(k-1) + i 
          m = j
          j = Id(j)
          GX(k)=GX(k)+zK(m)*FRX(X(k),Y(k),Z(k),X(j),Y(j),Z(j))
          GY(k)=GY(k)+zK(m)*FRY(X(k),Y(k),Z(k),X(j),Y(j),Z(j))
          GZ(k)=GZ(k)+zK(m)*FRZ(X(k),Y(k),Z(k),X(j),Y(j),Z(j))

          VTk(k)=VTk(k)+zK(m)*V(X(k),Y(k),Z(k),X(j),Y(j),Z(j))

        end do

        VTotal = VTotal + VTk(k)

      end do


      do k = 1, npar
        
        !write(32, *) Grdx(k)
        Grdx(k) = GX(k)
        Grdy(k) = GY(k)
        Grdz(k) = GZ(k)

      end do

      Epot = VTotal*0.5
      Epot = Epot/4.35988

      

      do k= 1, npar ! <---------
                                !
       fx= abs (Grdx(k))        !
       fy= abs (Grdy(k))        !
       fz= abs (Grdz(k))        !
                                !
       ff= sqrt (fx**2 + fy**2 + fz**2)
                                !
       if (fx.gt.xbig) then     ! x-component
       ixbig= k                 !
       xbig= fx                 !
       end if                   !
                                !
       if (fy.gt.ybig) then     ! y-component
       iybig= k                 !
       ybig= fy                 !
       end if                   !
                                !
       if (fz.gt.zbig) then     ! z-component
       izbig= k                 !
       zbig= fz                 !
       end if                   !
                                !
       if (ff.gt.big) then      ! magnitude
       ibig= k                  !
       big= ff                  !
       end if                   !
                                !
      end do ! <----------------

* reduce time steps for large forces

      if (big.gt.0.09) then
      dt= dtIni/2.
      else
      dt= dtIni
      end if

      write (*,*)  "largest positive forces in x,y,z"
      write (*,12) "for atoms:", ixbig, iybig, izbig
      write (*,14) "forces [Hartree/Bohr]", xbig, ybig, zbig
      write (*,12) "largest force magnitude for atom", ibig
      write (*,14) "force magnitude [Hartree/Bohr]", big

      if (dt.ne.dtIni) write (*,"(a,f9.2)") "the time step is", dt
      write (*,*)

   12 format (1x,a,3(1x,i4))
   14 format (1x,a,3(1x,f11.5))

      return

      end

* ***********************************
* get ab-initio energy & gradients
* ***********************************
      subroutine abEne (imethod,dtIni,dt,time,npar,
     &      Grdx,Grdy,Grdz,Epot,refPE)
      implicit double precision (a-h,o-z)
      include 'parametros'

      character qline*(maxchr)
      dimension itok(maxtok), jtok(maxtok)

      dimension xcrds(10)
      dimension Grdx(nmax), Grdy(nmax), Grdz(nmax)

      external kparse

* initialize variables

      xenergy= 0.d0
      Epot=  0.d0

      do i= 1, 10
      xcrds(i)= 0.d0
      end do

      do k= 1, npar
       Grdx(k)= 0.d0
       Grdy(k)= 0.d0
       Grdz(k)= 0.d0
      end do

* -------------------------
* find section of energy in TeraChem
* -------------------------

      if (imethod.eq.1) then    ! read TeraChem
      open (unit=16, file="terachem.out", status='unknown')

   10 read(16,20, end=50) qline ! <------
                                         !
      ntok= kparse(qline,maxchr,itok,jtok,maxtok)
                                         !
      if ((qline(itok(1):jtok(1)) .eq. 'FINAL') .and.
     &    (qline(itok(2):jtok(2)) .eq. 'ENERGY:')) then
                                         !
                          ! ---------------------------
                          ! ab-initio energy [Hartree]
                                         !
      read(qline(itok(3):jtok(3)),*) xenergy
                                         !
                                         ! take the energy of the first
      if (time.lt.0.d0) refPE= xenergy   ! geometry as the ref pot ene
                                         !
      Epot= (xenergy-refPE)              ! in Hartree
                                         !
      do k= 1, 25               ! junk files
      read(16,20) qline                  !
      ntok= kparse(qline,maxchr,itok,jtok,maxtok)
      if (qline(itok(1):jtok(1)).eq.'dE/dX') goto 15
      end do                             !
                                         !
   15 continue                           !
                                 ! --------------
                                 ! read gradients
                                         !
      do k= 1, npar               ! <----!--
      read(16,20) qline                  !
      ntok= kparse(qline,maxchr,itok,jtok,maxtok)
                                         !
      do j= 1, 3                         !
      read(qline(itok(j):jtok(j)),*) xcrds(j)
      end do                             !
                                         !
      Grdx(k)= xcrds(1)                  ! in Hartree/Bohr
      Grdy(k)= xcrds(2)                  !
      Grdz(k)= xcrds(3)                  !
      end do                      ! <----!--
                                         !
      end if                             !
      goto 10 ! <------------------------
      end if

* -------------------------
* find section of energy in NWChem
* -------------------------

      if (imethod.eq.2) then    ! read NWChem
      open (unit=16, file="nwchem.out", status='unknown')

   27 read(16,20, end=50) qline ! <------
                                         !
      ntok= kparse(qline,maxchr,itok,jtok,maxtok)
                                         !
      if ((qline(itok(1):jtok(1)) .eq. 'Total') .and.
     &    (qline(itok(2):jtok(2)) .eq. 'DFT')   .and.
     &    (qline(itok(3):jtok(3)) .eq. 'energy')) then
                                         !
                          ! ---------------------------
                          ! ab-initio energy [Hartree]
                                         !
      read(qline(itok(5):jtok(5)),*) xenergy
                                         !
                                         ! take the energy of the first
      if (time.lt.0.d0) refPE= xenergy   ! geometry as the ref pot ene
                                         !
      Epot= (xenergy-refPE)              ! in Hartree
      goto 18                            !
      end if                             !
                                         !
      goto 27 ! <------------------------

   18 continue

* -----------------------------------
* find section of energy gradients in NWChem
* -----------------------------------

   30 read(16,20, end=50) qline ! <------
                                         !
      ntok= kparse(qline,maxchr,itok,jtok,maxtok)
                                         !
      if ((qline(itok(1):jtok(1)) .eq. 'DFT') .and.
     &    (qline(itok(2):jtok(2)) .eq. 'ENERGY') .and.
     &    (qline(itok(3):jtok(3)) .eq. 'GRADIENTS')) then
                                         !
      do k= 1, 3                         ! junk files
      read(16,20) qline                  !
      end do                             !
                                         !
      do k= 1, npar               ! <----!--
      read(16,20) qline                  !
      ntok= kparse(qline,maxchr,itok,jtok,maxtok)
                                         !
      do j= 6, 8                         !
      read(qline(itok(j):jtok(j)),*) xcrds(j)
      end do                             !
                                         !
      Grdx(k)= xcrds(6)                  ! in Hartree/Bohr
      Grdy(k)= xcrds(7)                  !
      Grdz(k)= xcrds(8)                  !
      end do                      ! <----!--
                                         !
      goto 50                            !
      end if                             !
                                         !
      goto 30 ! <------------------------
      end if

   50 continue
      close (16)

* ---------------------------------
* obtain the largest forces in x,y,z
* ---------------------------------

      big=  -100.d0
      xbig= -100.d0
      ybig= -100.d0
      zbig= -100.d0

      do k= 1, npar ! <---------
                                !
       fx= abs (Grdx(k))        !
       fy= abs (Grdy(k))        !
       fz= abs (Grdz(k))        !
                                !
       ff= sqrt (fx**2 + fy**2 + fz**2)
                                !
       if (fx.gt.xbig) then     ! x-component
       ixbig= k                 !
       xbig= fx                 !
       end if                   !
                                !
       if (fy.gt.ybig) then     ! y-component
       iybig= k                 !
       ybig= fy                 !
       end if                   !
                                !
       if (fz.gt.zbig) then     ! z-component
       izbig= k                 !
       zbig= fz                 !
       end if                   !
                                !
       if (ff.gt.big) then      ! magnitude
       ibig= k                  !
       big= ff                  !
       end if                   !
                                !
      end do ! <----------------

* reduce time steps for large forces

      if (big.gt.0.09) then
      dt= dtIni/2.
      else
      dt= dtIni
      end if

      write (*,*)  "largest positive forces in x,y,z"
      write (*,12) "for atoms:", ixbig, iybig, izbig
      write (*,14) "forces [Hartree/Bohr]", xbig, ybig, zbig
      write (*,12) "largest force magnitude for atom", ibig
      write (*,14) "force magnitude [Hartree/Bohr]", big

      if (dt.ne.dtIni) write (*,"(a,f9.2)") "the time step is", dt
      write (*,*)

   12 format (1x,a,3(1x,i4))
   14 format (1x,a,3(1x,f11.5))

* ------------------------------
* broadcast problems of terachem
* ------------------------------

      if (xenergy.eq.0.0) then
        print*, "TERACHEM: got NO energy"
        print*, "possible lack of energy convergence"
        print*
c       stop
      end if

   20 format (a)

      return
      end
* *********************************
* ***********************************
* average the temp over the later ifreq steps
* ***********************************
      subroutine averageT (k,temp,ifreq,averTemp)
      implicit double precision (a-h,o-z)
      include 'parametros'

      dimension temp(maxt)

      icount= 0
      sum= 0.d0

      do i= k-ifreq, k
       icount= icount +1
       sum= sum + temp(i)
      end do

      averTemp= sum/dble(icount)

      return
      end
* ***********************************
* *****************************
* incident radiation upon a region
* stochastically breaking a couple of bonds
* more info in:
* http://www.epa.gov/radiation/understand/ionize_nonionize.html
* *****************************
      subroutine breakBond (npar,ncage,atsymbol,atomass,x,y,z)
      implicit double precision (a-h,o-z)
      include 'parametros'

      character atsymbol(nmax)*5
      dimension atomass(nmax)
      dimension x(nmax), y(nmax), z(nmax)

      integer itag(nmax)

      real*8 zbqlu01

* initialize labels

      icycle= 1

      do i= 1, npar  ! labels for breaking several
      itag(i)= 0     ! bonds stochastically
      end do

* ------------------------
* stochastically chose an atom
* ------------------------

   10 eta= zbqlu01()       ! zbqlu01 is in [0,1]
      nconf= npar-ncage    ! num of confined atoms
      n1= int(eta*nconf)   ! n1 is in [0,nconf]
      if (n1.eq.0) goto 10

      k= ncage + n1        ! this is the chosen atom

      if (itag(k).eq.1) goto 10        ! avoid repeated atoms
      if (atsymbol(k).eq."H") goto 10  ! avoid Hyd atoms

      write(*,"(a,f8.3)") "rnd number", eta
      itag(k)= 1           ! mark this atom

* ---------------------------
* find the two closest neighbors to atom k
* ---------------------------

      smalla= 1000.
      smallb= 1000.

      do i= ncage+1, npar ! <----  avoid cage atoms
                                 !
      if (i.ne.k) then           !
      v1= (x(k)-x(i))**2+ (y(k)-y(i))**2+ (z(k)-z(i))**2
      dist= sqrt(v1)             !
                                 !
      if (dist.lt.smalla) then   ! closest atom to K atom
      smalla= dist               !
      ismalla= i                 !
      end if                     ! also chose a
                                 ! heavy atom
      if ((dist.lt.smallb).and.(atsymbol(i).ne."H")) then
      smallb= dist               !
      ismallb= i                 !
      end if                     !
                                 !
      end if                     !
                                 !
      end do ! <-----------------

      j1= ismalla  ! this is the closest atom to k
      j2= ismallb  ! heavy atom closest to k

      itag(j1)= 1  ! mark these atoms
      itag(j2)= 1

* ------------------------
* find the closest neighbor to atom j2
* ------------------------

      small= 1000.

      do i= ncage+1, npar ! <----  avoid cage atoms
                                 !
      if (i.ne.j2) then          !
      v1= (x(j2)-x(i))**2+ (y(j2)-y(i))**2+ (z(j2)-z(i))**2
      dist= sqrt(v1)             !
                                 !
      if (dist.lt.small) then    ! closest atom
      small= dist                ! to j2
      ismall= i                  !
      end if                     !
                                 !
      end if                     !
      end do ! <-----------------

      j3= ismall   ! this is the closest atom to j2

      itag(j3)= 1  ! mark these atoms

* ------------------------
* eq of the line joining atoms j and k
* http://www.tec-digital.itcr.ac.cr/revistamatematica/cur
* sos-linea/Algebra-Lineal/algebra-vectorial-geova-walter/node5.html
* ------------------------

* points to use for the line
*   P= x(j),y(j),z(j)    ;    Q= x(k),y(k),z(k)

* director vector
*   vector(PQ)= x(k)-x(j), y(k)-y(j), z(k)-z(j)

* eq of the line joining P and Q
*   x,y,z= P + lambda * vector(PQ)

* -------------------------------------------
* break the molecular bond k-j1 weighting with atom masses
* -------------------------------------------

      dlam= 0.45 ! factor to increase the bond per atom

      a= x(k) +dlam *(x(k)-x(j1))/ atomass(k)
      b= y(k) +dlam *(y(k)-y(j1))/ atomass(k)
      c= z(k) +dlam *(z(k)-z(j1))/ atomass(k)

      d= x(j1) -dlam *(x(k)-x(j1))/ atomass(j1)
      e= y(j1) -dlam *(y(k)-y(j1))/ atomass(j1)
      f= z(j1) -dlam *(z(k)-z(j1))/ atomass(j1)

      write(*,"(a)") "1st broken bond"
      write(*,"(2(2x,a2,i4),a,f8.3)") atsymbol(k), k,
     &       atsymbol(j1), j1, " with distance", smalla

      dist1= sqrt((a-d)**2 + (b-e)**2 + (c-f)**2)
      write(*,"(a,f9.3)") "new distance between atoms", dist1

      dist2= sqrt((x(k)-a)**2 + (y(k)-b)**2 + (z(k)-c)**2)
      write(*,"(a,i4,1x,a,f9.3)") "increment for atom ", k,
     &    atsymbol(k), dist2

      dist3= sqrt((x(j1)-d)**2 + (y(j1)-e)**2 + (z(j1)-f)**2)
      write(*,"(a,i4,1x,a,f9.3)") "increment for atom ", j1,
     &    atsymbol(j1), dist3

      write(*,"(a,f9.3)") "total increment", dist2+dist3
      write(*,"(a,f9.3)") "increase factor", dist1/smalla

* assign new positions after opening the bond

      x(k)= a
      y(k)= b
      z(k)= c

      x(j1)= d
      y(j1)= e
      z(j1)= f

* -------------------------------------------
* in case that atoms j2 or j3 become k or j1 atoms
* re-assign coordinates
* -------------------------------------------

      if (j2.eq.k) then
      x(j2)= x(k)
      y(j2)= y(k)
      z(j2)= z(k)
      end if

      if (j2.eq.j1) then
      x(j2)= x(j1)
      y(j2)= y(j1)
      z(j2)= z(j1)
      end if

      if (j3.eq.k) then
      x(j3)= x(k)
      y(j3)= y(k)
      z(j3)= z(k)
      end if

      if (j3.eq.j1) then
      x(j3)= x(j1)
      y(j3)= y(j1)
      z(j3)= z(j1)
      end if

* re-compute the distance between j2-j3

      small= sqrt((x(j2)-x(j3))**2 +(y(j2)-y(j3))**2 +(z(j2)-z(j3))**2)

* -------------------------------------------
* break the molecular bond j2-j3 weighting with atom masses
* -------------------------------------------

      a= x(j2) +dlam *(x(j2)-x(j3))/ atomass(j2)
      b= y(j2) +dlam *(y(j2)-y(j3))/ atomass(j2)
      c= z(j2) +dlam *(z(j2)-z(j3))/ atomass(j2)

      d= x(j3) -dlam *(x(j2)-x(j3))/ atomass(j3)
      e= y(j3) -dlam *(y(j2)-y(j3))/ atomass(j3)
      f= z(j3) -dlam *(z(j2)-z(j3))/ atomass(j3)

      write(*,"(a)") "2nd broken bond"
      write(*,"(2(2x,a2,i4),a,f8.3)") atsymbol(j2), j2,
     &       atsymbol(j3), j3, " with distance", small

      dist1= sqrt((a-d)**2 + (b-e)**2 + (c-f)**2)
      write(*,"(a,f9.3)") "new distance between atoms", dist1

      dist2= sqrt((x(j2)-a)**2 + (y(j2)-b)**2 + (z(j2)-c)**2)
      write(*,"(a,i4,1x,a,f9.3)") "increment for atom ", j2,
     &    atsymbol(j2), dist2

      dist3= sqrt((x(j3)-d)**2 + (y(j3)-e)**2 + (z(j3)-f)**2)
      write(*,"(a,i4,1x,a,f9.3)") "increment for atom ", j3,
     &    atsymbol(j3), dist3

      write(*,"(a,f9.3)") "total increment", dist2+dist3
      write(*,"(a,f9.3)") "increase factor", dist1/small

      write(*,"(a,i3,a,i3,f9.3)") "distance between atoms ", k,
     &     " and ", j2, smallb 

* assign new positions after opening the bond

      x(j2)= a
      y(j2)= b
      z(j2)= c

      x(j3)= d
      y(j3)= e
      z(j3)= f

* -----------------------------
* write info of the molecular bond
* -----------------------------

c     write(*,"(a,3f9.3)") atsymbol(k),  x(k),  y(k),  z(k)
c     write(*,"(a,3f9.3)") atsymbol(j1), x(j1), y(j1), z(j1)

      icycle= icycle +1
c     if (icycle.le.1) goto 10    ! in case of breaking several bonds

      return
      end
* *****************************
* ***************************************
* get cartesian coordinates (a,b,c) from
* the polar coordinates (rho,teta,fi)
* ***************************************
      subroutine cartesian (rho,teta,fi,a,b,c)
      implicit double precision (a-h,o-z)

      pi= 3.1415926535898
      rad= pi/180.0

*  0 <= rho
*  0 <= teta <= 2pi
*  0 <= fi <= pi

      a= rho* sin(fi)* cos(teta)
      b= rho* sin(fi)* sin(teta)
      c= rho* cos(fi)

      return
      end
* ****************************
* ***********************************
* citations of IMD, TeraChem, etc.
* ***********************************
      subroutine citation (imethod)
      implicit double precision (a-h,o-z)
      include 'parametros'

      print*, "--------------------------------------------"
      print*, "  CITATION:"
      print*

      print*, "Interactive Molecular Dynamics Program"
      if (imethod.eq.1) then
      print*, "Module: Pressure Cooker (TeraChem version)"
      end if

      if (imethod.eq.2) then
      print*, "Module: Pressure Cooker (NWChem version)"
      end if

      print*, "Ruben Santamaria:   rso@fisica.unam.mx"
      print*, "IFUNAM, 2014"
      print*

c     print*, "TeraChem"
c     print*, "I.S. Ufimtsev and T.J. Martinez"
c     print*, "Quantum Chemistry on Graphical Processing Units. 3."
c     print*, "Analytical Energy Grad. & First Princ. Molec. Dyn."
c     print*, "Jour. Chem. Theory & Comput. Vol. 5, p2619, 2009"
c     print*

c     print*, "random-number-generator module from:"
c     print*, "Richard Chandler & Paul Northrop"
c     print*, "richard@stats.ucl.ac.uk, northrop@stats.ox.ac.uk"
c     print*

c     print*, "I, Ruben Santamaria, have taken every precaution for"
c     print*, "the writing of this program. Although I have made my"
c     print*, "best effort, it is not guaranteed to be free of error."
c     print*, "Therefore, it should not be relied on for solving"
c     print*, "problems where an error could result in injury or loss."
c     print*, "The author disclaims all liability in this regard."
c     print*, "Ruben Santamaria retains the rights to make changes to"
c     print*, "this program at any time, without previous notice. The"
c     print*, "program can be used, copied, modified and distributed"
c     print*, "under the condition of giving the acknowledgments to"
c     print*, "to the author, as indicated in the citation section"
c     print*, "of the program."
c     print*

      print*, "wall-clock time"
      call system ("date")

      return
      end
* ********************
* ************************************
* measure distances between particles
* ************************************
      subroutine distances (time, npar,x,y,z)
      implicit double precision (a-h,o-z)
      include 'parametros'

      dimension x(nmax), y(nmax), z(nmax)

* ------------------------------------
* find the smallest & largest distances
* ------------------------------------

      small= 1.d+12
      big=  -1.d+12

      do i= 1, npar ! <----------
      do j= 1, npar              !
                                 !
      if (i.lt.j) then ! <-------!------
                                 !
      v1= (x(i)-x(j))**2 + (y(i)-y(j))**2 + (z(i)-z(j))**2
      Rij= sqrt (v1)             !
                                 !
      if (Rij.lt.small) then     ! get the smallest
      small= Rij                 ! distance
      ismall= i                  !
      jsmall= j                  !
      end if                     !
                                 !
      if (Rij.gt.big) then       ! get the largest
      big= Rij                   ! distance
      ibig= i                    !
      jbig= j                    !
      end if                     !
                                 !
      end if ! <-----------------!------
                                 !
      end do                     !
      end do !<------------------

c     write(*,20) "small", time, ismall, jsmall, small
c     write(*,20) "big  ", time, ibig,   jbig,   big

   20 format (1x,a,f8.3,2i5,f8.3)

      return
      end
* ***********************************
* ***********************************
* get the potential energy of the springs
* ***********************************
      subroutine eneSpring (time,npar,ncage,x,y,z,x0,y0,z0,sK,Espring)
      implicit double precision (a-h,o-z)
      include 'parametros'

      dimension x0(nmax), y0(nmax), z0(nmax), sK(nmax)
      dimension x(nmax),  y(nmax),  z(nmax)

* ------------------------
* potential energy of springs
* ------------------------

      Espring= 0.d0
      sum= 0.d0

      do i= 1, ncage
      v1= (x(i)-x0(i))**2 +(y(i)-y0(i))**2 +(z(i)-z0(i))**2
      sum= sum + sK(i) * v1
      end do

      Espring= sum/2.  ! in (N/cm)*Angst^2

* Angst^2/cm= 1x10^{-20} m^2/ 1x10^{-02} m= 1x10^{-18} m

* Espring= Espring* 1.0d-18  ! in N*m=Joules

* 1 Hartree= 4.35988 x10^{-18} J

      Espring= Espring/ 4.35988   ! in Hartree

      return
      end
* *********************************
* ***********************************
* estimate CPU time from the first computations
* ***********************************
      subroutine estimateCPU (value)
      implicit double precision (a-h,o-z)
      include 'parametros'

      character qline*(maxchr)
      dimension itok(maxtok), jtok(maxtok)

      external kparse

      open (unit=16, file="terachem.out", status='unknown')

* --------------------------
* find the section of interest
* --------------------------

   10 read(16,20, end=50) qline ! <------
                                         !
      ntok= kparse(qline,maxchr,itok,jtok,maxtok)
                                         !
      if ((qline(itok(1):jtok(1)) .eq. 'Total') .and.
     &    (qline(itok(2):jtok(2)) .eq. 'processing') .and.
     &    (qline(itok(3):jtok(3)) .eq. 'time:')) then
                                         !
c     print*, "qline= ", qline           !
                                         ! read the time
      read(qline(itok(4):jtok(4)),*) value
                                         !
c     print*, "valor encontrado", value  !
                                         !
      end if                             !
      goto 10 ! <------------------------

   50 continue
   20 format (a)

      close (16)

      return
      end
* *********************************
* *********************
* draw a figure when the job is finished
* *********************
      subroutine figures

      write (*,*) "figure from: http://www.asciiworld.com"
      write (*,*)
      write (*,*) "          _,.-------.,_          "
      write (*,*) "      ,;~'             '~;,      "
      write (*,*) "    ,;                     ;,    "
      write (*,*) "   ;                         ;   "
      write (*,*) "  ,'                         ',  "
      write (*,*) " ,;                           ;, "
      write (*,*) " ; ;      .           .      ; ; "
      write (*,*) " | ;   ______       ______   ; | "

c     write (*,*) " |  `/~"     ~" . "~     "~\'  | "
      write (*,*) " |  `/~","     ~", char(34), "  .  ", char(34),
     &     "~          |"

      write (*,*) " |  ~  ,-~~~^~, | ,~^~~~-,  ~  | "
      write (*,*) "  |   |        }:{        |   |  "
      write (*,*) "  |   l       / | \       !   |  "

c     write (*,*) "  .~  (__,.--" .^. "--.,__)  ~.  "
      write (*,*) "  .~  (__,.--", char(34), " .^. ", char(34),
     & "--.,__)  ~."

      write (*,*) "  |     ---;' / | \ ';---     |  "
      write (*,*) "   \__.       \/^\/       .__/   "
      write (*,*) "    V| \                 / |V    "
      write (*,*) "     | |T~\___!___!___/~T| |     "
      write (*,*) "     | |`IIII_I_I_I_IIII'| |     "
      write (*,*) "     |  \,III I I I III,/  |     "
      write (*,*) "      \   `~~~~~~~~~~'    /      "
      write (*,*) "        \   .       .   /        "
      write (*,*) "          \.    ^    ./          "
      write (*,*) "            ^~~~^~~~^            "
      write (*,*)
      write (*,*) "         THIS IS THE END"
      write (*,*) "* ----------------------------- *"
      write (*,*)

      return
      end
* *********************

* ******************************
* collect charges of atom groups
* ******************************
      subroutine groupCharge (time, ngrp1,ngrp2,ngrp3)

      implicit real*8 (a-h,o-z)
      include 'parametros'

      character*3 name

* open file to read charges

      open (unit=13, file= './scr/charge_mull.xls', status= 'unknown')

* ------------------------------
* get charges of groups of atoms
* ------------------------------

      sum1= 0.d0
      sum2= 0.d0
      sum3= 0.d0

      do i= 1, ngrp1
      read (13,*) j, name, charge
      sum1= sum1 + charge
      end do

      do i= ngrp1+1, ngrp2
      read (13,*) j, name, charge
      sum2= sum2 + charge
      end do

      do i= ngrp2+1, ngrp3
      read (13,*) j, name, charge
      sum3= sum3 + charge
      end do

      write (*,12) time, sum1,sum2,sum3, sum1+sum2+sum3

      close (13)

   12 format (1x,"t= ",f8.3," charges=", 3(1x,f12.4),
     &   " total=",1x,f12.4)

      return
      end
* ************************
* *****************************
* gather info of the user, machine, etc
* *****************************
      subroutine infoUser
      implicit double precision (a-h,o-z)
      include 'parametros'

      character pwd*60

* ------------------------------------------
* identify initial cpu processes & directories
* ------------------------------------------

      print*
      call system ("whoami")
      call system ("hostname")

      call system ("pwd > t1")
      call workDir (pwd)
      write (*,"(3a)") "working directory ", pwd

c     print*, "wall-clock time"
      call system ("date")
c     call system ("date -R | cut -c1-25")

      print*, "fortran process ID   ", getpid()
      print*, "ab-initio process ID ", getpid()+3
      print*, "user ID  ", getuid()
      print*, "group ID ", getgid()

      call system ("rm t1 t2")
c     call system ("rm -rf permanent scratch")
c     call system ("mkdir permanent")
c     call system ("mkdir scratch")

      return
      end
* *****************************

* ***********************************
* impose initial conditions on positions and speeds
* ***********************************
      subroutine initConditOlla (npar,ncage, atomass, x,y,z, vx,vy,vz,
     &    x0,y0,z0,sK, cmx,cmy,cmz, tempEquil,
     &    nsteps,ibTime1,ibTime2,ibTime3)

      implicit double precision (a-h,o-z)
      include 'parametros'

      external gauss

      dimension atomass(nmax)
      dimension  x(nmax),  y(nmax),  z(nmax)
      dimension vx(nmax), vy(nmax), vz(nmax)

      dimension x0(nmax), y0(nmax), z0(nmax), sK(nmax)

      integer iBreakTime(3)
      real*8 zbqlu01, kB, iniTemp

* ------------------
* define constants
* ------------------

* Boltzmann constant

      kB= 1.3806504d-23  ! in Joule/Kelvin= Kg*(m/s)^2/Kelvin

* facKT *(Angst/fs)^2= (Joule/Kelvin)*(Kelvin/amu)

      facKT= 1.0d17/1.6605402

      iniTemp= tempEquil

* ------------------
* loop over the particles
* ------------------

      print*, "moving the origin to the CM of the cage"
      print*

      write(*,12) "required temp for the simulation", iniTemp
      print*, "giving maxwellian velocities to atoms"

      do i= 1, npar ! <----
                           !
                 ! -----------------
                 ! initial positions
                           !
         x(i)= x(i)-cmx    ! move the origin
         y(i)= y(i)-cmy    ! to the cage CM
         z(i)= z(i)-cmz    !
                           !
         x0(i)= x(i)       ! store initial
         y0(i)= y(i)       ! positions
         z0(i)= z(i)       !
                           !
                   ! --------------
                   ! initial speeds
                           !
      pKT= kB*iniTemp/ atomass(i)
                           !
      etax= gauss()        ! random numbers
      etay= gauss()        !
      etaz= gauss()        ! Maxwellian speeds
                           !
      vx(i)= etax *sqrt(pKT*facKT) ! in Angst/fs
      vy(i)= etay *sqrt(pKT*facKT)
      vz(i)= etaz *sqrt(pKT*facKT)
                           !
      end do ! <-----------

* ----------------------------
* get the speed of the mass center
* ----------------------------

      sumx= 0.d0
      sumy= 0.d0
      sumz= 0.d0
      sumass= 0.d0

      do i= ncage+1, npar   ! for the confined atoms
       sumass= sumass + atomass(i)
       sumx= sumx + atomass(i) *vx(i)
       sumy= sumy + atomass(i) *vy(i)
       sumz= sumz + atomass(i) *vz(i)
      end do

      v1= sumx/ sumass
      v2= sumy/ sumass
      v3= sumz/ sumass

* ---------------------------------
* refer the speeds to the speed of the CM
* ---------------------------------

      do i= ncage+1, npar ! <--
                               !
        vx(i)= vx(i)-v1        ! for the
        vy(i)= vy(i)-v2        ! confined
        vz(i)= vz(i)-v3        ! atoms
                               !
      end do ! <---------------

      print*, "making zero the speed of the CM of confined atoms"

* --------------------
* kinetic energy of confined atoms
* --------------------

      Ekin= 0.d0
      nn= 0    ! momentum degrees of freedom

      do i= ncage+1, npar ! <-------
                                    ! in amu*(Angst/fs)^2
       nn= nn + 1                   !
       v2= vx(i)**2 + vy(i)**2 + vz(i)**2
       Ekin= Ekin + atomass(i)*v2/ 2.d0
                                    !
      end do   ! <------------------

* facEkin * Joule= amu *(Angst/fs)^2

      facEkin= 1.6605402d-17

      EkinConf= Ekin*facEkin    ! in Joules

* compute the temperature from:  Ekin=3NkT/2 ==> T=2Ekin/(3Nk)

      tempConf= 2.*EkinConf/ (3.*dble(nn)*kB)  ! in Kelvin

* --------------------
* kinetic energy of cage atoms
* --------------------

      Ekin= 0.d0
      nn= 0    ! momentum degrees of freedom

      do i= 1, ncage ! <---------
                                 ! in amu*(Angst/fs)^2
       nn= nn + 1                !
       v2= vx(i)**2 + vy(i)**2 + vz(i)**2
       Ekin= Ekin + atomass(i)*v2/ 2.d0
                                  !
      end do   ! <----------------

      EkinCage= Ekin*facEkin    ! in Joules

      tempCage= 2.*EkinCage/ (3.*dble(nn)*kB)  ! in Kelvin

      write(*,12) "initial temp of cage atoms [K]", tempCage
      write(*,12) "initial temp of confined atoms [K]", tempConf

   12 format (1x,a,f11.3)

* ---------------------------
* chose the times to break bonds in stochastic form
* ---------------------------

      do k= 1, 3 ! <---------
                             !
      iBreakTime(k)= 0       ! default value
                             !
   10 eta= zbqlu01()         ! zbqlu01 is in [0,1]
      n1= int(eta*nsteps)    ! n1 is in [0,nsteps]
                             !
      if ((n1.lt.100).or.(n1.gt.nsteps-100)) goto 10
      iBreakTime(k)= n1      ! this break time is taken
                             !
      end do ! <-------------

      n1= iBreakTime(1)
      n2= iBreakTime(2)
      n3= iBreakTime(3)

* find the order of the numbers n1,n2,n3

      if ((n1.lt.n2).and.(n1.lt.n3)) then
      ibTime1= n1
         if (n2.lt.n3) then
         ibTime2= n2
         ibTime3= n3   ! the order is n1<n2<n3
         else
         ibTime2= n3
         ibTime3= n2   ! the order is n1<n3<n2
         end if
      end if

      if ((n2.lt.n1).and.(n2.lt.n3)) then
      ibTime1= n2
         if (n1.lt.n3) then
         ibTime2= n1
         ibTime3= n3   ! the order is n2<n1<n3
         else
         ibTime2= n3   
         ibTime3= n1   ! the order is n2<n3<n1
         end if
      end if

      if ((n3.lt.n1).and.(n3.lt.n2)) then
      ibTime1= n3
         if (n1.lt.n2) then
         ibTime2= n1
         ibTime3= n2   ! the order is n3<n1<n2
         else
         ibTime2= n2
         ibTime3= n1   ! the order is n3<n2<n1
         end if
      end if

* a minimum difference of 10 is required between
* consecutive numbers to give chance to the charge uptake

      if (ibTime2-ibTime1.le.10) then
      ibTime2= ibTime2+13   ! increase the value of ibTime2
      end if

      if (ibTime3-ibTime2.le.10) then
      ibTime3= ibTime3+13   ! increase the value of ibTime3
      end if

      write(*,"(/,1x,a,3i7)") "rnd times to break bonds",
     &   ibTime1, ibTime2, ibTime3

* ---------------------------
* make spring constants different
* ---------------------------

      do i= 1, ncage ! <-----
                             !
      if (i.eq.1) write(*,*) "spring constants were made different"
                             !
      eta= zbqlu01()         ! zbqlu01 is in [0,1]
                             !
      sK(i)= sK(1) + eta/10. ! maximum difference
                             ! up to 0.1 unit
      end do ! <-------------

      return
      end
* *******************************
* ***********************************
* create the nwchem input file
* ***********************************
      subroutine inptNwchem (ilabel,npar,atsymbol,x,y,z)
      implicit double precision (a-h,o-z)
      include 'parametros'

      character atsymbol(nmax)*5
      dimension x(nmax), y(nmax), z(nmax)

* ilabel= 1    generate the wavefunction (WF) from scratch
* ilabel= 2    reuse the WF from a previous computation
* ilabel= 3    print the nwchem input to the output file
* ilabel= 4    change the charge of the molecule & generate the WF
* ilabel= 5    change the charge of the molecule & reuse the WF

* -------------------------
* chose the appropriate run
* -------------------------

      if ((ilabel.eq.1).or.(ilabel.eq.4)) then
        open (14,file="nwchem.inp",status="unknown")
        write(*,"(/,a,i5)") "generating the wave function from scratch",
     &  ilabel
        ij= 14
      end if

      if ((ilabel.eq.2).or.(ilabel.eq.5)) then
        open (14,file="nwchem.inp",status="unknown")
        write(*,"(/,a,i5)") "reusing the wave function", ilabel
        ij= 14
      end if

      if (ilabel.eq.3) then
        print*, "the ab-initio input file may be"
        print*, "changed in module inptNwchem"
        ij= 6
      end if

* ---------------------
* write basic parameters
* ---------------------

      write(ij,*)
      write(ij,*) "title ", char(34), "mol computation", char(34)
      write(ij,*) "start molec"
c     write(ij,*) "echo"

      write(ij,*)

      write(ij,*) "memory total 4000 stack 1000 heap 1000 global 2000
     & mb"
      write(ij,*) "permanent_dir /home/rso/work/permanent"
      write(ij,*) "scratch_dir   /home/rso/work/scratch"

      write(ij,*)

      if (ilabel.eq.4) print*, "charge changed from -1 to 0"

      if (ilabel.eq.1) write(ij,*) "charge  0"  ! starting charge
      if (ilabel.eq.2) write(ij,*) "charge  0"
      if (ilabel.eq.4) write(ij,*) "charge  0"  ! changed charge
      if (ilabel.eq.5) write(ij,*) "charge  0"

* --------------
* write mol coords
* --------------

      write(ij,*)
      write(ij,*) "geometry units angstroms noautoz noautosym nocenter"

      if (ilabel.ne.3) then ! <---
      do i= 1, npar               !
      write(ij,"(2x,a,3f9.3)") atsymbol(i), x(i),y(i),z(i)
      end do                      !
      end if ! <------------------

      write(ij,*) "end"

* --------------
* write basis sets
* --------------

      write(ij,*)
      write(ij,*) "basis"
      write(ij,*) "  H library 6-31g"
      write(ij,*) "  C library 6-31g"
      write(ij,*) "  N library 6-31g"
      write(ij,*) "  O library 6-31g"
c     write(ij,*) "  S library 6-31g"
      write(ij,*) " He library 6-31g"
c     write(ij,*) " Pd library ", char(34), "LANL2DZ ECP", char(34)
      write(ij,*) "end"

c     write(ij,*)
c     write(ij,*) "set basis ", char(34), "cd basis", char(34)
c     write(ij,*) "basis ", char(34), "cd basis", char(34)
c     write(ij,*) "  H library 6-31g"
c     write(ij,*) "  C library 6-31g"
c     write(ij,*) "  N library 6-31g"
c     write(ij,*) "  O library 6-31g"
c     write(ij,*) "  S library 6-31g"
c     write(ij,*) " He library 6-31g"
c     write(ij,*) " Pd library ", char(34), "LANL2DZ ECP", char(34)
c     write(ij,*) "end"

c     write(ij,*)
c     write(ij,*) "ecp"
c     write(ij,*) " Pd library ", char(34), "LANL2DZ ECP", char(34)
c     write(ij,*) "end"

* --------------------
* write level of theory
* --------------------

      write(ij,*)

      write(ij,*) "dft"
      write(ij,*) " xc becke88 lyp"

      if (ilabel.eq.4) print*, "multiplicity changed from 2 to 1"

      if (ilabel.eq.1) write(ij,*) " mult 1"  ! starting multiplicity
      if (ilabel.eq.2) write(ij,*) " mult 1"
      if (ilabel.eq.4) write(ij,*) " mult 1"  ! changed multiplicity
      if (ilabel.eq.5) write(ij,*) " mult 1"

      write(ij,*) " iterations 1500"
      write(ij,*) " convergence energy  1e-5"
      write(ij,*) " convergence density 1e-5"
      write(ij,*) " grid medium"
      write(ij,*) " mulliken"

      if ((ilabel.eq.2).or.(ilabel.eq.5)) then
      write(ij,*)
     &  " vectors input /home/rso/work/permanent/molec.movecs"
      end if

      write(ij,*) "end"
      write(ij,*)

* ------------------
* write type of task
* ------------------

      write(ij,*) "task dft gradient"
      write(ij,*)

      if (ij.eq.14) close (ij)

      return
      end
* ********************
* ***********************************
* create the terachem input file
* ***********************************
      subroutine inptTera (ilabel)
      implicit double precision (a-h,o-z)

* ilabel= 1    generate the wavefunction (WF) from scratch
* ilabel= 2    reuse the WF from a previous computation
* ilabel= 3    print the terachem input to the output file
* ilabel= 4    change the charge of the molecule & generate the WF
* ilabel= 5    change the charge of the molecule & reuse the WF

* -------------------------
* chose the appropriate run
* -------------------------

      if (ilabel.eq.3) then
        print*, "the ab-initio input file may be"
        print*, "changed in module inptTera"
        ij= 6
      else
        open (14,file="terachem.inp",status="unknown")
        write(*,"(/,a,i2)") "WF type", ilabel
        ij= 14
      end if

* ---------------------
* write basic parameters
* ---------------------

      write (ij,*)

      if (ilabel.eq.4) print*, "charge changed from -1 to 0"
      if (ilabel.eq.4) print*, "multiplicity changed from 2 to 1"

      if (ilabel.eq.1) write(ij,*) "charge -1"  ! starting charge
      if (ilabel.eq.1) write(ij,*) "spinmult 2"
      if (ilabel.eq.2) write(ij,*) "charge -1"
      if (ilabel.eq.2) write(ij,*) "spinmult 2"

      if (ilabel.eq.4) write(ij,*) "charge  0"  ! changed charge
      if (ilabel.eq.4) write(ij,*) "spinmult 1"
      if (ilabel.eq.5) write(ij,*) "charge  0"
      if (ilabel.eq.5) write(ij,*) "spinmult 1"

      write (ij,*) "coordinates mol.xyz"
      write (ij,*) "units angstrom"
      write (ij,*) "method roblyp"
c     write (ij,*) "method ublyp"
      write (ij,*) "basis 6-31g"
      write (ij,*) "precision mixed"
      write (ij,*) "convthre 3.0e-5"
      write (ij,*) "maxit 1000"
      write (ij,*) "run gradient"
      write (ij,*) "gpus 4 0 1 2 3"  ! run 2 gpus, namely 0 & 3

* generate the wave function from scratch

      if (ilabel.eq.1)  write (ij,*) "guess generate"
      if (ilabel.eq.2)  write (ij,*) "guess ./scr/ca0   ./scr/cb0"
      if (ilabel.eq.4)  write (ij,*) "guess generate"
      if (ilabel.eq.5)  write (ij,*) "guess ./scr/ca0"

      write (ij,*) "end"
      write (ij,*)

      if (ij.eq.14) close (ij)

      return
      end
* ********************
* **********************************
* compute kinetic energies & temperature
* **********************************
      subroutine kinEne (npar,ncage, atomass, vx,vy,vz,
     &     EkinConf,EkinCage, tempConf,tempCage,Vol,dynPress)

      implicit double precision (a-h,o-z)
      include "parametros"

      dimension atomass(nmax)
      dimension vx(nmax), vy(nmax), vz(nmax)

      real*8 kB, mdof

* Boltzmann constant

      kB= 1.3806504d-23  ! in Joule/Kelvin= Kg*(m/s)^2/Kelvin

* --------------------
* kinetic energy of confined atoms
* --------------------

      Ekin= 0.d0
      mdof= 0.   ! momentum degrees of freedom

      do i= ncage+1, npar ! <-------
                                    ! in amu*(Angst/fs)^2
       mdof= mdof + 1.              !
       v2= vx(i)**2 + vy(i)**2 + vz(i)**2
       Ekin= Ekin + atomass(i)*v2/ 2.d0
                                    !
      end do   ! <------------------

* facEkin * Joule= amu *(Angst/fs)^2

      facEkin= 1.6605402d-17

      EkinConf= Ekin*facEkin    ! in Joules

* compute the temperature from:  Ekin=3NkT/2 ==> T=2Ekin/(3Nk)

      tempConf= 2.*EkinConf/ (3.*mdof*kB)  ! in Kelvin

* --------------------
* kinetic energy of cage atoms
* --------------------

      Ekin= 0.d0
      mdof= 0.   ! momentum degrees of freedom

      do i= 1, ncage ! <---------
                                 ! in amu*(Angst/fs)^2
       mdof= mdof + 1.           !
       v2= vx(i)**2 + vy(i)**2 + vz(i)**2
       Ekin= Ekin + atomass(i)*v2/ 2.d0
                                  !
      end do   ! <----------------

      EkinCage= Ekin*facEkin    ! in Joules

      tempCage= 2.*EkinCage/ (3.*mdof*kB)  ! in Kelvin

* Francis W. Sears \& Gerhard L. Salinger,
* {\it Termodin\'amica, teor\'ia cin\'etica, y termodin\'amica estad\'istica},
* 2nd ed., Revert\'e S. A. M\'exico, 1978, page 300.
* http://www.seykota.com/rm/gas/gas.htm

* --------------------
* dynamic pressure from the ideal gas:
* PV= NkT= Nk [2Ekin/(3Nk)]= 2Ekin/3
* refer to:
* http://hyperphysics.phy-astr.gsu.edu/hbase/kinetic/idegas.html
* --------------------

* from above T=2Ekin/(3Nk)  ==> PV=NkT=Nk(2Ekin/(3Nk))=2Ekin/3

* facPres2 * Pa= Joule/Angst^3

      facPres2= 1.0d30

      dynPress= 2./3.* EkinConf/Vol *facPres2 ! in Pascals

* 1 GPa= 1x10^{9} Pa
c     dynPress= dynPress/1.d09  ! GPa

* 1 pascal= (1/101325.0) atmosphere

      dynPress= dynPress/ 101325.0  ! atmospheres

* 1 Hartree= 4.35988 x 10^{-18} J

      EkinConf= EkinConf/4.35988d-18   ! in Hartree
      EkinCage= EkinCage/4.35988d-18   ! in Hartree

      return
      end
* *********************************
* ***********************************
* this part illustrates the call
* to the function kparse
* ***********************************
c     subroutine test
c     implicit double precision (a-h,o-z)
c     parameter (maxchr=128, maxtok=50, m1=200)

c     character qline*(maxchr)
c     dimension itok(maxtok),jtok(maxtok)
c     external kparse
c     external lenstr

* start debugging lines

c  10 read(11,900, end=50) qline
c     ntok= kparse(qline,maxchr,itok,jtok,maxtok)

c     if ((qline(itok(1):jtok(1)) .eq. 'DFT') .and.       ! 11111
c    &    (qline(itok(2):jtok(2)) .eq. 'ENERGY') .and.
c    &    (qline(itok(3):jtok(3)) .eq. 'GRADIENTS')) then

c     palabra= qline(itok(3):jtok(3))
c     nn= lenstr (palabra, maxchr)
c     print*, "longitud de la palabra ", palabra, "es", nn
c     print*, "longitud de la palabra ", jtok(3)-itok(3)+1

c  50 write(*,*) ' '
c 900 format(a)

* ************************************
* parte en palabras a la linea "qline"
* regresando los indices itok y jtok.
* ************************************
      function kparse(qline,lenq,itok,jtok,maxtok)
      implicit double precision (a-h,o-z)
      dimension itok(*),jtok(*)
      character qline*(*)
      character qspace*1,qtab*1
 
* qline(itok(1):jtok(1) <-- 1st word
* itok(1)   <-- first letter of 1st word
* jtok(1)   <-- last letter of 1st word

* qline(itok(2):jtok(2) <-- 2nd word
* itok(2)   <-- first letter of 2nd word
* jtok(2)   <-- last letter of 2nd word

      kparse= 0
      kount= 0
 
      qspace= ' '
      qtab= char(9) ! refers to a tab
 
* --------------------------------------------
* look for the first (non-blank or non-tab) token
* --------------------------------------------

   10 kount= kount + 1 ! <---------------------------
                                                     !
      if (kount.gt.lenq) goto 40                     !
                                                     !
* we are not interested on spaces or tabs            !
                                                     !
      if (qline(kount:kount) .eq. qspace .or.        !
     &    qline(kount:kount) .eq. qtab) goto 10 ! <--
 
* we've found a new token, check for maximum token count
 
      if (kparse .eq. maxtok) goto 40
 
* mark the beginning
 
      kparse= kparse + 1  ! label for the word
      itok(kparse)= kount ! first letter of the word

* -------------------------------
* look for the end of this token
* -------------------------------
 
   20 kount= kount + 1 ! <----------------------------
                                                      !
      if (kount .gt. lenq) goto 30                    !
                                                      !
      if (qline(kount:kount) .ne. qspace .and.        !
     &    qline(kount:kount) .ne. qtab) goto 20 ! <---
 
* end of the token, mark it
 
   30 jtok(kparse)= kount-1 ! last letter of the word

* ------------------------------- 
* continue searching for words
* ------------------------------- 
 
      if (kount.lt.lenq) goto 10
 
   40 return
      end
* ********************************
* computes the length of a string
* ********************************
      function lenq(qstring,max)
      implicit double precision (a-h,o-z)
      character qstring*(*)
 
      lenq= max

* find a blank space

   10 if (qstring(lenq:lenq).ne.' ') return

      lenq = lenq - 1
      if (lenq.le.0) return

      goto 10
 
      end
* ********************************
* computes the length of a string
* ********************************
      integer function lenstr(string, max)
      implicit none
      character*(*) string
      integer max

      lenstr= max

* find a blank space

   10 if (string(lenstr:lenstr).ne.' ') return
      print*,"aqui stoy ", string, lenstr
      print*,"aqui stoy 123456789012345678901234567890"

      lenstr= lenstr - 1

* if the lenght is 0 return

      if (lenstr.le.0) return

      goto 10

      return
      end
* ********************************
* ************************************
* get the effective radius of the confined particles
* ************************************
      subroutine largest (t,npar,ncage, x,y,z, radConf,radCage, Vol)
      implicit double precision (a-h,o-z)
      include 'parametros'

      dimension x(nmax), y(nmax), z(nmax)
      dimension dist(nmax)

      integer itag(nmax), jtag(nmax)

      pi= 3.1415926535898

* -----------------
* initialize tags
* -----------------

      knt= 0

      do j= ncage+1, npar ! <----
       knt= knt + 1              !
                                 !
       itag(knt)= 0              !
       jtag(knt)= 0              ! radial distances
                                 !
       dist(knt)= sqrt(x(j)**2 +y(j)**2 +z(j)**2)
                                 !
      end do ! <-----------------

      idim= knt

* ---------------------------------
* find the largest numbers of the array
* ---------------------------------

      do i= 1, idim ! <---------
                                !
      max= 0                    !
      dmax= 0.d0                !
                                !
* get the largest distance without a label
                                !
      do j= 1, idim             !
        if ((dist(j).ge.dmax).and.(itag(j).eq.0)) then
        dmax= dist(j)           !
        max= j                  !
        end if                  !
      end do                    !
                                !
      itag(max)= 1 ! switch to avoid examination of this num
      jtag(i)= max ! store the index of the largest number
                                !
      end do ! <----------------

* jtag(1)    ! <---- index of the largest num
* jtag(2)
* ......
* jtag(idim) ! <---- index of the smallest num

* ---------------------------
* radius of the confined particles
* ---------------------------

* 30% of atoms determine the confinement radius

      nn= (npar-ncage)/10 *3

* sum radial distances of most external atoms

      knt= 0
      sum= 0.d0

      do i= 1, nn ! <-----
                          !
       knt= knt + 1       !
       sum= sum + dist(jtag(i))
                          !
c      print*, i, jtag(i), dist(jtag(i))
                          !
      end do ! <----------

      radPart= sum /dble(knt)    ! average radius

      if (t.eq.-1.0) then
      write(*,14) "num of atoms defining the confinement radius",knt
      end if

* ------------------------
* get the effective confinement radius
* ------------------------

* radius of the cage

      sum= 0.d0

      do i= 1, ncage
      sum= sum +sqrt (x(i)**2 +y(i)**2 +z(i)**2)
      end do

      radCage= sum /dble(ncage) ! cage radius

* effective confinement radius

      radConf= radPart +(radCage-radPart)/ 2.

      area=   4.*pi* radConf**2    ! in Angst^2
      Vol= 4./3.*pi *radConf**3    ! in Angst^3

   14 format (1x,a,1x,i5)

      return
      end
* ************************************
* ***********************************
* make a movie from the atom coords
* ***********************************
      subroutine makeMovie (iframe,npar,atsymbol,x,y,z)
      implicit double precision (a-h,o-z)
      include 'parametros'

      character atsymbol(nmax)*5
      dimension x(nmax), y(nmax), z(nmax)

* -------------------------
* write coords for a movie
* -------------------------

      if (mod(iframe,1).eq.0) then     !
c     
      !write(12,*) npar + 1
      write(12,12) npar
      write(12,*) "frame", iframe

      !write(12,15) 'O', 0.0, 0.0, 0.0
      do i= 1, npar
      write(12,15) atsymbol(i), x(i),y(i),z(i)
      end do

      end if

* writing an additional atom

c     write(12,16) iframe, 3.0, 0.000, 0.000

   12 format (1x,i5)
   13 format (1x,a5,2x,i5)
   15 format (1x,a5,3(2x,f10.3))
   16 format (1x,i5,"X", 3(1x,f8.4))

      return
      end
* **************************************
* ******************************************
* get the center of mass from atom positions
* ******************************************
      subroutine massCenter (nk1,nk2,atmass,x,y,z,r1,r2,r3)
      implicit double precision (a-h,o-z) 
      include 'parametros'

      dimension atmass(nmax)
      dimension x(nmax),y(nmax),z(nmax)

* initialize summations

      sumx= 0.d0
      sumy= 0.d0
      sumz= 0.d0
      sumass= 0.d0

* get the mass center of the first nk1 particles

      do i= 1, nk1
       sumass= sumass + atmass(i)
       sumx= sumx + atmass(i) *x(i)
       sumy= sumy + atmass(i) *y(i)
       sumz= sumz + atmass(i) *z(i)
      end do

      r1= sumx/ sumass
      r2= sumy/ sumass
      r3= sumz/ sumass

* print the mass center position

      write(*,20) r1,r2,r3

   20 format(1x,"position of the CM",3(1x,f8.4))

      return
      end
* *************************************************
* **************************************
* collect information for the dynamics
* **************************************
      subroutine molInfo (npar,ncage, atsymbol,atomass, x,y,z,
     &   dtIni,dt,Time, rfac,vfac, tempEquil, sK,xi, radiusf,ifreqShr,
     &   ifreqStor,irestDyn,iwavFunct,nsteps,ibegin,trest,
     &   ngrp1,ngrp2,ngrp3, itest,imethod)

      implicit double precision (a-h,o-z)
      include 'parametros'

      character atsymbol(nmax)*5, atname*12

      dimension atomass(nmax)
      dimension x(nmax), y(nmax), z(nmax)
      dimension sK(nmax)

      character name, titleName*100
      integer atnum

      pi= 3.1415926535898
      speedC= 299792458.d0  ! m/s

* open the file with parameters for the dynamics

      open (unit=13, file= 'olla.inp', status= 'unknown')

* ------------------------------
* read parameters of the dynamics
* ------------------------------

      read(13,*) name, itest     ! make a test of the run
      read(13,*) name, imethod   ! chose the ab-initio program
      read(13,*) titleName       ! title of the dynamics
      read(13,*) name, dt        ! time step of the dynamics
      read(13,*) name, Time      ! total time of the dynamics
      read(13,*) name, rfac      ! scale factor of atom positions
      read(13,*) name, vfac      ! scale factor of atom velocities
      read(13,*) name, npar      ! num of particles
      read(13,*) name, ncage     ! number of cage atoms  [1,ncage]
      read(13,*) name, ngrp1,ngrp2,ngrp3  ! charge for group of atoms
      read(13,*) name, tempEquil ! required equilibrium temperature
      read(13,*) name, springK   ! spring constant K
      read(13,*) name, xi        ! friction constant of the solvent
      read(13,*) name, radiusf   ! factor to shrink the cage
      read(13,*) name, ifreqShr  ! freq to shrink the cage
      read(13,*) name, ifreqStor ! freq to store data in case of a crash
      read(13,*) name, irestDyn  ! restart the dynamics
      read(13,*) name, iwavFunct ! generate or reuse the wavefunction
      read(13,*) name

      write(*,*) "-----------------------------"

      if (irestDyn.eq.1) write(*,14) "dynamics starts for the 1st time"
      if (irestDyn.eq.2) write(*,14) "restarting the dynamics"

      if (itest.eq.1)   write(*,14) "running a test"

      if (itest.eq.2) then ! <---
                                 !
      write(*,14) "running ab-initio computation"
                                 !
      if (imethod.eq.1) write(*,14) "using TeraChem"
      if (imethod.eq.2) write(*,14) "using NWChem"
                                 !
      if (iwavFunct.eq.1) write(*,14)
     &  "creating the wave function from scratch"
      if (iwavFunct.eq.2) write(*,14) "reusing the wave function"
                                 !
      if ((iwavFunct.lt.1).or.(iwavFunct.gt.2)) then
      print*                     !
      print*, "BAD value for starting the wavefunction"
      print*                     !
      stop                       !
      end if                     !
                                 !
      end if ! <-----------------

      if (ncage.gt.npar) then
      print*
      print*, "N cage atoms > N conf atoms"
      print*, "!!!! be a smart guy !!!!"
      print*
      stop
      end if

      write(*,*)
      write(*,*)  titleName
      write(*,*)  "initial parameters of the dynamics"
      write(*,16) "time step [fs]  ", dt
      write(*,16) "total simulation time [fs]", Time
      write(*,16) "rfac      ", rfac
      write(*,16) "vfac      ", vfac
      write(*,14) "number of atoms", npar
      write(*,14) "ncage     ",      ncage
      write(*,18) "charge for group of atoms", ngrp1,ngrp2,ngrp3
      write(*,16) "tempEquil [K]",   tempEquil
      write(*,16) "spring constant [N/cm]",    springK
      write(*,16) "friction constant [ps^-1]", xi
      write(*,16) "factor to shrink the cage", radiusf
      write(*,14) "freq to shrink the cage",   ifreqShr
      write(*,14) "freq to store data",        ifreqStor

      ibegin= 0    ! default value to start the dynamics
      trest= 0.d0  ! default value to start the time
      dtIni= dt    ! store initial time step

      print*, "-----------------------------"
      print*, "reading:"
      print*, "coords [Angst], speeds [Angst/fs] & masses [amu]"

c     print*
c     print*, "num symbol  atnum  atname        atmass        x
c    &y         z"

* ------------------------------
* identify atoms & their features
* ------------------------------

      isum= 0

      do i= 1, npar ! <------
                             !
      isum= isum +1          ! define defaults
                             !
      atsymbol(i)= "NADA"    ! must NOT be an atom symbol
      atnum= 500000          ! must NOT be an atom number
      atname= "NOSE"         ! must NOT be an atom name
      atmass= 500000.0d0     ! must NOT be an atom mass
                             !
      read(13,*) atsymbol(i), x(i), y(i), z(i)
                             !
      call periodicTable (i,atsymbol,atnum,atname,atmass)
                             !
      atomass(i)= atmass     !
                             !
c     write(*,13) i, atsymbol(i), atnum, atname, atmass,
c    &    x(i), y(i), z(i)   !
                             !
      end do ! <-------------

      close (13)   ! file olla.inp

      write(*,14) "number of cage atoms", ncage
      write(*,14) "number of confined atoms", npar-ncage
      write(*,14) "total num of atoms", isum

* -------------------------------
* get initial radius of the cage
* -------------------------------

      big= -1.d0
      ibig= 0
      sum= 0.d0

      do i= 1, ncage ! <-----
                             !
        dist= sqrt (x(i)**2 + y(i)**2 + z(i)**2)
        sum= sum + dist      !
                             !
       if (dist.gt.big) then ! get the largest
        big= dist            ! distance deviation
        ibig= i              !
       end if                !
                             ! all spring constants
        sK(i)= springK       ! are the same
                             !
      end do ! <-------------

      cradius= sum/ dble(ncage) ! central radius

* reduction of the cage after a num of steps

      xnum= int(Time/ real(ifreqShr))
      rfin= cradius *radiusf**xnum

      print*
      write(*,16) "average cage radius [Angst]", cradius
      write(*,16) "largest distance deviation", big
      write(*,14) "of cage atom", ibig

      if (xnum.lt.0.0001) then
       write(*,16) "the cage shall not be reduced"
      else
       write(*,16) "the cage shall be reduced down to [Angst]", rfin
       write(*,16) "after a number of reductions", xnum
      end if

* ------------------------------
* vibration freq of cage atoms
* ------------------------------

* 1 fs= 1.x10^{-15} s
* sqrt(k/m)= sqrt[(N/cm)/amu]= facFreq * s^{-1}

      facFreq= 1./ sqrt(1.6605402d-29)

* omega= 2pi*freq= sqrt(k/m)  ==> freq= sqrt(k/m)/ 2pi

      freq1= facFreq * sqrt(springK/atomass(2)) /(2.*pi)  ! in 1/s
      freq2= freq1/ (100.*speedC)    ! in cm^-1
      period= 1./(freq1 *1.0d-15)    ! in fs

      print*
      write(*,*)  "vibrations of cage atoms (void space)"
      write(*,17) "vib freq. [1/s] ", freq1
      write(*,16) "vib freq. [1/cm]", freq2
      write(*,16) "vib period [fs]", period
      print*

      v1= Time/dt
      nsteps= int(v1)  ! number of time steps to compute
      print*, "energy calculations to perform", nsteps

* -----------------
* non-used variables
* -----------------

      rfac= 1.d0    ! non-used variables
      vfac= 1.d0    ! are given numbers

      print*
      write(*,*) "non-used variables have values:"
      write(*,16) "rfac", rfac
      write(*,16) "vfac", vfac

      print*, "* -------------------------------"
      print*, "   wall fluctuations in progress"
      print*, "* -------------------------------"

   13 format (1x,i3,1x,a,1x,i3,1x,a,4(1x,f9.4))
   14 format (1x,a,1x,i7)
   16 format (1x,a,1x,f11.3)
   17 format (1x,a,1x,e10.4)
   18 format (1x,a,3i6)

      return
      end
* *******************************
* **************************************************
* given atom symbol or atom number, this module finds
* (atom symbol, atom number, atom name, atom mass)
* using the periodic table of the elements
* **************************************************
      subroutine periodicTable (i,atsymbol,atnum,atname,atmass)
      implicit double precision (a-h,o-z)
      parameter (m1=1500)

      character atsymbol(m1)*5, atname*12
      integer atnum

* atom masses in au
* 1 atomic mass unit (amu) = 1.6605402 x 10^{-27} Kg

* WARNING: be careful in other programs with
* implicit real*8 (a-h,o-z)

      if ((atsymbol(i).eq.'H').or.(atnum.eq.1)) then
      atnum= 1
      atname= "Hydrogen"
      atsymbol(i)= 'H'
      atmass= 1.00794
      end if

      if ((atsymbol(i).eq.'He').or.(atnum.eq.2)) then
      atnum= 2
      atname= "Helium"
      atsymbol(i)= 'He'
      atmass= 4.002602
      end if

      if ((atsymbol(i).eq.'Li').or.(atnum.eq.3)) then
      atnum= 3
      atname= "Lithium"
      atsymbol(i)= 'Li'
      atmass= 6.941
      end if

      if ((atsymbol(i).eq.'Be'.or.(atnum.eq.4))) then
      atnum= 4
      atname= "Beryllium"
      atsymbol(i)= 'Be'
      atmass= 9.012182
      end if

      if ((atsymbol(i).eq.'B').or.(atnum.eq.5)) then
      atnum= 5
      atname= "Boron"
      atsymbol(i)= 'B'
      atmass= 10.811
      end if

      if ((atsymbol(i).eq.'C').or.(atnum.eq.6)) then
      atnum= 6
      atname= "Carbon"
      atsymbol(i)= 'C'
      atmass= 12.0107
      end if

      if ((atsymbol(i).eq.'N').or.(atnum.eq.7)) then
      atnum= 7
      atname= "Nitrogen"
      atsymbol(i)= 'N'
      atmass= 14.003070
      end if

      if ((atsymbol(i).eq.'O').or.(atnum.eq.8)) then
      atnum= 8
      atname= "Oxygen"
      atsymbol(i)= 'O'
      atmass= 15.99491
      end if

      if ((atsymbol(i).eq.'F').or.(atnum.eq.9)) then
      atnum= 9
      atname= "Fluorine"
      atsymbol(i)= 'F'
      atmass= 18.9984032
      end if

      if ((atsymbol(i).eq.'Ne').or.(atnum.eq.10)) then
      atnum= 10
      atname= "Neon"
      atsymbol(i)= 'Ne'
      atmass= 20.1797
      end if

      if ((atsymbol(i).eq.'Na').or.(atnum.eq.11)) then
      atnum= 11
      atname= "Sodium"
      atsymbol(i)= 'Na'
      atmass= 22.989770
      end if

      if ((atsymbol(i).eq.'Mg').or.(atnum.eq.12)) then
      atnum= 12
      atname= "Magnesium"
      atsymbol(i)= 'Mg'
      atmass= 23.985040
      end if

      if ((atsymbol(i).eq.'Al').or.(atnum.eq.13)) then
      atnum= 13
      atname= "Aluminum"
      atsymbol(i)= 'Al'
      atmass= 26.981540
      end if

      if ((atsymbol(i).eq.'Si').or.(atnum.eq.14)) then
      atnum= 14
      atname= "Silicon"
      atsymbol(i)= 'Si'
      atmass= 27.976930
      end if

      if ((atsymbol(i).eq.'P').or.(atnum.eq.15)) then
      atnum= 15
      atname= "Phosphorous"
      atsymbol(i)= 'P'
      atmass= 30.973760
      end if

      if ((atsymbol(i).eq.'S').or.(atnum.eq.16)) then
      atnum= 16
      atname= "Sulfur"
      atsymbol(i)= 'S'
      atmass= 32.066
      end if

      if ((atsymbol(i).eq.'Cl').or.(atnum.eq.17)) then
      atnum= 17
      atname= "Chlorine"
      atsymbol(i)= 'Cl'
      atmass= 35.4527
      end if

      if ((atsymbol(i).eq.'Ar').or.(atnum.eq.18)) then
      atnum= 18
      atname= "Argon"
      atsymbol(i)= 'Ar'
      atmass= 39.948
      end if

      if ((atsymbol(i).eq.'K').or.(atnum.eq.19)) then
      atnum= 19
      atname= "Potassium"
      atsymbol(i)= 'K'
      atmass= 39.0983
      end if

      if ((atsymbol(i).eq.'Ca').or.(atnum.eq.20)) then
      atnum= 20
      atname= "Calcium"
      atsymbol(i)= 'Ca'
      atmass= 40.078
      end if

      if ((atsymbol(i).eq.'Sc').or.(atnum.eq.21)) then
      atnum= 21
      atname= "Scandium"
      atsymbol(i)= 'Sc'
      atmass= 44.955910
      end if

      if ((atsymbol(i).eq.'Ti').or.(atnum.eq.22)) then
      atnum= 22
      atname= "Titanium"
      atsymbol(i)= 'Ti'
      atmass= 47.867
      end if

      if ((atsymbol(i).eq.'V').or.(atnum.eq.23)) then
      atnum= 23
      atname= "Vanadium"
      atsymbol(i)= 'V'
      atmass= 50.9415
      end if

      if ((atsymbol(i).eq.'Cr').or.(atnum.eq.24)) then
      atnum= 24
      atname= "Chromium"
      atsymbol(i)= 'Cr'
      atmass= 51.9961
      end if

      if ((atsymbol(i).eq.'Mn').or.(atnum.eq.25)) then
      atnum= 25
      atname= "Manganese"
      atsymbol(i)= 'Mn'
      atmass= 54.938049
      end if

      if ((atsymbol(i).eq.'Fe').or.(atnum.eq.26)) then
      atname= "Iron"
      atsymbol(i)= 'Fe'
      atmass= 55.845
      end if

      if ((atsymbol(i).eq.'Co').or.(atnum.eq.27)) then
      atnum= 27
      atname= "Cobalt"
      atsymbol(i)= 'Co'
      atmass= 58.933200
      end if

      if ((atsymbol(i).eq.'Ni').or.(atnum.eq.28)) then
      atnum= 28
      atname= "Nickel"
      atsymbol(i)= 'Ni'
      atmass= 58.6934
      end if

      if ((atsymbol(i).eq.'Cu').or.(atnum.eq.29)) then
      atnum= 29
      atname= "Copper"
      atsymbol(i)= 'Cu'
      atmass= 63.546
      end if

      if ((atsymbol(i).eq.'Zn').or.(atnum.eq.30)) then
      atnum= 30
      atname= "Zinc"
      atsymbol(i)= 'Zn'
      atmass= 65.39
      end if

      if ((atsymbol(i).eq.'Ga').or.(atnum.eq.31)) then
      atnum= 31
      atname= "Gallium"
      atsymbol(i)= 'Ga'
      atmass= 69.723
      end if

      if ((atsymbol(i).eq.'Ge').or.(atnum.eq.32)) then
      atnum= 32
      atname= "Germanium"
      atsymbol(i)= 'Ge'
      atmass= 72.61
      end if

      if ((atsymbol(i).eq.'As').or.(atnum.eq.33)) then
      atnum= 33
      atname= "Arsenic"
      atsymbol(i)= 'As'
      atmass= 74.92160
      end if

      if ((atsymbol(i).eq.'Se').or.(atnum.eq.34)) then
      atnum= 34
      atname= "Selenium"
      atsymbol(i)= 'Se'
      atmass= 78.96
      end if

      if ((atsymbol(i).eq.'Br').or.(atnum.eq.35)) then
      atnum= 35
      atname= "Bromine"
      atsymbol(i)= 'Br'
      atmass= 79.904
      end if

      if ((atsymbol(i).eq.'Kr').or.(atnum.eq.36)) then
      atnum= 36
      atname= "Kripton"
      atsymbol(i)= 'Kr'
      atmass= 83.80
      end if

      if ((atsymbol(i).eq.'Rb').or.(atnum.eq.37)) then
      atnum= 37
      atname= "Rubidium"
      atsymbol(i)= 'Rb'
      atmass= 85.4678
      end if

      if ((atsymbol(i).eq.'Sr').or.(atnum.eq.38)) then
      atnum= 38
      atname= "Strontium"
      atsymbol(i)= 'Sr'
      atmass= 87.62
      end if

      if ((atsymbol(i).eq.'Y').or.(atnum.eq.39)) then
      atnum= 39
      atname= "Yttrium"
      atsymbol(i)= 'Y'
      atmass= 88.90585
      end if

      if ((atsymbol(i).eq.'Zr').or.(atnum.eq.40)) then
      atnum= 40
      atname= "Zirconium"
      atsymbol(i)= 'Zr'
      atmass= 91.224
      end if

      if ((atsymbol(i).eq.'Nb').or.(atnum.eq.41)) then
      atnum= 41
      atname= "Niobium"
      atsymbol(i)= 'Nb'
      atmass= 92.90638
      end if

      if ((atsymbol(i).eq.'Mo').or.(atnum.eq.42)) then
      atnum= 42
      atname= "Molybdenum"
      atsymbol(i)= 'Mo'
      atmass= 95.94
      end if

      if ((atsymbol(i).eq.'Tc').or.(atnum.eq.43)) then
      atnum= 43
      atname= "Technetium"
      atsymbol(i)= 'Tc'
      atmass= 98
      end if

      if ((atsymbol(i).eq.'Ru').or.(atnum.eq.44)) then
      atnum= 44
      atname= "Ruthenium"
      atsymbol(i)= 'Ru'
      atmass= 101.07
      end if

      if ((atsymbol(i).eq.'Rh').or.(atnum.eq.45)) then
      atnum= 45
      atname= "Rhodium"
      atsymbol(i)= 'Rh'
      atmass= 102.90550
      end if

      if ((atsymbol(i).eq.'Pd').or.(atnum.eq.46)) then
      atnum= 46
      atname= "Palladium"
      atsymbol(i)= 'Pd'
      atmass= 106.42
      end if

      if ((atsymbol(i).eq.'Ag').or.(atnum.eq.47)) then
      atnum= 47
      atname= "Silver"
      atsymbol(i)= 'Ag'
      atmass= 107.8682
      end if

      if ((atsymbol(i).eq.'Cd').or.(atnum.eq.48)) then
      atnum= 48
      atname= "Cadmium"
      atsymbol(i)= 'Cd'
      atmass= 112.411
      end if

      if ((atsymbol(i).eq.'In').or.(atnum.eq.49)) then
      atnum= 49
      atname= "Indium"
      atsymbol(i)= 'In'
      atmass= 114.818
      end if

      if ((atsymbol(i).eq.'Sn').or.(atnum.eq.50)) then
      atnum= 50
      atname= "Tin"
      atsymbol(i)= 'Sn'
      atmass= 118.710
      end if

      if ((atsymbol(i).eq.'Sb').or.(atnum.eq.51)) then
      atnum= 51
      atname= "Antimony"
      atsymbol(i)= 'Sb'
      atmass= 121.76
      end if

      if ((atsymbol(i).eq.'Te').or.(atnum.eq.52)) then
      atnum= 52
      atname= "Tellurium"
      atsymbol(i)= 'Te'
      atmass= 127.60
      end if

      if ((atsymbol(i).eq.'I').or.(atnum.eq.53)) then
      atnum= 53
      atname= "Iodine"
      atsymbol(i)= 'I'
      atmass= 126.90447
      end if

      if ((atsymbol(i).eq.'Xe').or.(atnum.eq.54)) then
      atnum= 54
      atname= "Xenon"
      atsymbol(i)= 'Xe'
      atmass= 131.29
      end if

      return
      end
* ************************************************
* *************************
* perturb the confined atoms
* *************************
      subroutine perturb (npar,nstatic, x,y,z, vx,vy,vz,
     &    xnew,ynew,znew, vxnew,vynew,vznew, rfac,vfac)

      implicit double precision (a-h,o-z)
      include 'parametros'

      dimension  x(nmax),  y(nmax),  z(nmax)
      dimension vx(nmax), vy(nmax), vz(nmax)

      dimension  xnew(nmax),  ynew(nmax),  znew(nmax)
      dimension vxnew(nmax), vynew(nmax), vznew(nmax)

      do i= 1, npar ! <-------
                              !
      if (i.le.nstatic) then  !
       xnew(i)= x(i)          !
       ynew(i)= y(i)          ! in Angst
       znew(i)= z(i)          !
                              !
       vxnew(i)= vx(i)        !
       vynew(i)= vy(i)        ! in Angst/fs
       vznew(i)= vz(i)        !
      endif                   !
                              !
                              !
      if (i.gt.nstatic) then  !
       xnew(i)= x(i) *rfac    !
       ynew(i)= y(i) *rfac    ! in Angst
       znew(i)= z(i) *rfac    !
                              !
       vxnew(i)= vx(i) *vfac  !
       vynew(i)= vy(i) *vfac  ! in Angst/fs
       vznew(i)= vz(i) *vfac  !
      endif                   !
                              !
      end do ! <--------------

      return
      end
* ****************************
* **********************************
* get polar coordinates (rho,teta,fi)
* from the cartesian coordinates (a,b,c)
* **********************************
      subroutine polar (a,b,c,rho,teta,fi)
      implicit double precision (a-h,o-z)

      pi= 3.1415926535898
      rad= pi/180.0

*  0 <= rho
*  0 <= teta <= 2pi    teta is in the (x,y) plane
*  0 <= fi <= pi       goes from +z to -z

      rho= sqrt(a**2 + b**2 + c**2)
      fi= acos(c/rho)

* get teta in the 1st quadrant

      if ((a.ge.0.).and.(b.eq.0.)) then
      teta= 0.d0
      end if

      if ((a.gt.0.).and.(b.gt.0.)) then
      teta= atan(b/a)
      end if

      if ((a.eq.0.).and.(b.gt.0.)) then
      teta=  pi/2.d0
      end if

* get teta in the 2nd quadrant

      if ((a.lt.0.).and.(b.gt.0.)) then
      teta= pi - atan(abs(b/a))
      end if

      if ((a.lt.0.).and.(b.eq.0.)) then
      teta= pi
      end if

* get teta in the 3rd quadrant

      if ((a.lt.0.).and.(b.lt.0.)) then
      teta= pi + atan(abs(b/a))
      end if

      if ((a.eq.0.).and.(b.lt.0.)) then
      teta= 3./2.* pi
      end if

* get teta in the 4th quadrant

      if ((a.gt.0.).and.(b.lt.0.)) then
      teta= 2.*pi - atan(abs(b/a))
      end if

* the point has no teta angle

      if ((a.eq.0.).and.(b.eq.0.)) then
      print*, "ERROR: x=0 and y=0"
      end if

      return
      end
* **********************************
* ***********************************
* store data of the dynamics for a restart
*
* WARNING:
* all arrays & summations in memory are lost after a crash,
* only information written to output files is preserved
* make a backup of all files before making a restart,
* ***********************************
      subroutine restartOlla (ilabel,ibegin,dt,time,npar,ncage,x,y,z,
     &   vx,vy,vz,Gradx,Grady,Gradz, x0,y0,z0,sK, refPE,
     &    ibTime1,ibTime2,ibTime3)

      implicit double precision (a-h,o-z)
      include 'parametros'

      character*25 name

      dimension x0(nmax), y0(nmax), z0(nmax), sK(nmax)
      dimension  x(nmax),  y(nmax),  z(nmax)
      dimension vx(nmax), vy(nmax), vz(nmax)
      dimension Gradx(nmax), Grady(nmax), Gradz(nmax)

* ilabel= 1 ==> write information
* ilabel= 2 ==> read  information

* ------------------------
* open the file for storing data
* ------------------------

      open (unit=15, file= 'restart.dat',status= 'unknown')

* ------------------------
* write parameters & variables
* ------------------------

      if (ilabel.eq.1) then ! <---
                                  !
      write(15,"(a,i6)") "step", ibegin
      write(15,"(a,f9.3)") "time", time
      write(15,"(a,f9.3)") "time-step", dt
      write(15,*) "referencePE", refPE
      write(15,"(a,3i6)") "rnd-times", ibTime1,ibTime2,ibTime3
                                  !
      write(15,*) "positions"     !
                                  !
      do i= 1, npar               !
      write(15,12) i, x(i), y(i), z(i)
      end do                      !
                                  !
      write(15,*) "speeds"        !
                                  !
      do i= 1, npar               !
      write(15,12) i, vx(i), vy(i), vz(i)
      end do                      !
                                  !
      write(15,*) "gradients"     !
                                  !
      do i= 1, npar               !
      write(15,12) i, Gradx(i), Grady(i), Gradz(i)
      end do                      !
                                  !
      write(15,*) "cage-reference-positions"
                                  !
      do i= 1, ncage              !
      write(15,12) i, x0(i), y0(i), z0(i)
      end do                      !
                                  !
      write(15,*) "spring-constants"
                                  !
      do i= 1, ncage              !
      write(15,12) i, sK(i)       !
      end do                      !
                                  !
      end if ! <------------------

* ------------------------
* read parameters & variables
* ------------------------

      if (ilabel.eq.2) then ! <---
                                  !
      write(*,"(a)") "reading stored variables"
                                  !
      read(15,*) name, ibegin     !
      read(15,*) name, time       !
      read(15,*) name, dt         !
      read(15,*) name, refPE      !
      read(15,*) name, ibTime1,ibTime2,ibTime3
                                  !
      read(15,*) name             !
      print*, name                !
                                  !
      do i= 1, npar               !
      read(15,*) k, x(i), y(i), z(i)
      end do                      !
                                  !
      read(15,*) name             !
      print*, name                !
                                  !
      do i= 1, npar               !
      read(15,*) k, vx(i), vy(i), vz(i)
      end do                      !
                                  !
      read(15,*) name             !
      print*, name                !
                                  !
      do i= 1, npar               !
      read(15,*) k, Gradx(i), Grady(i), Gradz(i)
      end do                      !
                                  !
      read(15,*) name             !
      print*, name                !
                                  !
      do i= 1, ncage              !
      read(15,*) k, x0(i), y0(i), z0(i)
      end do                      !
                                  !
      read(15,*) name             !
      print*, name                !
                                  !
      do i= 1, ncage              !
      read(15,*) k, sK(i)         !
      end do                      !
                                  !
      write(*,"(1x,a,3i5)") "rnd times to break bonds",
     &  ibTime1,ibTime2,ibTime3   !
                                  !
      end if ! <------------------

      close (15)

   12 format(i6,3f10.5)

      return
      end
* *******************************
* ******************
* chose the wavefunction (WF)
* sequence of events in breaking a bond
*
* time intervals         events
* T0       t<=T1   T1    T1+1<=t<=T1+dT    T1+dT      T1+dT+1<=t<=T2
* normal   reuse  break    excited        down to         normal
* state     WF    bond      state      normal state       state
* WF=1      2      4          5              1              2
*
*                  T2    T2+1<=t<=T2+dT    T2+dT      T2+dT+1<=t<=T3
*                 break    excited        down to         normal
*                 bond      state      normal state       state
*         WF=      4          5               1              2
*
*                  T3    T3+1<=t<=T3+dT    T3+dT        T2+dT+1<=t
*                 break    excited        down to         normal
*                 bond      state      normal state       state
*         WF=      4          5               1              2
* ************************************
      subroutine setWFs (kk,imethod,ibTime1,ibTime2,ibTime3,iAft,
     &   npar,atsymbol,x,y,z)
      implicit double precision (a-h,o-z)
      include 'parametros'

      character atsymbol(nmax)*5
      dimension x(nmax),y(nmax),z(nmax)


* normal state <= T1

      if (kk.le.ibTime1) then
      if (imethod.eq.1) call inptTera (2)
      if (imethod.eq.2) call inptNwchem (2,npar,atsymbol,x,y,z)
      end if

* T1 <= excited state <= T1+dT

      if ((kk.ge.ibTime1+1).and.(kk.le.ibTime1+iAft)) then
      if (imethod.eq.1) call inptTera (5)
      if (imethod.eq.2) call inptNwchem (5,npar,atsymbol,x,y,z)
      end if

* T1+dT <= down to normal state <= T2

      if ((kk.ge.ibTime1+iAft+1).and.(kk.le.ibTime2)) then
      if (imethod.eq.1) call inptTera (2)
      if (imethod.eq.2) call inptNwchem (2,npar,atsymbol,x,y,z)
      end if

* T2 <= excited state <= T2+dT

      if ((kk.ge.ibTime2+1).and.(kk.le.ibTime2+iAft)) then
      if (imethod.eq.1) call inptTera (5)
      if (imethod.eq.2) call inptNwchem (5,npar,atsymbol,x,y,z)
      end if

* T2+dT <= down to normal state <= T3

      if ((kk.ge.ibTime2+iAft+1).and.(kk.le.ibTime3)) then
      if (imethod.eq.1) call inptTera (2)
      if (imethod.eq.2) call inptNwchem (2,npar,atsymbol,x,y,z)
      end if

* T3 <= excited state <= T3+dT

      if ((kk.ge.ibTime3+1).and.(kk.le.ibTime3+iAft)) then
      if (imethod.eq.1) call inptTera (5)
      if (imethod.eq.2) call inptNwchem (5,npar,atsymbol,x,y,z)
      end if

* T3+dT <= down to normal state

      if  (kk.ge.ibTime3+iAft+1) then
      if (imethod.eq.1) call inptTera (2)
      if (imethod.eq.2) call inptNwchem (2,npar,atsymbol,x,y,z)
      end if

      if (imethod.eq.1) call writeMol (npar,atsymbol,x,y,z)

      return
      end
* ******************
* ********************
* module to write the ab-initio script
* the file created here shall be called from
* the main program to run the ab-initio program
* ********************
      subroutine scriptAB (imethod)
      implicit double precision (a-h,o-z)
      include 'parametros'

      character name1*44, name2*38, name3*32, name4*33

* char(34)= "
* char(35)= #

* large paths for nwchem

      name1= "/usr/local/openmpi-1.6.4/bin/mpirun -np 1 "
      name2= "/home/rso/nwchem/bin/LINUX64/nwchem "
      name3= "/home/rso/work/nwchem.inp "
      name4= "> /home/rso/work/nwchem.out"

* -------------------
* script for TeraChem
* -------------------

      if (imethod.eq.1) then ! <-----------
                                           !
      print*, "creating the script to run terachem"
      open(14,file="script.tera",status="unknown")
                                           !
      write(14,*)                          !
      write(14,*) char(35), " script to submit terachem"
      write(14,*) char(35), " commands may be changed in module"
      write(14,*) char(35), " scriptTera"  !
      write(14,*)                          !
c     write(14,*) "  echo  "               !
      write(14,*) "  echo  ", char(34),    !
     &      "!!!!!! terachem is working now !!!!!!", char(34)
      write(14,*)                          !
      write(14,*) char(35), " submit the terachem job"
                                           !
 1000 write(14,*) " /usr/local/TeraChem/terachem",
     &              "  terachem.inp  >  terachem.out"
                                           !
      write(14,*)                          !
      write(14,*) char(35), " get the dft ene"
      write(14,*) "  egrep  ", char(34),   !
     &   "FINAL ENERGY|CENTER OF|DIPOLE", char(34), "   terachem.out"
      write(14,*)                          !
      write(14,*) "  echo  ", char(34),    !
     &      "!!!!!! terachem finished this step !!!!!!", char(34)
      write(14,*) "  echo  ", char(34), "   ", char(34)
      write(14,*) "  tail -4 terachem.out" !
      write(14,*)                          !
                                           !
      call system ("chmod u+x script.tera")! make it executable
                                           !
      end if ! <---------------------------

* ----------------
* script for NWChem
* ----------------

      if (imethod.eq.2) then ! <-------------
                                             !
      print*, "creating the script to run nwchem"
      open(14,file="script.nwchem",status="unknown")
                                             !
      write(14,*)                            !
      write(14,*) char(35), " script to submit nwchem"
      write(14,*) char(35), " commands may be changed in module"
      write(14,*) char(35), " scriptNwchem"  !
      write(14,*)                            !
c     write(14,*) "  echo  "                 !
      write(14,*) "  echo  ", char(34),      !
     &      "!!!!!! nwchem is working now !!!!!!", char(34)
      write(14,*)                            !
      write(14,*) char(35), " submit the nwchem job"
                                             !
      write(14,*)        name2, name3, name4 !
c     write(14,*) name1, name2, name3, name4 !
c     write(14,*) " mpirun -np 10 /bin/nwchem  nwchem.inp > nwchem.out"
      write(14,*)                            !
                                             !
      write(14,*) char(35), " get the dft ene"
      write(14,*) "  egrep  ", char(34), "DFT energy", char(34),
     &   " nwchem.out"                       !
      write(14,*)                            !
                                             !
      write(14,*) "  echo  ", char(34),      !
     &      "!!!!!! nwchem finished this step !!!!!!", char(34)
      write(14,*)  "  tail -1  ", " nwchem.out"
      write(14,*) "  echo  "                 !
      write(14,*)                            !
                                             !
      call system ("chmod u+x script.nwchem")! make it executable
                                             !
      end if ! <-----------------------------

      close (14)

      return
      end
* ********************
* ***********************************
* shrink the cage and atoms inside the cage
* ***********************************
      subroutine shrink (npar,ncage,x0,y0,z0,x,y,z,radiusf)
      implicit double precision (a-h,o-z)
      include 'parametros'

      dimension  x(nmax),  y(nmax),  z(nmax)
      dimension x0(nmax), y0(nmax), z0(nmax)

* --------------------------------------------
* resize the reference radius of cage particles
* --------------------------------------------

      sum1= 0.d0
      sum2= 0.d0

      do i= 1, ncage ! <----
                            !
      a= x0(i)              !
      b= y0(i)              !
      c= z0(i)              !
                            !
      dold= sqrt (a**2 +b**2 +c**2)
                            !
                            ! transform to polar coords
      call polar (a,b,c,rho,teta,fi)
                            !
      rho= rho *radiusf     ! scale radial distances
                            !
                            ! back to cartesian coords
      call cartesian (rho,teta,fi,a,b,c)
                            !
      x0(i)= a              ! store reference
      y0(i)= b              ! positions
      z0(i)= c              !
                            !
      x(i)= a               ! present positions
      y(i)= b               !
      z(i)= c               !
                            !
      dnew= sqrt (a**2 +b**2 +c**2)
                            !
      sum1= sum1 + dold     ! sum distances
      sum2= sum2 + dnew     ! of cage atoms
                            !
      end do ! <------------

* average radius

      rold= sum1 /dble(ncage)
      rnew= sum2 /dble(ncage)

      write(*,"(1x,a,1x,f9.4)") "old reference rad [Angst]", rold
      write(*,"(1x,a,1x,f9.4)") "new reference rad [Angst]", rnew
      write(*,"(1x,a,1x,f9.4)") "shrinking factor", radiusf

* --------------------------------------------
* resize the present radii of confined particles
* --------------------------------------------

      do i= ncage+1, npar ! <--
                               !
      a= x(i)                  ! present positions
      b= y(i)                  !
      c= z(i)                  !
                               !
      call polar (a,b,c,rho,teta,fi)
                               !
      rho= rho *radiusf        !
                               !
      call cartesian (rho,teta,fi,a,b,c)
                               !
      x(i)= a                  ! scaled positions
      y(i)= b                  !
      z(i)= c                  !
                               !
      end do ! <---------------

      return
      end
* ***********************************
* ***********************************
* slow down speeds of confined atoms close to the cage
* to simulate stochastic loss of kinetic energy, thus
* simulating a cooling process
* ***********************************
      subroutine slowSpeeds (time,npar,ncage,x,y,z,vx,vy,vz,
     &     heatup,cooldown)
      implicit double precision (a-h,o-z)
      include 'parametros'

      dimension  x(nmax),  y(nmax),  z(nmax)
      dimension vx(nmax), vy(nmax), vz(nmax)

      character*3 heatup, cooldown

      integer ilabel(nmax)

* --------------------
* initialize parameters
* --------------------

      do i= 1, npar ! labels for atoms whose
      ilabel(i)= 0  ! speeds shall be reduced
      end do

      xlambda= 0.93 ! factor to reduce atom speeds
      dthresh= 3.50 ! threshold distance to the cage

      if (time.eq.1.) then
      write(*,12) "factor to reduce atom speeds", xlambda
      write(*,12) "distance from cage for cooling down", dthresh
      end if

   12 format (1x,a,2x,f6.3)

c     print*
c     print*, "atoms with scaled speeds at time=", t
c     print*, "atom number         old speed         new speed"

* --------------------
* classify atoms
* --------------------

      icont= 0

      do j= 1, ncage ! <-----  atoms of the cage
      do i= ncage+1, npar    ! atoms in the cage
                             !
       dx= (x(i)-x(j))**2    !
       dy= (y(i)-y(j))**2    !
       dz= (z(i)-z(j))**2    !
       dist= sqrt(dx+dy+dz)  ! find atoms close
                             ! to the cage
                             !
      if ((dist.lt.dthresh).and.(ilabel(i).eq.0)) then
       icont= icont +1       !
       Vold= sqrt(vx(i)**2 + vy(i)**2 + vz(i)**2)
                             !
       vx(i)= vx(i) *xlambda !
       vy(i)= vy(i) *xlambda ! scale down
       vz(i)= vz(i) *xlambda ! these speeds
                             !
       ilabel(i)= i          ! tag this atom
       Vnew= xlambda* Vold   ! new speed
                             !
c     write(*,14) i, Vold, Vnew
      end if                 !
                             !
      end do                 !
      end do ! <-------------

      print*, "num of cooled atoms", icont

* reset switches

      heatup=   "off"
      cooldown= "off"

   14 format (1x,i5,2(1x,f9.5))

      return
      end
* *****************************
* ***********************************
* evolve positions & speeds using the Taylor expansion
* ***********************************
      subroutine taylor (npar, nstatic, atomass,
     &         x,y,z, vx,vy,vz, Gradx,Grady,Gradz,
     &         xnew,ynew,znew, vxnew,vynew,vznew, dt)

      implicit double precision (a-h,o-z)
      include 'parametros'

      dimension  atomass(nmax)
      dimension  x(nmax),  y(nmax),  z(nmax)
      dimension vx(nmax), vy(nmax), vz(nmax)
      dimension Gradx(nmax), Grady(nmax), Gradz(nmax)

      dimension  xnew(nmax),  ynew(nmax),  znew(nmax)
      dimension vxnew(nmax), vynew(nmax), vznew(nmax)

      print*, "evolving coords and speeds with Taylor"

* nstatic son las particulas "estaticas" del contenedor 

* -------------------------
*  accel in terms of dE/dx
* -------------------------

* Fx= -dE/dx

*  Xnew= x + vx*dt +ax*dt^2/2
*      = x + vx*dt -(1/m)*Gradx*dt^2/2

* Vxnew= vx +ax*dt
*      = vx -(1/m)*Gradx*dt

* factor to change units

      faccl= 0.49616184

* ----------------------------------------
* progress positions & speeds of confined atoms
* ----------------------------------------

      do i= nstatic+1, npar ! <--
                                 ! in Angst/fs^2
      ax= -faccl *Gradx(i) /atomass(i)
      ay= -faccl *Grady(i) /atomass(i)
      az= -faccl *Gradz(i) /atomass(i)
                                 !
      vxnew(i)= vx(i) + ax*dt    ! in Angst/fs
      vynew(i)= vy(i) + ay*dt    !
      vznew(i)= vz(i) + az*dt    !
                                 ! in Angst
      xnew(i)= x(i) + vx(i)*dt + ax *dt**2 /2.d0
      ynew(i)= y(i) + vy(i)*dt + ay *dt**2 /2.d0
      znew(i)= z(i) + vz(i)*dt + az *dt**2 /2.d0
                                 !
      end do ! <-----------------

      return
      end
* *******************************

* **************************
* informative module of
* units transformations
* **************************
      subroutine units

      write (*,*) "-------------------------"
      write (*,*) "physical variables & their unis"
      write (*,*) "-------------------------"
      write (*,*) "[position]= Angst              [mass]= Kg"
      write (*,*) "[speeds]=   Angst/fs           [time]= fs"
      write (*,*) "[acceleration]= Angst/fs^2     [Force]= Newton"
      write (*,*) "                               [Energy]= Joules"
      write (*,*)

      write (*,*) "1 Angst= 1x10^{-10} mt"
      write (*,*) "1 mt=    1x10^{+10} Angst"
      write (*,*) "1 Bohr= 0.529178 Angst= 0.529178 x10^{-10} mt"
      write (*,*) "1 amu=   1.6605402 x10^{-27} Kg"
      write (*,*) "1 Hartree= 4.35988 x10^{-18} J"
      write (*,*) "1 Hartree= 27.211396132 eV"
      write (*,*) "1 fs= 1x10^{-15} s"
      write (*,*) "1 atmosphere= 101325.0 pascals"
      write (*,*) "1 pascal= 1/101325.0 atmosphere"
      write (*,*)

      write (*,*) "Newton= Kg *mt/ s^2"
      write (*,*) "Joule= Newton * mt= Kg *(mt/s)^2"
      write (*,*)

      write (*,*) "kB= 1.3806504 x10^{-23}    in Joule/Kelvin"
      write (*,*) "pi= 3.1415926535898"
      write (*,*) "c= 2.99792458 x10^8 m/s= 2.99792458 x10^{10} cm/s",
     &          " speed of light"
      write (*,*) "mole= 6.02214179 x10^{23}"
      write (*,*) "R= 8.3144621 J/(mol*Kelvin)   universal constant",
     &        " of gases"
      write (*,*) "Planck constant= hbar=6.62606957x10^{-34} (kg*m^2)/s"
      write (*,*)
      write (*,*) "--------------------------------------"
      write (*,*) "units of positon= acceleration * dt^2"
      write (*,*) "--------------------------------------"
      write (*,*) "[atom mass]= atomass x amu"
      write (*,*) "[Force]= [dE/dx]= Hartree/Bohr"
      write (*,*) "[dt]= fs"
      write (*,*)

      write (*,*) "rx= ax * dt^2 /2"
      write (*,*) " = (Fx/m) * dt^2 /2"
      write (*,*) " = -(1/m) (dE/dx) *dt^2 /2"
      write (*,*) " = -(1/m) (Gradx) *dt^2 /2 [(1/amu)*Hartree/Bohr",
     &            " *fs^2]"
      write (*,*)

      write (*,*) "[(1/amu) *Hartree/Bohr *fs^2]"
      write (*,*) "= 1/(1.6605402 x10^{-27} Kg) *"
      write (*,*) "    4.35988 x10^{-18} J/ 0.529178 x10^{-10} mt",
     &             " *1x10^{-30} s^2"
      write (*,*) "= 4.35988 /(1.6605402 * 0.529178)x10^{+27-18+10-30}",
     &            " J/(Kg*mt)*s^2"
      write (*,*) "= 4.9616184 x10^{-11} m"
      write (*,*)

      write (*,*) "transformation factor for positions derived from",
     &           " the force"
      write (*,*)

      write (*,*) "faccl= 0.49616184    dimensionless"
      write (*,*)

      write (*,*) "faccl*Angst= (1/amu)*(Hartree/Bohr)*(fs^2)"
      write (*,*)

      write (*,*) "----------------------------------"
      write (*,*) "units of speed= acceleration * dt"
      write (*,*) "----------------------------------"
      write (*,*) "from the result for positions we have:"
      write (*,*) "fvel= 0.49616184"
      write (*,*) "fvel *(Angst/fs)= (1/amu) *(Hartree/Bohr) *fs"
      write (*,*)

      write (*,*) "------------------------"
      write (*,*) "units of kinetic energy"
      write (*,*) "------------------------"
      write (*,*) "[Ek]= [mass *vx(i)^2]= amu *(Angst/fs)^2"
      write (*,*)

      write (*,*) "1 amu = 1.6605402 x10^{-27} Kg"
      write (*,*) "1 Angst= 1x10^{-10} mt"
      write (*,*) "1 fsec= 1x10^{-15} sec"
      write (*,*)

      write (*,*) "[Ek]= [1.6605402x10^{-27} Kg] [1x10^{-10} mt]^2/",
     &            " [1x10^{-15} sec]^2"
      write (*,*)

      write (*,*) "Joule= Newton*mt= Kg*(mt/s)^2"
      write (*,*)

      write (*,*) "transformation factor between units of kinetic",
     &            " energy"

      write (*,*) "facEkin= 1.6605402d-17"
      write (*,*) "facEkin * Joule= amu *(Angst/fs)^2"
      write (*,*)

      write (*,*) "-----------------------------"
      write (*,*) "units of pressure= force/area"
      write (*,*) "-----------------------------"
      write (*,*) "radial force= F . R/|R|= -(dE/dx,dE/dy,dE/dz) .",
     &            " R/|R| [Hartree/Bohr]"
      write (*,*)

      write (*,*) "[Hartree/Bohr]= 4.35988 x10^{-18} J/",
     &            " 0.529178x10^{-10} mt"
      write (*,*) "= (4.35988/ 0.529178) x10^{-8} [J/m=Newton]"
      write (*,*) "= 8.2389668504738d-08 [N]"
      write (*,*)

      write (*,*) "area [Angst^2]= 1x10^{-20} mt^2"
      write (*,*)

      write (*,*) "pressure= radial force/area= 8.2389668504738d+12",
     &            " [N/m^2=Pa]"
      write (*,*)

      write (*,*) "facPres= 8.2389668504738d+12"
      write (*,*) "facPres * Pa= (Hartree/Bohr)/(Angst^2)"
      write (*,*)

      write (*,*) "Joule/Angst^3= Joule/[1x10^{-30} mt^3]"
      write (*,*) "             = 1x10^{30} Newton/mt^2"
      write (*,*) "             = 1x10^{30} Pascal"
      write (*,*)

      write (*,*) "facPres2= 1.x10^{30}"
      write (*,*) "facPres2 * Pa= Joule/Angst^3"
      write (*,*)

      write (*,*) "-------------"
      write (*,*) "units of kT"
      write (*,*) "-------------"
      write (*,*) "1 Joule= 2.29371248 x 10^{17}   Hartree"
      write (*,*) "kB= 1.3806504 x 10^{-23}       Joule/Kelvin"
      write (*,*) "  = 1.3806504x10^{-23} * 2.29371248x10^{17}",
     &            " Hartree/Kelvin"
      write (*,*) "  = 3.16681505 x 10^{-6}     Hartree/Kelvin"
      write (*,*)

      write (*,*) "[Temp]= Kelvin"
      write (*,*) "[kB*Temp]= Joule"
      write (*,*)

      write (*,*) "1 amu= 1.6605402 x10^{-27} Kg"
      write (*,*) "m/s= 1x10^{+10)/1x10^{+15} Angst/fs= 1x10^{-5}",
     &            " Angst/fs"
      write (*,*)

      write (*,*) "[kB*Temp/amu]= J/ (1.6605402x10^{-27} Kg)"
      write (*,*) "  = (m/s)^2/ (1.6605402x10^{-27})"
      write (*,*) "  = (Angst/fs)^2 * 1x10^{-10})/ (1.6605402x10^{-27})"
      write (*,*) "  = 1x10^{+17}/1.6605402 (Angst/fs)^2"
      write (*,*)

      write (*,*) "facKT= 1.0d17 /1.6605402"
      write (*,*) "facKT *(Angst/fs)^2= (Joule/Kelvin)*(Kelvin/amu)"
      write (*,*)

      write (*,*) "-------------------------"
      write (*,*) "units of spring constant K"
      write (*,*) "-------------------------"
      write (*,*) "[spring K]= Newton/cm"
      write (*,*) "1 amu= 1.6605402 x10^{-27} Kg"
      write (*,*)

      write (*,*) "Newton/cm= Kg*mt/(s^2*cm)"
      write (*,*) "         = Kg*100cm/ (1x10^{30} fs^2 *cm)"
      write (*,*) "         = Kg/(1x10^{28} fs^2)"
      write (*,*) "         = amu/(1.6605402 x10^{-27} *1x10^{28} fs^2)"
      write (*,*) "         = (amu/fs^2) *1.0/16.605402"
      write (*,*)

      write (*,*) "fspring= 1.0/ 16.605402"
      write (*,*) "Newton/cm = fspring *(amu/fs^2)"
      write (*,*)

      write (*,*) "-------------------------"
      write (*,*) "units of thermostat mass"
      write (*,*) "-------------------------"
      write (*,*) "kB= 1.3806504 x 10^{-23}       Joule/Kelvin"
      write (*,*) "1 Joule= 2.29371248 x 10^{17} Hartree"
      write (*,*) "1 mole= 6.02214179 x 10^{23}"
      write (*,*)

      write (*,*) "Mt= (KJoule/Mole) x ps^2"
      write (*,*) "  = 1000 x (2.29371248 x 10^{17} Hartree) x",
     &            " (10^{-12} sec)^2 /6.02214179 x 10^{23}"
      write (*,*) "  = (2.29371248/ 6.02214179) x 10^(3+17-24-23)",
     &            " Hartree * sec^2"
      write (*,*) "  = 3.808798530 x 10^{-27} Hartree * sec^2"
      write (*,*)

      write (*,*) "fmt= 3.808798530d-27"
      write (*,*) "(KJoule/Mole) x ps^2= fmt * (Hartree x sec^2)"
      write (*,*)

      write (*,*) "-------------------------"
      write (*,*) "units of frequency"
      write (*,*) "-------------------------"
      write (*,*) "(N/cm)/amu= (Kg mt/s^2)/ (cm *1.6605402x10^{-27} Kg)"
      write (*,*) "          = (100 cm/s^2)/(cm *1.6605402x10^{-27})"
      write (*,*) "          = (1/s^2)/1.6605402x10^{-29}"
      write (*,*)

      write (*,*) "[freq]= sqrt ([spring K/mass]"
      write (*,*) "      = sqrt [(N/cm)/amu]"
      write (*,*) "      = 1/sqrt(1.6605402x10^{-29}) s^{-1}"
      write (*,*)

      write (*,*) "facFreq= 1./sqrt(1.6605402d-29)"
      write (*,*) "sqrt[(N/cm)/amu]= facFreq * s^{-1}"
      write (*,*)

      write (*,*) "-------------------------------------------"
      write (*,*) "units for speed introducing the universal constant",
     &            " of gases"
      write (*,*) "-------------------------------------------"
      write (*,*) "universal constant of gases   R=8.3144621 J/",
     &            " (mol*Kelvin)"
      write (*,*) "Joule= newton*m= Kg*m/s^2 *m= Kg*m^2/s^2"
      write (*,*) "mole= 6.02214179 *10^{23}"
      write (*,*)

      write (*,*) "amu= 1.6605402 x10^{-27} Kg= 1.6605402*10^{-27} Kg",
     &            " *6.02214179*10^{23}/ mol"
      write (*,*) "     = 1.6605402 *6.02214179 *10^{-4} Kg/mol"
      write (*,*) "     = 1*10^{-3} Kg/mol"
      write (*,*)

      write (*,*) "[sqrt(RT/mass)]= sqrt [J/(mol*Kelvin) * Kelvin/",
     &            " (1*10^{-3}*Kg/mol)]"
      write (*,*) "               = sqrt [10^3 (m/s)^2]"
      write (*,*) "               = sqrt [10^3] m/s"
      write (*,*)

      write (*,*) "m/s= 1x10^{+10)/1x10^{+15) Angst/fs= 1x10^{-5}",
     &            " Angst/fs"
      write (*,*)
 
      write (*,*) "fac= 10^{-5}* sqrt(1000.)"
      write (*,*)

      write (*,*) "[sqrt(RT/m)]= sqrt [J/(mol*Kelvin) * Kelvin/amu]=",
     &            " fac *Angst/fs"
      write (*,*)

      write (*,*) "-----------------------------------"
      write (*,*) "convert energy to wavelength units"
      write (*,*) "-----------------------------------"
      write (*,*) "energy= hbar*freq"
      write (*,*) "lambda*freq= c"
      write (*,*)

      write (*,*) "energy= hbar*c/lambda"
      write (*,*) "lambda= hbar*c/energy"
      write (*,*)

      write (*,*) "hbar= 6.62606957 x10^{-34} (kg*m^2)/s"
      write (*,*) "c= 2.99792458 x10^8 m/s"
      write (*,*) "1 Hartree= 4.35988 x10^{-18} J=N*m=(kg*m/s^2)*m=",
     &            " kg*m^2/s^2"
      write (*,*) "alpha [hartree]= 4.35988 x10^{-18} *alpha [J]"
      write (*,*)

      write (*,*) "lambda= c*hbar/energy [(m/s)*((kg*m^2)/s) /J]= mt"
      write (*,*) "      = (2.99792458* 6.62606957/4.35988)",
     &            " x10^{8-34+18}/alpha  mt"
      write (*,*) "      = [(4.55619348)/alfa] x10^{-8} mt"
      write (*,*)

      write (*,*) "1 mt= 1x10^{+10} Angst"
      write (*,*) "lambda= [(455.619348)/alfa]  Angst"
      write (*,*)

      write (*,*) "calculator: frequency-wavelength-energy"
      write (*,*) "http://www2.chemistry.msu.edu/faculty/reusch",
     &            "/virttxtjml/cnvcalc.htm"
      write (*,*)

      write (*,*) "ionizing & non-ionizing radiation information"
      write (*,*) "http://www.epa.gov/radiation/understand",
     &            "/ionize_nonionize.html"
      write (*,*) "------------------"

      return
      end
* ****************
* ***********************************
* get the path of the working directory
* ***********************************
      subroutine workDir (name)
      implicit double precision (a-h,o-z)
      include 'parametros'

      character  qline*(maxchr), name*60
      dimension  itok(maxtok), jtok(maxtok)

      external kparse

* --------------------------
* parse the present working directory
* --------------------------

      open (unit=16, file="t1", status='unknown')
      open (unit=17, file="t2", status='unknown')

      read (16,12) qline

      ntok= kparse (qline,maxchr,itok,jtok,maxtok)
      write (17,12) char(34), qline (itok(1):jtok(1)), char(34)

      close (16)
      close (17)

* ---------------------
* read the directory name
* ---------------------

      open (unit=16, file="t2", status='unknown')
      read (16,12) name
      close (16)

   12 format (3a)

      return
      end
* *********************************
* ***********************************
* molec to be computed by terachem
* ***********************************
      subroutine writeMol (npar,atsymbol,x,y,z)
      implicit double precision (a-h,o-z)
      include 'parametros'

      character atsymbol(nmax)*5
      dimension  x(nmax), y(nmax), z(nmax)

      open (unit=14, file= 'mol.xyz', status= 'unknown')

      write(14,*) npar
      write(14,*) "molec"

      do m= 1, npar
      write(14,35) atsymbol(m), x(m),y(m),z(m)
c     write(*,35)  atsymbol(m), x(m),y(m),z(m)
      end do

   35 format (1x,a5, 3(2x,f10.3))
c  35 format (1x,a, 3(1x,f9.4))

      close (14)

      return
      end
* **********************************

* *************************
* file: randgen.f
* richard chandler
* (richard@stats.ucl.ac.uk)
* paul northrop
* (northrop@stats.ox.ac.uk)
* last modified: 26/8/03
* see file randgen.txt for details
* *************************
      block data zbqlbd01
 
* initializes seed array etc. for random number generator.
* the values below have themselves been generated using the
* nag generator.
 
      common /zbql0001/ zbqlix,b,c
      double precision zbqlix(43),b,c
      integer i

      data (zbqlix(i),i=1,43) /8.001441d7,5.5321801d8,
     +1.69570999d8,2.88589940d8,2.91581871d8,1.03842493d8,
     +7.9952507d7,3.81202335d8,3.11575334d8,4.02878631d8,
     +2.49757109d8,1.15192595d8,2.10629619d8,3.99952890d8,
     +4.12280521d8,1.33873288d8,7.1345525d7,2.23467704d8,
     +2.82934796d8,9.9756750d7,1.68564303d8,2.86817366d8,
     +1.14310713d8,3.47045253d8,9.3762426d7 ,1.09670477d8,
     +3.20029657d8,3.26369301d8,9.441177d6,3.53244738d8,
     +2.44771580d8,1.59804337d8,2.07319904d8,3.37342907d8,
     +3.75423178d8,7.0893571d7 ,4.26059785d8,3.95854390d8,
     +2.0081010d7,5.9250059d7,1.62176640d8,3.20429173d8,
     +2.63576576d8/

      data b / 4.294967291d9 /
      data c / 0.0d0 /

      end
* *************************
      subroutine zbqlini (seed)
* *************************
* to initialize the random number generator - either
* repeatably or nonrepeatably. need double precision
* variables because integer storage can't handle the
* numbers involved

* arguments
* =========
* seed     (integer, input). user-input number which generates
* elements of the array zbqlix, which is subsequently used
* in the random number generation algorithm. if seed=0,
* the array is seeded using the system clock if the
* fortran implementation allows it.

* parameters
* ==========
* lflno     (integer). number of lowest file handle to try when
* opening a temporary file to copy the system clock into.
* default is 80 to keep out of the way of any existing
* open files (although the program keeps searching till
* it finds an available handle). if this causes problems,
* (which will only happen if handles 80 through 99 are
* already in use), decrease the default value.

      integer lflno
      parameter (lflno=80)

* variables
* =========
* seed     see above
* zbqlix     seed array for the random number generator. defined
* in zbqlbd01
* b,c     used in congruential initialisation of zbqlix
* ss,mm,\}     system clock secs, mins, hours and days
* hh,dd \}
* filno     file handle used for temporary file
* init     indicates whether generator has already been initialised
     
      integer seed,ss,mm,hh,dd,filno,i
      integer init
      double precision zbqlix(43),b,c
      double precision tmpvar1,dss,dmm,dhh,ddd

      common /zbql0001/ zbqlix,b,c
      save init

     
* ensure we don't call this more than once in a program    
     
cc      if (init.ge.1) then
cc       if(init.eq.1) then
cc        write(*,1)
cc        init = 2
cc       endif
cc       return
cc      else
cc       init = 1
cc      endif
     
* if seed = 0, cat the contents of the clock into a file    
* and transform to obtain zqblix(1), then use a congr.    
* algorithm to set remaining elements. otherwise take    
* specified value of seed.    
     
*---------------------------------------------  
* nb for systems which do not support the             
* (non-standard) 'call system' command,               
* this will not work, and the first clause            
* of the following if block should be                 
* commented out.                             
*---------------------------------------------  
      if (seed.eq.0) then
*---------------------------------------------  
* comment out from here if you don't have             
* 'call system' capability ...                     
*---------------------------------------------  
       call system(' date +%s%m%h%j > zbql1234.tmp')
     
*       try all file numbers for lflno to 999     
     
       filno = lflno
 10    open(filno,file='zbql1234.tmp',err=11)
       goto 12
 11    filno = filno + 1
       if (filno.gt.999) then
        write(*,2)
        return
       endif
       goto 10
 12    read(filno,'(3(i2),i3)') ss,mm,hh,dd
       close(filno)
       call system('rm zbql1234.tmp')
       dss = dint((dble(ss)/6.0d1) * b)
       dmm = dint((dble(mm)/6.0d1) * b)
       dhh = dint((dble(hh)/2.4d1) * b)
       ddd = dint((dble(dd)/3.65d2) * b)
       tmpvar1 = dmod(dss+dmm+dhh+ddd,b)
*---------------------------------------------   
* ... to here (end of commenting out for                   
* tab users without 'call system' capability                  
*---------------------------------------------   
      else
       tmpvar1 = dmod(dble(seed),b)
      endif
      zbqlix(1) = tmpvar1
      do 100 i = 2,43
       tmpvar1 = zbqlix(i-1)*3.0269d4
       tmpvar1 = dmod(tmpvar1,b)      
       zbqlix(i) = tmpvar1
 100  continue

 1    format(//5x,'****warning**** you have called routine zbqlini ',
     +'more than',/5x,'once. i''m ignoring any subsequent calls.',//)
 2    format(//5x,'**** error **** in routine zbqlini, i couldn''t',
     +' find an',/5x,
     +'available file number. to rectify the problem, decrease the ',
     +'value of',/5x,
     +'the parameter lflno at the start of this routine (in file ',
     +'randgen.f)',/5x,
     +'and recompile. any number less than 100 should work.')
      end
*---------------------------------------------
      function zbqlu01()
 
*       returns a uniform random number between 0 & 1, using
*       a marsaglia-zaman type subtract-with-borrow generator.
*       uses double precision, rather than integer, arithmetic
*       throughout because mz's integer constants overflow
*       32-bit integer storage (which goes from -2^31 to 2^31).
*       ideally, we would explicitly truncate all integer
*       quantities at each stage to ensure that the double
*       precision representations do not accumulate approximation
*       error; however, on some machines the use of dnint to
*       accomplish this is *seriously* slow (run-time increased
*       by a factor of about 3). this double precision version
*       has been tested against an integer implementation that
*       uses long integers (non-standard and, again, slow) -
*       the output was identical up to the 16th decimal place
*       after 10^10 calls, so we're probably ok ...
 
      double precision zbqlu01,b,c,zbqlix(43),x,b2,binv
      integer curpos,id22,id43

      common /zbql0001/ zbqlix,b,c
      save /zbql0001/
      save curpos,id22,id43
      data curpos,id22,id43 /1,22,43/

      b2 = b
      binv = 1.0d0/b
 5    x = zbqlix(id22) - zbqlix(id43) - c
      if (x.lt.0.0d0) then
       x = x + b
       c = 1.0d0
      else
       c = 0.0d0
      endif
      zbqlix(id43) = x
 
* update array pointers. do explicit check for bounds of each to
* avoid expense of modular arithmetic. if one of them is 0 the others
* won't be
 
      curpos = curpos - 1
      id22 = id22 - 1
      id43 = id43 - 1
      if (curpos.eq.0) then
       curpos=43
      elseif (id22.eq.0) then
       id22 = 43
      elseif (id43.eq.0) then
       id43 = 43
      endif
     
* the integer arithmetic there can yield x=0, which can cause     
* problems in subsequent routines (e.g. zbqlexp). the problem    
* is simply that x is discrete whereas u is supposed to     
* be continuous - hence if x is 0, go back and generate another    
* x and return x/b^2 (etc.), which will be uniform on (0,1/b).     
     
      if (x.lt.binv) then
       b2 = b2*b
       goto 5
      endif

      zbqlu01 = x/b2

      end
* ***********************************

      DOUBLE PRECISION  FUNCTION V(X1, Y1, Z1,
     &  X2, Y2, Z2)


              DOUBLE PRECISION X1, X2, Y1, Y2, Z1, Z2
              DOUBLE PRECISION k, lo, l


              k = 20
c             lo = 1.5  

              lo = 0.9584

              l = sqrt((X2-X1)**2 + (Y2-Y1)**2 + 
     & (Z2-Z1)**2)


              V = (l-lo)**2

              RETURN

              END


       DOUBLE PRECISION  FUNCTION FRX(X1, Y1, Z1,
     & X2, Y2, Z2)


              DOUBLE PRECISION X1, X2, Y1, Y2, Z1, Z2
              DOUBLE PRECISION k, lo, l


c             lo = 1.5
              lo = 0.9584

              l = sqrt((X2-X1)**2 + (Y2-Y1)**2 + 
     & (Z2-Z1)**2)


              FRX = (l-lo)*(X2-X1)/l 

              RETURN

        END

 
 
       DOUBLE PRECISION  FUNCTION FRY(X1, Y1, Z1,
     & X2, Y2, Z2)


              DOUBLE PRECISION X1, X2, Y1, Y2, Z1, Z2
              DOUBLE PRECISION k, lo, l


              lo = 0.9584
c              lo = 1.5

              l = sqrt((X2-X1)**2 + (Y2-Y1)**2 + 
     & (Z2-Z1)**2)


              FRY = (l-lo)*(Y2-Y1)/l 

              RETURN

        END

              

       DOUBLE PRECISION  FUNCTION FRZ(X1, Y1, Z1,
     & X2, Y2, Z2)


              DOUBLE PRECISION X1, X2, Y1, Y2, Z1, Z2
              DOUBLE PRECISION k, lo, l


c             lo = 1.5
c             0.9584
              lo = 0.9584

              l = sqrt((X2-X1)**2 + (Y2-Y1)**2 + 
     & (Z2-Z1)**2)


              FRZ = (l-lo)*(Z2-Z1)/l 

              RETURN

        END










