program twostars2_tst
!----------------------------------------------------------------------------------------
! sample program showing the use of the "twostars_m" module to compute the top of the
! atmosphere insolation on a circumbinary planet 
!
! same as twostars1, but with direct input of effective temperatures and stellar radii
!
! written by Siegfried Eggl 20140610
!
! last modified by Siegfried Eggl 20180409 
!   > modified output to csv
!
!
! dependencies: module twostars_m
!
! Input: file "twostars2.inn"
!
! Output: file "twostars2_out.csv"
!
!----------------------------------------------------------------------------------------
use twostars_m

implicit none
!define
integer, parameter :: dp = selected_real_kind(15, 307) 	  !precision parameter
integer::nbands						  !number of spectral bands
real(kind=dp)::e1						!eccentricity of the stellar binary
real(kind=dp),dimension(1:2)::teff,rs,semia,man		!stellar effective temperatures, stellar radii, semimajor axes of binary and planetary orbits, initial mean longitudes 
real(kind=dp),dimension(1:3)::mass				!system masses: m(1:2) stars, m(3) planet
real(kind=dp),dimension(1:4)::spin				!spin quantities: (1) obliquity [deg] (2) angle of precession [deg] (3) rotation period [D] (4) initial hour angle [deg] (=initial position of continents wrt x-axis) 
real(kind=dp),dimension(:,:,:),allocatable::bands		!spectral band information (upper and lower limits per band per star for each band)
real(kind=dp)::t,tend,dt					!time, end time, output timestep
 
!results
integer::ok							!configuration check parameter: 0-> everything is ok, else...
real(kind=dp),dimension(1:2)::decl,phi				!declination and hour angle wrt each star
real(kind=dp),dimension(:,:),allocatable::bndflux		!band fluxes for each star and each band

character(len=200)::fmt_out !changed 20180409 

 !-----------------------------------------------------------------
 !read input file
 !----------------------------------------------------------------
 open(21,file="twostars2.inn",status='old')
 read(21,*)tend,dt	  !end time, output timestep
 read(21,*)mass           !masses of star1,star2,planet [Msun]
 read(21,*)teff           !THIS IS DIFFERENT: effective temperatures of both stars [K]
 read(21,*)rs		  !THIS IS DIFFERENT: stellar radii [Rsun]
 read(21,*)semia 	  !binary and planetary semimajor axes [au]
 read(21,*)e1		  !binary eccentricity
 read(21,*)man		  !initial mean longitudes [rad] (starting positions of the planet and the stars on their orbits)
 read(21,*)spin           !planetary spin quantities: (1) obliquity [deg] (2) angle of precession [deg] (3) rotation period [D] (4) initial hour angle [deg] (=initial position of continents wrt x-axis) 
 read(21,*)nbands	  !number of spectral bands
 
 allocate(bands(2,nbands,2),bndflux(2,nbands))
 read(21,*)bands(1,:,1)	  	!... [1/cm] (1,1:nbands,1) lower spectral band limits for star 1
 read(21,*)bands(1,:,2)    	!... [1/cm] (1,1:nbands,2) upper spectral band limits for star 1
 read(21,*)bands(2,:,1)		!... [1/cm] (2,1:nbands,1) lower spectral band limits for star 2
 read(21,*)bands(2,:,2)		!... [1/cm] (2,1:nbands,1) upper spectral band limits for star 2
 close(21)
 
 !-----------------------------------------------------------------
 !check whether the dynamical configuration is possible (stellar radii and effective temperatures are not checked)
 !----------------------------------------------------------------
 call checkconf(mass,semia,e1,ok)
 write(*,*)'configuration ok?',ok
 if(ok.gt.0.and.ok.le.4) STOP


 !produce the output file with the insolation2 subroutine
 open(22,file="twostars2_out.csv",status='replace')                 !changed 20180409 
 !write header into ouput file.
 write(22,*)'1: time [D], 2: declination wrt star1 [rad], 3: declination wrt star2 [rad],',&
 '4: hour angle wrt star 1 [rad], 5: hour angle wrt star2 [rad],', &
 '6: sum flux star1 [W/m^2], 7: sum flux star2 [W/m^2]'

fmt_out="(6 (F21.15, ', '), F21.15)"                                        !changed 20180409 

 !do time evolution to see the change in insolation
 t=0.d0
 do while(t.lt.tend)
  !call module function
  call insolation2(mass,teff,rs,semia,e1,man,spin,t,nbands,bands,bndflux,decl,phi) 
  !write output at each timestep
  write(22,fmt=fmt_out)t,decl(:),phi(:),sum(bndflux(1,:)),sum(bndflux(2,:))     !changed 20180409 
  t=t+dt 
 end do
 
close(22)
end program 
