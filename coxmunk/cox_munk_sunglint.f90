Module Cox_Munk


contains

subroutine sunglint( SZA,VZA,AZI,Tdown,Tup,m,wind_speed,AZI_wind,I_glint)
!Calcul la probailité de reflection speculaire pour un configuration
!d'observation
!Direct transmittances: Tdown,Tup

implicit none
real*8,intent(in) :: SZA,VZA,AZI,Tdown,Tup,m,wind_speed,AZI_wind
real*8,intent(out) :: I_glint
real*8 :: Q_glint,U_glint 
real*8 :: wind_speed_in,AZI_corr,sigma2_wind,thetaN,Pdist,theta_scat_in
real*8,parameter :: pi=3.141592653589793238462643d0,degrad=pi/180d0
real*8 :: omega,sigma1,sigma2,cos_sigma1,cos_sigma2
real*8,dimension(4,4) ::  Rf_pol,Tf_wa,L1,L2,Tup_mat
!real*8,parameter :: ind_refrac=1.334
real*8 :: ind_refrac
integer :: i_angle,i,j
real*8 :: xi,eta,sc2,su2,sc,su,c21,c03,c40,c22,c04
real*8 :: zx,zy,zx_prime,zy_prime
!---for SHADOW
logical :: SHADOW,breon
real*8 :: xnormg,snorm2,snorm,gamma0,gamma,cotan0,cotan,sdivc0,sdivc,arg0,arg,SH
real*8 :: delta
real*8 :: muv,mu0,muazi,sinv,sin0,sinazi

!!!!ATTENTION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! chgt aziangle pour convention 0 sol et sat opposés, 180 sol et sat meme cote
if(AZI>=0..and.AZI<180.)then
   AZI_corr=180-AZI
else
   AZI_corr=360+180-AZI
endif !(AZI)


muv=dcos(VZA*degrad)
sinv=dsin(VZA*degrad)
mu0=dcos(SZA*degrad)
sin0=dsin(SZA*degrad)
muazi=dcos(AZI*degrad)
sinazi=dsin(AZI*degrad)

ind_refrac=m

SHADOW=.true.
!SHADOW=.false.

breon=.true.
breon=.false.

wind_speed_in=wind_speed
!!_______Cox and Munk not defined for wind=0
if(wind_speed<0.1)then
  wind_speed_in=0.1
! I_glint=0d0
! Q_glint=0d0
! U_glint=0d0
! return
endif


!*---------------------------------------------------------------------*
!*           Compute scattered direction and surface normal
!*---------------------------------------------------------------------*

!Calcul de proba de reflex speculaire pour direction de visee thetaN
theta_scat_in=theta_scat_cm(SZA,VZA,AZI_corr)
omega=(180.-theta_scat_in)/2
!angle de reflexion speculaire = (scatt_angl-PI)/2
!print*,'angle reflex',((pi-theta_scat_in*degrad)/2.)/degrad
thetaN=dacos((mu0+muv)&
&           /(2.*dcos(omega*degrad)))


!*---------------------------------------------------------------------*
!*  Computation of the Fresnel matrix of reflection on the sea surface
!*---------------------------------------------------------------------*

call Fresnel(omega,ind_refrac,Rf_pol,Tf_wa)
!print*,'Rf_pol before rotation',Rf_pol(1,1)


!*---------------------------------------------------------------------*
!* Compute ocean surface slope distribution 'Pdist'                    *
!*---------------------------------------------------------------------*

!------------
select case (1)
case(1)
 ! WU formula [Wu, 1990]
  if(wind_speed_in<7)then
   sigma2_wind=0.0276*dlog10(wind_speed_in)+0.009
  else
   sigma2_wind=0.138*dlog10(wind_speed_in)-0.084
  endif

case(2)
 ! Original Cox Munk
  sigma2_wind=3d-3+512d-5*wind_speed_in

case(3)
 ! Hu et al 2008 based on CALIPSO / AMSR-E comparison
  if(wind_speed_in<7)then
   sigma2_wind=1.46d-2*dsqrt(wind_speed_in)
  elseif(wind_speed_in>=7 .and. wind_speed_in<13.3)then
   sigma2_wind=3d-3+512d-5*wind_speed_in
  else
   sigma2_wind=0.138*dlog10(wind_speed_in)-0.084
  endif
end select



if(.false.)then
!--- sigma2 cox Munk f(wind_speed_in)
 Pdist=1./(pi*sigma2_wind)*&
&                 exp(-1.*(dtan(thetaN)**2/sigma2_wind))

else
if(breon)then
! from Breon Henriot 2006 JGR
 sc2    = 3d-3 + 1.85d-3*wind_speed_in
 su2    = 1d-3 + 3.16d-3*wind_speed_in
 sc     = dsqrt(sc2)
 su     = dsqrt(su2)
 c21    = -9d-4*wind_speed_in**2d0
 c03    = -0.45d0/(1d0+dexp(7d0-wind_speed_in))
 c40    = 0.30d0
 c22    = 0.12d0
 c04    = 0.4d0
else
!historical values from COX MUNK
 sc2    = 0.003d0 + 1.92d-3*wind_speed_in
 su2    =           3.16d-3*wind_speed_in
 sc     = dsqrt(sc2)
 su     = dsqrt(su2)
 c21    = 0.01d0 - 8.6d-3*wind_speed_in
 c03    = 0.04d0 - 33.d-3*wind_speed_in
 c40    = 0.40d0
 c22    = 0.12d0
 c04    = 0.23d0
endif
 sigma2_wind=sc2+su2

 zx=(-sinv*sinazi)/(mu0+muv)
 zy=(sinv*muazi+sin0)/(mu0+muv)

 zx_prime=dcos(AZI_wind*degrad)*zx+dsin(AZI_wind*degrad)*zy
 zy_prime=-dsin(AZI_wind*degrad)*zx+dcos(AZI_wind*degrad)*zy

 xi = zx_prime/sc
 eta= zy_prime/su
 
 Pdist = dexp(-5d-1*(xi**2+eta**2))/(2.d0*pi*sc*su) *    &
&                ( 1.d0     -                            &
&                  c21*(xi**2-1.d0)*eta/2.d0 -           &
&                  c03*(eta**3-3.d0*eta)/6.d0 +          &
&                  c40*(xi**4-6.d0*eta**2+3.d0)/24.d0 +  &
&                  c04*(eta**4-6.d0*eta**2+3.d0)/24.d0 + &
&                  c22*(xi**2-1.d0)*(eta**2-1.d0)/4.d0 )
! print*,'xi,eta',xi,eta,sc*su,Pdist

 if(Pdist<=epsilon(Pdist))Pdist=0d0
!warning very crude regularization
! if(Pdist>1d0)Pdist=1d0 
endif
!print*,'P1',Pdist,dexp(-5d-1*(xi**2+eta**2)),(2.d0*pi*sc*su),( 1.d0     -                            &
!&                  c21*(xi**2-1.d0)*eta/2.d0 -           &
!&                  c03*(eta**3-3.d0*eta)/6.d0 +          &
!&                  c40*(xi**4-6.d0*eta**2+3.d0)/24.d0 +  &
!&                  c04*(eta**4-6.d0*eta**2+3.d0)/24.d0 + &
!&                  c22*(xi**2-1.d0)*(eta**2-1.d0)/4.d0 )

!*----------------------------------------------------------------------*
!*                           Direct Transmittance
!*---------------------------------------------------------------------*
!Tdown=exp(-1*(OptiThick)/dcos(SZA*degrad))
!Tup=exp(-1*(OptiThick)/dcos(VZA*degrad))

!*---------------------------------------------------------------------*
!*                    Rotation in the reference plane
!*---------------------------------------------------------------------*

cos_sigma1=(dcos(pi-SZA*degrad)-muv*dcos(theta_scat_in*degrad))/(sinv*dsin(theta_scat_in*degrad))
cos_sigma2=(muv-dcos(pi-SZA*degrad)*dcos(theta_scat_in*degrad))/(dsin(pi-SZA*degrad)*dsin(theta_scat_in*degrad))

sigma1=1d0*dacos(cos_sigma1)
if(AZI>=0..and.AZI<180.)then
   sigma2=1d0*dacos(cos_sigma2)
!   sigma1=-1d0*dacos(cos_sigma1)
else
   sigma2=-1d0*dacos(cos_sigma2)
   sigma1=-1d0*dacos(cos_sigma1)
endif !(AZI)


L1=0.
L2=0.
L1(1,1)=1
L1(4,4)=1
L2(1,1)=1
L2(4,4)=1

L1(2,2)=dcos(-2d0*sigma1)
L1(3,3)=L1(2,2)
L2(2,2)=dcos(2d0*(pi-sigma2))
L2(3,3)=L2(2,2)

L1(2,3)=dsin(-2d0*sigma1)
L1(3,2)=-1d0*L1(2,3)
L2(2,3)=dsin(2d0*(pi-sigma2))
L2(3,2)=-1d0*L2(2,3)

Tup_mat=0
!print*,'Tup before rot ',Tup
!Tup_mat(1,1)=1d0-Tup
!Tup_mat(2,2)=1d0-Tup
!Tup_mat(3,3)=1d0-Tup
!Tup_mat(4,4)=1d0-Tup
Tup_mat(1,1)=Tup
Tup_mat(2,2)=Tup
Tup_mat(3,3)=Tup
Tup_mat(4,4)=Tup


!print*,'L1 ',L1
!print*,'L2 ',L2
if(AZI/=0..and.AZI/=180..and.AZI/=360..and.VZA/=0.)then
 Rf_pol=matmul(Rf_pol,L1)
 Rf_pol=matmul(L2,Rf_pol)
! Tup_mat=matmul(Tup_mat,L1)
! Tup_mat=matmul(L2,Tup_mat)
endif
! Conversion of extinction matrix into transmission matrix
!Tup_mat=Tup_mat
!print*,'Tup after rot ',Tup_mat(1,1),Tup_mat(1,2),Tup_mat(1,3)

!print*,'rot Tup12=',Tup_mat(3,:)

!----------------------------------------------------------------------
!               Compute the shadowing factor 'SH'                      
!----------------------------------------------------------------------
      if (SHADOW) then
         xnormg = dsqrt(2.D0/pi)
         snorm2 = sigma2_wind/2.D0
         snorm  = dsqrt(snorm2)
         if (SZA.eq.0.D0) then
            gamma0 = 0.D0
         else
            cotan0 = mu0/sin0
            sdivc0 = snorm/cotan0
            arg0   = cotan0/(dsqrt(2.D0)*snorm)
            gamma0 =(xnormg*sdivc0*dexp(-arg0*arg0)-derfc(arg0))/2.D0
         endif
         if (VZA.eq.0.D0) then
            gamma  = 0.D0
         else
            cotan  = muv/sinv
            sdivc  = snorm/cotan
            arg    = cotan/(dsqrt(2.D0)*snorm)
            gamma  =(xnormg*sdivc*dexp(-arg*arg) - derfc(arg))/2.D0
         endif
         SH       = 1.D0/( gamma + gamma0 + 1.D0)
      else
         SH       = 1.D0
      endif
 
!*---------------------------------------------------------------------*
!*        Forward Scattering Component for Sun glint Stokes Single scatt 
!*---------------------------------------------------------------------*
!Rayleigh approximation
!Pforward=0d0
!Pforward(1,1)=1.5d0
!Pforward(2,2)=1.5d0
!Pforward(3,3)=1.5d0
!Pforward(4,4)=1.5d0


!*---------------------------------------------------------------------*
!*        Sun glint Stokes component (in reflectance unit) at TOA
!*---------------------------------------------------------------------*
Pdist=SH*(pi*Pdist)/(4.*mu0*muv*dcos(thetaN)**4)
!print*,'P,R',Pdist,Rf_pol(1,1),(4.*mu0*muv*dcos(thetaN)**4)
! I_glint=Tdown*(Tup+0.015*(1+delta/2))*Rf_pol(1,1)&
 I_glint=Tdown*Tup*Rf_pol(1,1)&
! I_glint=Tdown*(Tup_mat(1,1)*Rf_pol(1,1)+Tup_mat(2,1)*Rf_pol(2,1)+Tup_mat(3,1)*Rf_pol(3,1))&
! I_glint=Tdown*(Tup_mat(1,1)*Rf_pol(1,1)+Tup_mat(1,2)*Rf_pol(2,1)+Tup_mat(1,3)*Rf_pol(3,1))&
&        *Pdist

 if(I_glint<0)then
  print*,'**** **** **** ****'
  print*,'I_glint = ',I_glint
 endif

! Q_glint=Tdown*(Tup+0.015*3/2*delta)*Rf_pol(1,2)&
 Q_glint=Tdown*Tup*Rf_pol(2,1)&
! Q_glint=Tdown*(Tup_mat(1,2)*Rf_pol(1,1)+Tup_mat(2,2)*Rf_pol(2,1)+Tup_mat(3,2)*Rf_pol(3,1))&
! Q_glint=Tdown*(Tup_mat(2,1)*Rf_pol(1,1)+Tup_mat(2,2)*Rf_pol(2,1)+Tup_mat(2,3)*Rf_pol(3,1))&
&        *Pdist

! U_glint=Tdown*(Tup+0.015*3/2*delta)*Rf_pol(1,3)&
 U_glint=Tdown*Tup*Rf_pol(3,1)&
! U_glint=Tdown*(Tup_mat(1,3)*Rf_pol(1,1)+Tup_mat(2,3)*Rf_pol(2,1)+Tup_mat(3,3)*Rf_pol(3,1))&
! U_glint=Tdown*(Tup_mat(3,1)*Rf_pol(1,1)+Tup_mat(3,2)*Rf_pol(2,1)+Tup_mat(3,3)*Rf_pol(3,1))&
&        *Pdist


!-- to avoid artifact in computation of derivatives (e.g., Jacobian)
!if(I_glint < 1d-8)I_glint=0d0
!if(Q_glint < 1d-8)Q_glint=0d0
!if(U_glint < 1d-8)U_glint=0d0

!print*,'Rf_pol',Rf_pol

return
end subroutine


subroutine Fresnel(angle_deg,m,Rf_pol,Tf)
!Computation of elements of the Relection Fresnel Matrix in Stokes
! formalism (V component is assumed negligible)

! input: angle: refraction angle in degrees; m : refractive index
! Output:
!        Rf_pol = reflexion matrix for illumination from above
!        Tf = Transmission matrix for illumination from below

implicit none
!!!!!!!!angle en degres
real*8,intent(in)  :: angle_deg,m
real*8,dimension(4,4),intent(out)  :: Rf_pol,Tf
real*8             :: mu,racine,rl,rr,angle_below,cos_above,tl,tr
real*8             :: degrad,angle
degrad=datan(1d0)/45.


angle=angle_deg*degrad
racine=sqrt(m**2-(dsin(angle))**2);
rl=(racine-m**2*dcos(angle))/(racine+m**2*dcos(angle));
rr=(dcos(angle)-racine)/(dcos(angle)+racine);

!        Rf_pol = reflexion matrix for illumination from above
Rf_pol=0.;
Rf_pol(1,1)=rl**2+rr**2;
Rf_pol(2,2)=Rf_pol(1,1);
Rf_pol(1,2)=rl**2-rr**2;
Rf_pol(2,1)=Rf_pol(1,2);
Rf_pol(3,3)=2*rl*rr;
Rf_pol(4,4)=2*rl*rr;

Rf_pol=0.5*Rf_pol;


!        Tf_wa = Transmission matrix for illumination from below
angle_below=dasin(dsin(angle)/m);
cos_above=sqrt(1-(dsin(angle_below))**2*m**2);
tl=(2*m*dcos(angle_below))/(dcos(angle_below)+m*cos_above);
tr=(2*m*dcos(angle_below))/(m*dcos(angle_below)+cos_above);

!Tf(1:4,1:4)=0.;
Tf(1,1)=tl**2+tr**2;
Tf(2,2)=tl**2+tr**2;
Tf(1,2)=tl**2-tr**2;
Tf(2,1)=tl**2-tr**2;
Tf(3,3)=2*tl*tr;
Tf(4,4)=2*tl*tr;
Tf=1/(4*m*dcos(angle_below))*Tf;

Tf=2*dcos(angle)/m**2*Tf;
return
end subroutine



function theta_scat_cm(SZA_in,VZA_in,AZI_in)
implicit none


real*8,intent(in)  :: SZA_in,VZA_in,AZI_in
real*8             :: degrad,theta_scat_cm
degrad=datan(1d0)/45.

         theta_scat_cm=dcos((180d0-SZA_in)*degrad)*dcos(VZA_in*degrad)+dsin((180d0-SZA_in)*degrad)*sin(VZA_in*degrad)&
&                  *dcos(AZI_in*degrad)
         theta_scat_cm=dacos(1d0*theta_scat_cm)/degrad
return
end function theta_scat_cm

end Module


