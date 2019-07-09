

      subroutine lsdsr(rs,z,mu,excsr,vxcsrup,vxcsrdown)
ccc Hartree atomic units used
ccc for given density parameter 'rs', spin polarization 'z'
ccc and cutoff parameter 'mu' 
ccc gives the complementary  short-range exchange-correlation
ccc energy  (i.e., xc energy of jellium minus xc energy of long-range
ccc interacting electron gas) => 'excsr'
ccc and the corresponding exchange-correlation potentials for
ccc spin-up and spin-down electrons => 'vxcsrup','vxcsrdown'
ccc from Paziani, Moroni, Gori-Giorgi, and Bachelet, cond-mat/0601353
      implicit none
      double precision rs,z,mu,excsr,vxcsrup,vxcsrdown
      double precision eclr,exlr,ec,ecd,ecz,ex
      double precision vclrup,vclrdown,vxlrup,vxlrdown
      double precision vxup,vxdown,vcup,vcdown
      double precision pi,alpha,cf
      pi=dacos(-1.d0)
      alpha=(4.d0/9.d0/pi)**(1.d0/3.d0)
      cf=1.d0/alpha

      ex=-3.d0*cf/rs/8.d0/pi*((1.d0+z)**(4.d0/3.d0)+
     $     (1.d0-z)**(4.d0/3.d0))
      ex=0d0

      vxup=-(1.d0+z)**(1.d0/3.d0)*(3.d0/2.d0/pi)**(2.d0/3.d0)/rs
      vxdown=-(1.d0-z)**(1.d0/3.d0)*(3.d0/2.d0/pi)**(2.d0/3.d0)/rs
      vxup = 0d0
      vxdown = 0d0

      call ecPW(rs,z,ec,ecd,ecz)
      vcup=ec-rs/3.d0*ecd-(z-1.d0)*ecz
      vcdown=ec-rs/3.d0*ecd-(z+1.d0)*ecz

      call exchangelr(rs,z,mu,exlr)
      exlr = 0d0
      call vexchangelr(rs,z,mu,vxlrup,vxlrdown)
      vxlrup = 0d0
      vxlrdown = 0d0
      call ecorrlr(rs,z,mu,eclr)
      call vcorrlr(rs,z,mu,vclrup,vclrdown)
      excsr=ex+ec-(exlr+eclr)
      vxcsrup=vxup+vcup-(vxlrup+vclrup)
      vxcsrdown=vxdown+vcdown-(vxlrdown+vclrdown)

      return
      end



      subroutine ecorrlr(rs,z,mu,eclr)
ccc Hartree atomic units used
ccc for given density parameter rs, spin polarization z
ccc and cutoff parameter mu 
ccc gives the correlation energy of the LR gas
ccc  => eclr
      implicit none
      double precision rs,z,mu,eclr,ec,ecd,ecz
      double precision pi,alpha,cf,phi
      double precision g0,dpol,d2anti,d3anti,Qrpa
      double precision coe2,coe3,coe4,coe5
      double precision a1,a2,a3,a4,b0
      double precision q1a,q2a,q3a,t1a,t2a,t3a,adib
      pi=dacos(-1.d0)
      alpha=(4.d0/9.d0/pi)**(1.d0/3.d0)
      cf=1.d0/alpha

      phi=((1.d0+z)**(2.d0/3.d0)+(1.d0-z)**(2.d0/3.d0))/2.d0
cc parameters from the fit
      adib   = 0.784949d0   
      q1a    = -0.388d0   
      q2a    = 0.676d0   
      q3a    = 0.547d0   
      t1a    = -4.95d0   
      t2a    = 1.d0    
      t3a    = 0.31d0   

      b0=adib*rs

      d2anti=(q1a*rs+q2a*rs**2)*exp(-abs(q3a)*rs)/rs**2
      d3anti=(t1a*rs+t2a*rs**2)*exp(-abs(t3a)*rs)/rs**3

      coe2=-3.d0/8.d0/rs**3*(1.d0-z**2)*(g0(rs)-0.5d0)

      coe3=-(1.d0-z**2)*g0(rs)/(sqrt(2.d0*pi)*rs**3)

      if(abs(z).eq.1.d0) then

        coe4=-9.d0/64.d0/rs**3*(dpol(rs)
     $        -cf**2*2**(5.d0/3.d0)/5.d0/rs**2) 
        coe5=-9.d0/40.d0/(sqrt(2.d0*pi)*rs**3)*dpol(rs)

      else

         coe4=-9.d0/64.d0/rs**3*(((1.d0+z)/2.d0)**2*
     $        dpol(rs*(2/(1.d0+z))**(1.d0/3.d0))+((1.d0-z)/2.d0)**2
     $        *dpol(rs*(2.d0/(1.d0-z))**(1.d0/3.d0))+
     $        (1.-z**2)*d2anti-cf**2/10.d0*((1.d0+z)**(8.d0/3.d0)
     $        +(1.-z)**(8.d0/3.d0))/rs**2)

         coe5=-9.d0/40.d0/(sqrt(2.d0*pi)*rs**3)*(((1.d0+z)/2.d0)**2
     $        *dpol(rs*(2.d0/(1.d0+z))**(1.d0/3.d0))+((1.d0-z)/2.d0)**2
     $        *dpol(rs*(2.d0/(1.d0-z))**(1.d0/3.d0))+(1.d0-z**2)*
     $        d3anti)
      endif

      call ecPW(rs,z,ec,ecd,ecz)

      a1=4.d0*b0**6*coe3+b0**8*coe5
      a2=4.d0*b0**6*coe2+b0**8*coe4+6.d0*b0**4*ec
      a3=b0**8*coe3
      a4=b0**6*(b0**2*coe2+4.d0*ec)
      
      eclr=(phi**3*Qrpa(mu*sqrt(rs)/phi)+a1*mu**3+a2*mu**4+a3*mu**5+
     $     a4*mu**6+b0**8*mu**8*ec)/((1.d0+b0**2*mu**2)**4)

      return
      end

      subroutine vcorrlr(rs,z,mu,vclrup,vclrdown)
ccc Hartree atomic units used
ccc for given density parameter rs, spin polarization z
ccc and cutoff mu it gives the correlation LSD potential for LR interaction
ccc  => vclrup (spin-up electrons), vclrdown (spin-down electrons)
      implicit none
      double precision rs,z,mu,eclr,eclrrs,eclrz,vclrup,vclrdown
      double precision ec,ecd,ecz
      double precision pi,alpha,cf,phi
      double precision g0,dpol,d2anti,d3anti,Qrpa
      double precision g0d,dpold,d2antid,d3antid,Qrpad,x
      double precision coe2,coe3,coe4,coe5
      double precision coe2rs,coe3rs,coe4rs,coe5rs
      double precision coe2z,coe3z,coe4z,coe5z
      double precision a1,a2,a3,a4,a5,b0,a1rs,a2rs,a3rs,a4rs,a5rs,
     $     b0rs,a1z,a2z,a3z,a4z,a5z,b0z
      double precision q1a,q2a,q3a,t1a,t2a,t3a,adib
      
      pi=dacos(-1.d0)
      alpha=(4.d0/9.d0/pi)**(1.d0/3.d0)
      cf=1.d0/alpha

      phi=((1.d0+z)**(2.d0/3.d0)+(1.d0-z)**(2.d0/3.d0))/2.d0
cc parameters from the fit
      adib   = 0.784949d0   
      q1a    = -0.388d0   
      q2a    = 0.676d0   
      q3a    = 0.547d0   
      t1a    = -4.95d0   
      t2a    = 1.d0    
      t3a    = 0.31d0   

      b0=adib*rs

      d2anti=(q1a+q2a*rs)*exp(-q3a*rs)/rs
      d3anti=(t1a+t2a*rs)*exp(-t3a*rs)/rs**2

      d2antid=-((q1a + q1a*q3a*rs + q2a*q3a*rs**2)/
     -    rs**2)*exp(-q3a*rs)
      d3antid=-((rs*t2a*(1 + rs*t3a) + t1a*(2 + rs*t3a))/
     -    rs**3)*exp(-rs*t3a)

      coe2=-3.d0/8.d0/rs**3*(1.d0-z**2)*(g0(rs)-0.5d0)
      coe2rs=-3.d0/8.d0/rs**3*(1.d0-z**2)*g0d(rs)+
     $     9.d0/8.d0/rs**4*(1.d0-z**2)*(g0(rs)-0.5d0)
      coe2z=-3.d0/8.d0/rs**3*(-2.d0*z)*(g0(rs)-0.5d0)

      coe3=-(1.d0-z**2)*g0(rs)/(sqrt(2.d0*pi)*rs**3)
      coe3rs=-(1.d0-z**2)*g0d(rs)/(sqrt(2.d0*pi)*rs**3)+
     $    3.d0*(1.d0-z**2)*g0(rs)/(sqrt(2.d0*pi)*rs**4) 
      coe3z=2.d0*z*g0(rs)/(sqrt(2.d0*pi)*rs**3)

      if(abs(z).eq.1.d0) then

        coe4=-9.d0/64.d0/rs**3*(dpol(rs)
     $        -cf**2*2**(5.d0/3.d0)/5.d0/rs**2)
        coe4rs=-3.d0/rs*coe4-9.d0/64.d0/rs**3*(dpold(rs)
     $        +2.d0*cf**2*2**(5.d0/3.d0)/5.d0/rs**3)
        coe4z=-9.d0/64.d0/rs**3*(dpol(rs)-rs/6.d0*dpold(rs)-2.d0*d2anti
     $       -4.d0/15.d0/rs**2*cf**2*2.d0**(5.d0/3.d0))*z
        coe5=-9.d0/40.d0/(sqrt(2.d0*pi)*rs**3)*dpol(rs)
        coe5rs=-3.d0/rs*coe5-9.d0/40.d0/(sqrt(2.d0*pi)*rs**3)*dpold(rs)
        coe5z=-9.d0/40.d0/(sqrt(2.d0*pi)*rs**3)*(dpol(rs)-rs/6.d0*
     $       dpold(rs)-2.d0*d3anti)*z

      else

         coe4=-9.d0/64.d0/rs**3*(((1.d0+z)/2.d0)**2*
     $        dpol(rs*(2/(1.d0+z))**(1.d0/3.d0))+((1.d0-z)/2.d0)**2
     $        *dpol(rs*(2.d0/(1.d0-z))**(1.d0/3.d0))+
     $        (1.-z**2)*d2anti-cf**2/10.d0*((1.d0+z)**(8.d0/3.d0)
     $        +(1.-z)**(8.d0/3.d0))/rs**2)
         coe4rs=-3.d0/rs*coe4-9.d0/64.d0/rs**3*(
     $        ((1.d0+z)/2.d0)**(5.d0/3.d0)*dpold(rs*(2/(1.d0+z))**
     $        (1.d0/3.d0))+((1.d0-z)/2.d0)**(5.d0/3.d0)*
     $        dpold(rs*(2/(1.d0-z))**(1.d0/3.d0))+(1.d0-z**2)*
     $        d2antid+cf**2/5.d0*((1.d0+z)**(8.d0/3.d0)
     $        +(1.d0-z)**(8.d0/3.d0))/rs**3)
         coe4z=-9.d0/64.d0/rs**3*(1.d0/2.d0*(1.d0+z)*
     $        dpol(rs*(2/(1.d0+z))**(1.d0/3.d0))-1.d0/2.d0*(1.d0-z)*
     $        dpol(rs*(2/(1.d0-z))**(1.d0/3.d0))-rs/6.d0*
     $        ((1.d0+z)/2.d0)**(2.d0/3.d0)*dpold(rs*(2/(1.d0+z))
     $        **(1.d0/3.d0))+rs/6.d0*((1.d0-z)/2.d0)**(2.d0/3.d0)
     $        *dpold(rs*(2/(1.d0-z))**(1.d0/3.d0))-2.d0*z*d2anti-
     $        4.d0/15.d0/rs**2*cf**2*((1.d0+z)**(5.d0/3.d0)-
     $        (1.d0-z)**(5.d0/3.d0)))

         coe5=-9.d0/40.d0/(sqrt(2.d0*pi)*rs**3)*(((1.d0+z)/2.d0)**2
     $        *dpol(rs*(2.d0/(1.d0+z))**(1.d0/3.d0))+((1.d0-z)/2.d0)**2
     $        *dpol(rs*(2.d0/(1.d0-z))**(1.d0/3.d0))+(1.d0-z**2)*
     $        d3anti)
         coe5rs=-3.d0/rs*coe5-9.d0/(40.d0*sqrt(2.d0*pi)*rs**3)*(
     $        ((1.d0+z)/2.d0)**(5.d0/3.d0)*dpold(rs*(2/(1.d0+z))**
     $        (1.d0/3.d0))+((1.d0-z)/2.d0)**(5.d0/3.d0)*
     $        dpold(rs*(2/(1.d0-z))**(1.d0/3.d0))+(1.d0-z**2)*
     $        d3antid)
         coe5z=-9.d0/40.d0/(sqrt(2.d0*pi)*rs**3)*(1.d0/2.d0*(1.d0+z)*
     $        dpol(rs*(2/(1.d0+z))**(1.d0/3.d0))-1.d0/2.d0*(1.d0-z)*
     $        dpol(rs*(2/(1.d0-z))**(1.d0/3.d0))-rs/6.d0*
     $        ((1.d0+z)/2.d0)**(2.d0/3.d0)*dpold(rs*(2/(1.d0+z))
     $        **(1.d0/3.d0))+rs/6.d0*((1.d0-z)/2.d0)**(2.d0/3.d0)
     $        *dpold(rs*(2/(1.d0-z))**(1.d0/3.d0))-2.d0*z*d3anti)

      endif

      call ecPW(rs,z,ec,ecd,ecz)

      a1=4.d0*b0**6*coe3+b0**8*coe5
      a1rs=24.d0*adib*b0**5*coe3+4.d0*b0**6*coe3rs+8.d0*adib*b0**7*
     $     coe5+b0**8*coe5rs
      a1z=4.d0*b0**6*coe3z+b0**8*coe5z

      a2=4.d0*b0**6*coe2+b0**8*coe4+6.d0*b0**4*ec
      a2rs=24.d0*adib*b0**5*coe2+4.d0*b0**6*coe2rs+8.d0*adib*b0**7*
     $     coe4+b0**8*coe4rs+24.d0*adib*b0**3*ec+6.d0*b0**4*ecd
      a2z=4.d0*b0**6*coe2z+b0**8*coe4z+6.d0*b0**4*ecz

      a3=b0**8*coe3
      a3rs=8.d0*adib*b0**7*coe3+b0**8*coe3rs
      a3z=b0**8*coe3z

      a4=b0**6*(b0**2*coe2+4.d0*ec)
      a4rs=8.d0*adib*b0**7*coe2+b0**8*coe2rs+24.d0*adib*b0**5*ec+
     $     4.d0*b0**6*ecd
      a4z=b0**6*(b0**2*coe2z+4.d0*ecz)

      a5=b0**8*ec
      a5rs=8.d0*adib*b0**7*ec+b0**8*ecd
      a5z=b0**8*ecz

      x=mu*sqrt(rs)/phi

      eclr=(phi**3*Qrpa(x)+a1*mu**3+a2*mu**4+a3*mu**5+
     $     a4*mu**6+a5*mu**8)/((1.d0+b0**2*mu**2)**4)
      
      eclrrs=-4.d0/(1.d0+b0**2*mu**2)*2.d0*adib*b0*mu**2*eclr+
     $     1.d0/((1.d0+b0**2*mu**2)**4)*(phi**2*mu/(2.d0*sqrt(rs))
     $     *Qrpad(x)+
     $     a1rs*mu**3+a2rs*mu**4+a3rs*mu**5+a4rs*mu**6+a5rs*mu**8)


      if(z.eq.1.d0) then
         vclrup=eclr-rs/3.d0*eclrrs
         vclrdown=0.d0
      elseif(z.eq.-1.d0) then
         vclrup=0.d0
         vclrdown=eclr-rs/3.d0*eclrrs
      else

         eclrz=(phi**2*((1.d0+z)**(-1.d0/3.d0)-(1.d0-z)**(-1.d0/3.d0))
     $        *Qrpa(x)-phi*Qrpad(x)*mu*sqrt(rs)*((1.d0+z)**(-1.d0/3.d0)
     $        -(1.d0-z)**(-1.d0/3.d0))/3.d0+
     $        a1z*mu**3+a2z*mu**4+a3z*mu**5+
     $        a4z*mu**6+a5z*mu**8)/((1.d0+b0**2*mu**2)**4)

         vclrup=eclr-rs/3.d0*eclrrs-(z-1.d0)*eclrz
         vclrdown=eclr-rs/3.d0*eclrrs-(z+1.d0)*eclrz
      endif
      return
      end


      double precision function g0(x)
ccc on-top pair-distribution function
ccc Gori-Giorgi and Perdew, PRB 64, 155102 (2001)
ccc x -> rs
      implicit none
      double precision C0f,D0f,E0f,F0f,x
      C0f             = 0.0819306d0    
      D0f             = 0.752411d0     
      E0f             = -0.0127713d0   
      F0f             = 0.00185898d0   
      g0=(1.d0-(0.7317d0-D0f)*x+C0f*x**2+E0f*x**3+
     $     F0f*x**4)*exp(-abs(D0f)*x)/2.d0
      return
      end

      double precision function g0d(rs)
ccc derivative of on-top pair-distribution function
ccc Gori-Giorgi and Perdew, PRB 64, 155102 (2001)
      implicit none
      double precision Bg0,Cg0,Dg0,Eg0,Fg0,rs
      Cg0             = 0.0819306d0    
      Fg0             = 0.752411d0     
      Dg0             = -0.0127713d0   
      Eg0             = 0.00185898d0
      Bg0             =0.7317d0-Fg0
      g0d=(-Bg0+2*Cg0*rs+3*Dg0*rs**2+4*Eg0*rs**3)/2.d0*exp(-Fg0*rs)
     -   - (Fg0*(1 - Bg0*rs + Cg0*rs**2 + Dg0*rs**3 + Eg0*rs**4))/
     -   2.d0*exp(-Fg0*rs)
      return
      end


      double precision function dpol(rs)
      implicit none
      double precision cf,pi,rs,p2p,p3p
      pi=dacos(-1.d0)
      cf=(9.d0*pi/4.d0)**(1.d0/3.d0)
      p2p    = 0.04d0 
      p3p    = 0.4319d0   
      dpol=2.d0**(5.d0/3.d0)/5.d0*cf**2/rs**2*(1.d0+(p3p-0.454555d0)*rs)
     $     /(1.d0+p3p*rs+p2p*rs**2)
      return
      end

      double precision function dpold(rs)
      implicit none
      double precision cf,pi,rs,p2p,p3p
      pi=dacos(-1.d0)
      cf=(9.d0*pi/4.d0)**(1.d0/3.d0)
      p2p    = 0.04d0 
      p3p    = 0.4319d0   
      dpold=2.d0**(5.d0/3.d0)/5.d0*cf**2*
     - (-2. + (0.454555 - 4.*p3p)*rs + 
     -    (-4.*p2p + 
     -       (0.90911 - 2.*p3p)*p3p)*rs**2
     -      + p2p*(1.363665 - 3.*p3p)*
     -     rs**3)/
     -  (rs**3*(1. + p3p*rs + p2p*rs**2)**2)
      return
      end

      double precision function Qrpa(x)
      implicit none
      double precision pi,a2,b2,c2,d2,x,Acoul
      pi=dacos(-1.d0)
      Acoul=2.d0*(log(2.d0)-1.d0)/pi**2
      a2              = 5.84605d0 
      c2              = 3.91744d0 
      d2              = 3.44851d0
      b2=d2-3.d0/(2.d0*pi*Acoul)*(4.d0/(9.d0*pi))**(1.d0/3.d0)
      Qrpa=Acoul*log((1.d0+a2*x+b2*x**2+c2*x**3)/(1.d0+a2*x+d2*x**2))
      return
      end

      double precision function Qrpad(x)
      implicit none
      double precision pi,a2,b2,c2,d2,x,Acoul
      pi=dacos(-1.d0)
      Acoul=2.d0*(log(2.d0)-1.d0)/pi**2
      a2              = 5.84605d0 
      c2              = 3.91744d0 
      d2              = 3.44851d0
      b2=d2-3.d0/(2.d0*pi*Acoul)*(4.d0/(9.d0*pi))**(1.d0/3.d0)
      Qrpad=Acoul*((x*(b2*(2.d0 + a2*x) + 
     -      c2*x*(3.d0 + 2.d0*a2*x) + 
     -      d2*(-2.d0 - a2*x + c2*x**3)))/
     -  ((1.d0 + a2*x + d2*x**2)*
     -    (1.d0 + a2*x + b2*x**2 + c2*x**3)))
      return
      end

      subroutine exchangelr(rs,z,mu,exlr)
ccc Hartree atomic units used
ccc for given density parameter rs, spin polarization z
ccc and cutoff mu it gives the exchange energy of the LR gas
ccc  => exlr
      implicit none
      double precision rs,z,mu,exlr
      double precision pi,alpha,fx,y
      pi=dacos(-1.d0)
      alpha=(4.d0/9.d0/pi)**(1.d0/3.d0)
      if(abs(z).eq.1.d0) then
         y=mu*alpha*rs/2.d0/2.d0**(1.d0/3.d0)
         fx=-((-3*y + 4*y**3 +(2*y - 4*y**3)*exp(-1./(4.*y**2)) + 
     $        sqrt(pi)*erf(1/(2.*y)))/pi)
         exlr=mu*fx
      else
         y=mu*alpha*rs/2.d0/(1.+z)**(1.d0/3.d0)
         fx=-((-3*y + 4*y**3 +(2*y - 4*y**3)*exp(-1./(4.*y**2)) + 
     $        sqrt(pi)*erf(1/(2.*y)))/pi)
         exlr=(1.d0+z)*mu*fx/2.d0
         y=mu*alpha*rs/2.d0/(1.-z)**(1.d0/3.d0)
         fx=-((-3*y + 4*y**3 +(2*y - 4*y**3)*exp(-1./(4.*y**2)) + 
     $        sqrt(pi)*erf(1/(2.*y)))/pi)
         exlr=exlr+(1.d0-z)*mu*fx/2.d0
      endif
      return
      end

      subroutine vexchangelr(rs,z,mu,vxlrup,vxlrdown)
ccc Hartree atomic units used
ccc for given density parameter rs, spin polarization z
ccc and cutoff mu it gives the exchange LSD potential for LR interaction
ccc  => vxlrup (spin-up electrons), vxlrdown (spin-down electrons)
      implicit none
      double precision rs,z,mu,vxlrup,vxlrdown
      double precision pi,alpha,fx,fx1,y,exlr,derrs,derz
      pi=dacos(-1.d0)
      alpha=(4.d0/9.d0/pi)**(1.d0/3.d0)
      if(z.eq.1.d0) then
         vxlrup=(rs*alpha*mu**2)/
     -   (2**(1.d0/3.d0)*pi) - (rs*alpha*mu**2)/(2**(1.d0/3.d0)*pi)*
     -     exp(-2**(2.d0/3.d0)/(rs**2*alpha**2*mu**2)) - 
     -  (mu*erf(2**(1.d0/3.d0)/(rs*alpha*mu)))/sqrt(Pi)
         vxlrdown=0.d0
      elseif(z.eq.-1.d0) then
         vxlrdown=(rs*alpha*mu**2)/
     -   (2**(1.d0/3.d0)*pi) - (rs*alpha*mu**2)/(2**(1.d0/3.d0)*pi)*
     -     exp(-2**(2.d0/3.d0)/(rs**2*alpha**2*mu**2)) - 
     -  (mu*erf(2**(1.d0/3.d0)/(rs*alpha*mu)))/sqrt(Pi)
         vxlrup=0.d0
      else       
         y=mu*alpha*rs/2.d0/(1.+z)**(1.d0/3.d0)
         fx=-((-3*y + 4*y**3 +(2*y - 4*y**3)*exp(-1./(4.*y**2)) + 
     $        sqrt(pi)*erf(1/(2.*y)))/pi)
         fx1=(3.d0*(1 + (-4.d0 + 4.d0*exp(-1.d0/(4.d0*y**2)))*y**2))/pi
         derrs=1.d0/4.d0*(1.d0+z)**(2.d0/3.d0)*mu**2*alpha*fx1
         derz=1.d0/2.d0*mu*fx-1.d0/6.d0*fx1*mu*y
         vxlrup=rs/3.d0*derrs+(z-1.d0)*derz
         vxlrdown=rs/3.d0*derrs+(z+1.d0)*derz
         
         y=mu*alpha*rs/2.d0/(1.-z)**(1.d0/3.d0)
         fx=-((-3*y + 4*y**3 +(2*y - 4*y**3)*exp(-1./(4.*y**2)) + 
     $        sqrt(pi)*erf(1/(2.*y)))/pi)
         fx1=(3.d0*(1 + (-4.d0 + 4.d0*exp(-1.d0/(4.d0*y**2)))*y**2))/pi
         derrs=1.d0/4.d0*(1.d0-z)**(2.d0/3.d0)*mu**2*alpha*fx1
         derz=-1.d0/2.d0*mu*fx+1.d0/6.d0*fx1*mu*y
         vxlrup=vxlrup+rs/3.d0*derrs+(z-1.d0)*derz
         vxlrdown=vxlrdown+rs/3.d0*derrs+(z+1.d0)*derz
      
         call exchangelr(rs,z,mu,exlr)
         vxlrup=exlr-vxlrup
         vxlrdown=exlr-vxlrdown
      endif
      return
      end

