!+++++++++++++++++++++++++++++++++++++++++
!    Ternary-EQE_EL
!+++++++++++++++++++++++++++++++++++++++++ 
program eigenequation 
implicit doubleprecision(a-h,o-z)
parameter(kk=4)
dimension z(kk,kk), en(kk), eddi(kk)
real(8)      x, kr_LED, knr_LED, kr_CTF, knr_CTF, kr_LENFA, knr_LENFA, kr_CTNFA, knr_CTNFA, kr, knr, EQE_EL
real(8)      Bolz, kbT

common /para/E_LED, E_CTF, E_LENFA, E_CTNFA, t_LED_CTF, t_LED_LENFA, t_LED_CTNFA, t_CTF_CTNFA, t_LENFA_CTNFA

E_LED             =   13600.0d0
E_CTF             =   10000.0d0
E_LENFA           =   11200.0d0
E_CTNFA           =   11200.0d0
t_LED_CTF         =   150.0d0
t_LED_LENFA       =   0.0d0
t_LED_CTNFA       =   0.0d0
t_CTF_CTNFA       =   150.0d0
t_LENFA_CTNFA     =   150.0d0


kr_LED            =   1.0d+8
knr_LED           =   1.0d+10
kr_CTF            =   1.0d+4
knr_CTF           =   1.0d+11
kr_LENFA          =   1.0d+8
knr_LENFA         =   1.0d+10
kr_CTNFA          =   1.0d+4
knr_CTNFA         =   1.0d+11
              
kbT               =   200

 open(22,file='eqe-tenary.dat',status='unknown')
        


        call hami(kk,z)
        call TRED2(kk,kk,EN,EDDI,Z)
        call IMTQL2(kk,kk,EN,EDDI,Z,IERR)
        
!        write(22,100) (z(m,3), m=1,kk)
!100     format(4f18.9) 

Bolz = 0.0d0
do m = 1, kk
        Bolz = Bolz + exp(-en(m)/kbT)
enddo

x = -0.1
do n = 1, 11
   x = x + 0.1   
    
kr     = 0.0d0
knr    = 0.0d0
EQE_EL = 0.0d0
do m = 1, kk
        kr = kr + (exp(-en(m)/kbT)/Bolz)*(z(1,m)*z(1,m)*kr_LED+x*z(2,m)*z(2,m)*kr_CTF+ &
                   & (1-x)*z(3,m)*z(3,m)*kr_LENFA+(1-x)*z(4,m)*z(4,m)*kr_CTNFA)
        knr= knr+ (exp(-en(m)/kbT)/Bolz)*(z(1,m)*z(1,m)*knr_LED+x*z(2,m)*z(2,m)*knr_CTF+ & 
                   & (1-x)*z(3,m)*z(3,m)*knr_LENFA+(1-x)*z(4,m)*z(4,m)*knr_CTNFA)
enddo

        EQE_EL = 0.2*(kr/(kr+knr))
        write(*,*) EQE_EL
enddo

stop
end

!+++++++++subroutine hami+++++++++++++++++
subroutine hami(kk,a)
implicit doubleprecision(a-h,o-z)
dimension a(kk,kk)
common /para/E_LED, E_CTF, E_LENFA, E_CTNFA, t_LED_CTF, t_LED_LENFA, t_LED_CTNFA, t_CTF_CTNFA, t_LENFA_CTNFA

a=0.0

a(1,1)  =  E_LED
a(1,2)  =  t_LED_CTF
a(1,3)  =  t_LED_LENFA
a(1,4)  =  t_LED_CTNFA

a(2,1)  =  t_LED_CTF
a(2,2)  =  E_CTF
a(2,3)  =  0.0d0
a(2,4)  =  t_LENFA_CTNFA

a(3,1)  =  t_LED_LENFA
a(3,2)  =  0.0d0
a(3,3)  =  E_LENFA
a(3,4)  =  t_LENFA_CTNFA

a(4,1)  =  t_LED_CTNFA
a(4,2)  =  t_LENFA_CTNFA
a(4,3)  =  t_LENFA_CTNFA
a(4,4)  =  E_CTNFA

return
end
!++++++++++++++++++++++++++++++++++++++++++	
       SUBROUTINE TRED2(NM,N,D,E,Z)  
       IMPLICIT DOUBLE PRECISION(A-H,O-z)
       DIMENSION D(N),E(N),Z(NM,N)
       IF(N.EQ.1) GOTO 320
       DO 300 II=2,N
         I=N+2-II
         L=I-1
         H=0.0D0
         SCALE=0.0D0
         IF(L.LT.2) GOTO 130
         DO 120 K=1,L
120      SCALE=SCALE+DABS(Z(I,K))
         IF(SCALE.ne.0.0D0) GOTO 140
130      E(I)=Z(I,L)
         GOTO 290
140      DO 150 K=1,L
           Z(I,K)=Z(I,K)/SCALE
           H=H+Z(I,K)**2
150      CONTINUE
         F=Z(I,L)
         G=-DSIGN(DSQRT(H),F)
         E(I)=SCALE*G
         H=H-F*G
         Z(I,L)=F-G
         F=0.0D0
         DO 240 J=1,L
           Z(J,I)=Z(I,J)/(SCALE*H)
           G=0.0D0
           DO 180 K=1,J
180        G=G+Z(J,K)*Z(I,K)
           JP1=J+1
           IF(L.LT.JP1) GOTO 220
           DO 200 K=JP1,L
200        G=G+Z(K,J)*Z(I,K)
220        E(J)=G/H
           F=F+E(J)*Z(I,J)
240      CONTINUE
         HH=F/(H+H)
         DO 260 J=1,L
           F=Z(I,J)
           G=E(J)-HH*F
           E(J)=G
           DO 260 K=1,J
             Z(J,K)=Z(J,K)-F*E(K)-G*Z(I,K)
260      CONTINUE
         DO 280 K=1,L
280      Z(I,K)=SCALE*Z(I,K)
290      D(I)=H
300   CONTINUE
320   d(1)=0.0D0
      E(1)=0.0D0
      DO 500 I=1,N
        L=I-1
        IF(D(I).eq.0.0d0) GOTO 380
        DO 360 J=1,L
          G=0.0D0
          DO 340 K=1,L
340       G=G+Z(I,K)*Z(K,J)
          DO 360 K=1,L
            Z(K,J)=Z(K,J)-G*Z(K,I)
360     CONTINUE
380     D(I)=Z(I,I)
        Z(I,I)=1.0D0
        IF(L.LT.1) GOTO 500
        DO 400 J=1,L
          Z(I,J)=0.0D0
          Z(J,I)=0.0D0
400     CONTINUE
500   CONTINUE
      RETURN
      END
!****************************************************************
      subroutine imtql2(nm,n,d,e,z,ierr)   
      implicit double precision(a-h,o-z)
      dimension d(n),e(n),z(nm,n)
      dmach=2.0d0**(-37)
      ierr=0
      if(n.eq.1)goto 1011
      do 101 i=2,n
 101  e(i-1)=e(i)
      e(n)=0.0d0
      do 240 l=1,n
        j=0
 105    do 110 m=l,n
          if(m.eq.n)goto 120
          if(dabs(e(m)).le.dmach*(dabs(d(m))+dabs(d(m+1))))goto 120
 110    continue
 120    p=d(l)
        if(m.eq.l)goto 240
        if(j.eq.30)goto 1010
        j=j+1
        g=(d(l+1)-p)/(2.0d0*e(l))
        r=dsqrt(g*g+1.0d0)
        g=d(m)-p+e(l)/(g+dsign(r,g))
        s=1.0d0
        c=1.0d0
        p=0.0d0
        mml=m-l
        do 200 ii=1,mml
          i=m-ii
          f=s*e(i)
          b=c*e(i)
          if(dabs(f).lt.dabs(g)) goto 150
          c=g/f
          r=dsqrt(c*c+1.0d0)
          e(i+1)=f*r
          s=1.0d0/r
          c=c*s
          goto 160
 150      s=f/g
          r=dsqrt(s*s+1.0d0)
          e(i+1)=g*r
          c=1.0d0/r
          s=s*c
 160      g=d(i+1)-p
          r=(d(i)-g)*s+2.0d0*c*b
          p=s*r
          d(i+1)=g+p
          g=c*r-b
          do 180 k=1,n
            f=z(k,i+1)
            z(k,i+1)=s*z(k,i)+c*f
            z(k,i)=c*z(k,i)-s*f
 180      continue
 200    continue
        d(l)=d(l)-p
        e(l)=g
        e(m)=0.0d0
        goto 105
 240  continue
      do 300 ii=2,n
        i=ii-1
        k=i
        p=d(I)
        do 260 j=ii,n
          if(d(j).ge.p)goto 260
          k=j
          p=d(j)
 260    continue
        if(k.eq.i) goto 300
        d(k)=d(i)
        d(i)=p
        do 280 ll=1,n
          p=z(ll,i)
          z(ll,i)=z(ll,k)
          z(ll,k)=p
 280    continue
 300  continue
      goto 1011
1010  ierr=l
1011  return
      end
!*****************************
