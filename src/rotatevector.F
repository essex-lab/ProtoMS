      subroutine rotateVector(newvec,vec,theta)
      
C*******This code has been converted from a piece of
C*******C code, placed under the GPL. I am unable to 
C*******find the name of the original author, though
C*******if you recognise it, please email me your name
C*******and I will be happy to acknowledge you.
      
      double precision newvec(3),vec(3)
      double precision t1,t2,t3,t4,t5,t6,t7,t8,t9,t10
      double precision t11,t12,t13,t14,t15,t16,t17,t18,t19
      double precision ia,ib,ic,id,ie,ff,ig,ih,ii
      double precision theta
      double precision x1,y1,z1,x2,y2,z2
      double precision sintheta,costheta
c##############################################

      sintheta = dsin(theta)
      costheta = dcos(theta)
      
      x1 = newvec(1)
      y1 = newvec(2)
      z1 = newvec(3)
      
      x2 = vec(1)
      y2 = vec(2)
      z2 = vec(3)
      
      t1 = x2*x2
      t2 = y2*y2
      t3 = z2*z2
      
      t4 = 1.0-t1
      t5 = 1.0-t2
      t6 = 1.0-t3
      
      t7 = 1.0-costheta
      
      t8 = x2*y2
      t9 = x2*z2
      t10 = y2*z2
      
      t11 = x2*sintheta
      t12 = y2*sintheta
      t13 = z2*sintheta
      
      t14 = t4*costheta
      t15 = t5*costheta
      t16 = t6*costheta
      
      t17 = t8*t7
      t18 = t9*t7
      t19 = t10*t7
      
      ia = t1+t14
      ib = t17+t13
      ic = t18-t12
      
      id = t17-t13
      ie = t2+t15
      ff = t19+t11
      
      ig = t18+t12
      ih = t19-t11
      ii = t3+t16
      
      newvec(1) = x1*ia + y1*ib + z1*ic
      newvec(2) = x1*id + y1*ie + z1*ff
      newvec(3) = x1*ig + y1*ih + z1*ii
      
      end
