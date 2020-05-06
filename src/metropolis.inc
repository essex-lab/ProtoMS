
c     This is the metropolis test common block - used to 
c     add variables that may be needed for any metropolis
c     test. Its a bit messy, but it allows changes to
c     occur to the function without changing the way it
c     is called...

c     old and new volumes for NPT metropolis test (movtyp = 1)
      double precision mt_vold,mt_vnew
c     1 / RT
      double precision mt_beta_T
c     Pcav
      double precision Pcav
c     Volume
      double precision V
c     biasing probability
      double precision mt_bias
c      for GCMC test
      logical gcmc   

     
      common /METTESTDATA/ V,Pcav,mt_vold,mt_vnew,mt_beta_T,mt_bias,gcmc
                      save /METTESTDATA/
