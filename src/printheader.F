      subroutine printHeader
      implicit none
      include 'printer.inc'
      include 'constants.inc'
      include 'version.inc'
c###############################################################
c
c     This routine prints the header. It includes
c     'idme.inc', which is an automatically generated file
c     that enables the version and compile date to be
c     known by the program.
c
c     (C) Christopher Woods, October 2004
c
c###############################################################
      
      call printBoxBar(HEADER)
      call printBoxText(HEADER,"")
      write(printstring,*) progname," ",progversion
      call printBoxText(HEADER,printstring)
      call printBoxText(HEADER,"")
      call printBoxText(HEADER,"Prototype Molecular Simulation")
      call printBoxText(HEADER,"")
      write(printstring,*) "This is a program to perform simulations on molecular systems."
      call printBoxText(HEADER,printstring) 
      call printBoxText(HEADER,"")
      
      write(printstring,*) "This program is free software; you can redistribute ",
     .              "it and/or modify it under the terms of the GNU Public License ",
     .              "as published by the Free Software Foundation; either version 2 of the ",
     .              "License, or (at your option) any later version."
      call printBoxText(HEADER,printstring)
      call printBoxText(HEADER,"")
      
      write(printstring,*) "This program is distributed in the hope that it will be useful, ",
     .              "but without any warranty; without even the implied warranty of ",
     .              "merchantability or fitness for a particular purpose.  See the GNU General ",
     .              "Public License for more details."
      call printBoxText(HEADER,printstring)
      call printBoxText(HEADER,"")

      write(printstring,*) "You should have received a copy of the GNU General Public License along ",
     .              "with this program; if not, write to the Free Software Foundation, Inc., ",
     .              "59 Temple Place, Suite 330, Boston, MA  02111-1307  USA."
      call printBoxText(HEADER,printstring)
      call printBoxText(HEADER,"")
      
c     include dynamically generated compilation information
      include 'idme.inc'      
            
      call printBoxText(HEADER,"")

      call printBoxText(HEADER,"Using physical constants / conversion factors...")
      write(printstring,*) "PI = ",PI
      call printBoxText(HEADER,printstring)
      write(printstring,*) "Gas Constant (kcal mol-1 K-1) = ",GASR
      call printBoxText(HEADER,printstring)
      write(printstring,*) "1/4*Pi*Epsilon0 (UnitCharge-2 kcal A mol-1) = ",FAC4PIEPS0
      call printBoxText(HEADER,printstring)
      write(printstring,*) "Atmospheric pressure (kcal mol-1 A-3) = ",P2INTERNALP
      call printBoxText(HEADER,printstring)
      call printBoxText(HEADER,"")

      call printBoxText(HEADER,"(C) Samuel Genheden, Gregory A. Ross, Richard Bradshaw,",
     .    " Chris Cave-Ayland, Ana I. Cabedo Martinez 2014-2015")
      call printBoxText(HEADER,"(C) Christopher Woods, Julien Michel, Samuel Genheden 2002-2014")
      call printBoxText(HEADER,"(C) JAWS and GCMC protocols by Michael Bodnarchuk, 2012")
      call printBoxText(HEADER,"(C) OpenMP/MPI RETI code by Samuel Genheden, 2012-2014")
      call printBoxText(HEADER,"")
      
      call printBoxBar(HEADER)
      
      end
