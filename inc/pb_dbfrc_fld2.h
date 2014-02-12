
      ! inner loop over atoms...
      ! compute reaction field energy, forces, and coulomb field
      ! between the current surface charge and natom atomic charges

      call get_coulomb(natom,crg,crd,acrg,ax,ay,az,qex,qey,qez,eelrf,coulomb)

      ! now it is time to compute dbf
      ! part a: compute E_in and E_out ...
      !    compute E on the inner side of surface position crd(1:3)

      call gradu(xm,ym,zm,-ONE,8,4,fx0,fy0,fz0,up,dudxi0,dudyi0,dudzi0,phi,zv)

      dudxi0 = -dudxi0*FOURPI*eps0*AMBER_ELECTROSTATIC
      dudyi0 = -dudyi0*FOURPI*eps0*AMBER_ELECTROSTATIC
      dudzi0 = -dudzi0*FOURPI*eps0*AMBER_ELECTROSTATIC

      !    add the coulomb field to get the total E of inner side
      !    convert to displacement, D, for consistency with other methods

      dudxi0 = (dudxi0 + coulomb(1))*epsin/eps0
      dudyi0 = (dudyi0 + coulomb(2))*epsin/eps0
      dudzi0 = (dudzi0 + coulomb(3))*epsin/eps0

      !    get the normal direction of position crd(1:3)
      !    get normal field component, which is continuous
      !    get D of the outer side based on the jump conditions

      rn = rn/sqrt( sum(rn*rn) )

      dudni = dudxi0*rn(1)+dudyi0*rn(2)+dudzi0*rn(3)
      !dudno = dudni
      !dudxo0 = (dudxi0 - dudni*rn(1))/epsin*epsout + dudno*rn(1)
      !dudyo0 = (dudyi0 - dudni*rn(2))/epsin*epsout + dudno*rn(2)
      !dudzo0 = (dudzi0 - dudni*rn(3))/epsin*epsout + dudno*rn(3)

      ! part b: apply the normal field approximation
      !         or use the total field

      E2 = dudni*dudni
      !E2 = dudxi0*dudxo0 + dudyi0*dudyo0 + dudzi0*dudzo0

      ! part c: compute the surface force element

      fbnd = HALF*INV_FOURPI*(eps0*(epsin-epsout)/(epsin*epsout))*E2*ds
      dum(1) = fbnd*rn(1)
      dum(2) = fbnd*rn(2)
      dum(3) = fbnd*rn(3)

      call dbfrc(natom,iatm,x,dum,cnx,cny,cnz,rnx,rny,rnz)
