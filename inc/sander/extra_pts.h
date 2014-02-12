! frames are centered on an atom atcenter that "owns" extra points
! numep are number of extra points attached to atcenter
! epframe are pointers to (at most 2) extra points attached to atcenter
! type is 1 if atom has at least 2 other bonds, pref to heavy atoms
! type is 2 if atom has only one other bond e.g. O6 in guanine
! first middle and third are atom nums defining frame
! loc_frame are local coords

#define BC_EXTRA_PT 12

integer ifrtyp,iatcen,inumep,iepfr,ifrst,imid,ithrd
integer leploc
integer numfr,numextra
integer numnb14,inb_14
common/fr_ptr/ifrtyp,iatcen,inumep,iepfr,ifrst,imid,ithrd, &
      leploc,inb_14,numfr,numextra,numnb14
