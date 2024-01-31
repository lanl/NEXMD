module multitmd

#include "assert.fh"
#include "copyright.h"
#include "dprec.fh"

    public mtmdread, mtmdenergy, mtmdprint

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine mtmdread(refcrds, mtmdrmsdarr, mtmdfrc, mtmdmaskarr, mtmdsteps, mtmdtgtatms, &
        x, name, irsnam, ipres, isymbl, natom, nres, tgtnum, maxtgt, mtmdlun, tgtin, reflun)

        ! Subroutine Multiple Targeted Molecular Dynamics READ.

        ! This subroutine reads 1) The mtmd file containing the tgt namelist
        !                       2) All the reference files identifed in that file

        ! Author: Matthew Seetin
        ! Date: 5/2008

        ! INPUT VARIABLES:
        ! ---------------

        !        X(*) : Coordinate Array
        !     NAME(*) : The name corresponding to atom I.
        !    IPRES(*) : Residue pointer array. IPRES(I) points to the first atom
        !               of residue I.
        !       TGTIN : Unit from which the restraints/weight changes will be read.
        !        IOUT : Unit for informational prints.
        !       NATOM : The number of atoms.
        !        NRES : The number of residues.
        !      MAXTGT : The maximum number of restraints for which array space is avail.
        !       TGTIN : The file name of the MTMD file given in the command line
        !     MTMDLUN : The logical unit number to be used for opening the MTMD file
        !      REFLUN : The logical unit number to be used for opening the reference files.

        ! OUTPUT:

        !  MTMDRMSDARR(2,I) : The slope and intercept, respectively, defining the dependency
        !                 the target RMSD relative to reference I.
        !   MTMDFRC(2,I): The slope and intercept for the force constant to reference I.
        !     REFCRDS(*): The coordinates of all reference structures.  The coordinates
        !                 for each structure are listed one after another in this array
        ! MTMDMASKARR(*):  An integer array identifying which atoms to use for fitting
        !                 the RMSD calculation, and where forces will be applied.  The
        !                 mask for each reference structure are listed one after another.
        ! MTMDSTEPS(3,I): The beginning (1) and ending (2) steps over which this
        !                 restraint is to be applied. MTMDSTEPS(3,I) gives the increment
        !                 between changes in the target values. If MTMDSTEPS(3,I)=0,
        !                 the target values are varied continuously.
        ! MTMDTGTATMS(I): The number of atoms in the mask for each reference.
        !        TGTNUM : The total number of reference structures defined.

        use findmask, only : atommask
        use file_io_dat, only : MAX_FN_LEN

        implicit none

        character(len=4) name(*), isymbl(*), irsnam(*)
        _REAL_ x
        integer ipres, mtmdtgtatms
        dimension x(*), ipres(*), mtmdtgtatms(*)
        integer natom
        integer nres
        integer tgtnum
        integer maxtgt
        integer mtmdlun
        character(len=MAX_FN_LEN) tgtin
        integer reflun

        ! Contents of tgt namelist
        character(MAX_FN_LEN) refin
        character(256) mtmdmask
        integer mtmdstep1
        integer mtmdstep2
        _REAL_ mtmdforce
        _REAL_ mtmdforce2
        _REAL_ mtmdrmsd
        _REAL_ mtmdrmsd2
        integer mtmdform
        integer mtmdninc
        integer mtmdmult
        integer mtmdvari

        _REAL_ refcrds
        dimension refcrds(*)

        _REAL_ mtmdrmsdarr
        dimension mtmdrmsdarr(2, *)

        _REAL_ mtmdfrc
        dimension mtmdfrc(2, *)

        integer mtmdmaskarr, mtmdsteps
        dimension mtmdmaskarr(*), mtmdsteps(3, *)

        integer ifind
        integer step1, step2
        integer i, j
        integer crd_ptr
        integer mask_ptr
        integer ntu
        _REAL_ rntu

        namelist /tgt/ refin, mtmdmask, mtmdstep1, mtmdstep2, mtmdforce, mtmdforce2, &
            mtmdrmsd, mtmdrmsd2, mtmdform, mtmdninc, mtmdmult, mtmdvari

        ! Defaults for the following variables will be the value given in the previous namelist
        ! If they have not been defined in any namelist, the defaults are as given below.

        call amopen(mtmdlun, tgtin, 'O', 'F', 'R')

        mtmdmask = '*'
        mtmdstep1 = 0
        mtmdstep2 = 0
        mtmdforce = 0.0d0
        mtmdforce2 = 0.0d0
        mtmdrmsd = 0.0d0
        mtmdrmsd2 = 0.0d0
        mtmdform = 1
        mtmdninc = 0
        mtmdvari = 0
        mtmdmult = 0

        do i = 1, maxtgt ! loop over all targets
            refin = ''

            call nmlsrc('tgt', mtmdlun, ifind)
            if (ifind == 0) exit ! exit from loop over all targets if none are found.

            read (mtmdlun, nml=tgt, end=100)

100         if (len_trim(refin) <= 0) then
                tgtnum = i - 1
                exit ! exit from loop over all targets if no reference file is given.
            end if

            crd_ptr = (i - 1)*3*natom + 1
            mask_ptr = (i - 1)*natom + 1

            ! Most reference files are formatted inpcrds.  If not, the user can specify MTMDFORM=0
            ! in the mtmd file to use an unformatted (binary) reference file.
            ! DRR - File opening/closing is now done in rdrest. reflun is set to 10
            !       inside rdrest (although it should probably be defined in files.h).

            call rdrest(natom, mtmdform, refin, refcrds(crd_ptr))

            write (6, '(a,a)') '     Read in coords from ', refin(1:len_trim(refin))

            ! Parse the mask to find the atom numbers of the atoms that make it up

            call atommask(natom, nres, 0, name, isymbl, &
                ipres, irsnam, x, mtmdmask, mtmdmaskarr(mask_ptr))

            mtmdtgtatms(i) = 0

            do j = 1, natom
                if (mtmdmaskarr(mask_ptr - 1 + j) <= 0) cycle
                mtmdtgtatms(i) = mtmdtgtatms(i) + 1
                mtmdmaskarr(mask_ptr - 1 + mtmdtgtatms(i)) = j
            end do

            write (6, '(a,a,a,i5,a)') &
                '     Mask "', mtmdmask(1:len_trim(mtmdmask) - 1), &
                '" matches ', mtmdtgtatms(i), ' atoms'

            ! Read step information into mtmdsteps

            if (mtmdstep2 > mtmdstep1 .or. mtmdstep2 <= 0) then
                mtmdsteps(1, i) = mtmdstep1
                mtmdsteps(2, i) = mtmdstep2
            else
                mtmdsteps(1, i) = mtmdstep2
                mtmdsteps(2, i) = mtmdstep1
            end if
            mtmdsteps(3, i) = max(1, mtmdninc)

            ! read in steps, time evolution, etc.  Calculate slopes and intercepts.

            if (mtmdstep2 == 0 .or. mtmdvari <= 0 .or. mtmdstep1 == mtmdstep2) then
                mtmdrmsdarr(1, i) = 0.0d0
                mtmdrmsdarr(2, i) = mtmdrmsd
                mtmdfrc(1, i) = 0.0d0
                mtmdfrc(2, i) = mtmdforce
            else
                step1 = mtmdstep1
                step2 = mtmdstep2
                mtmdrmsdarr(1, i) = (mtmdrmsd2 - mtmdrmsd)/(step2 - step1)
                mtmdrmsdarr(2, i) = (mtmdrmsd2) - mtmdrmsdarr(1, i)*step2

                ! If MTMDMULT > 0, then use a constant multipier every MTMDNINC steps, instead
                ! of a linear interpolation for the force constants.

                if (mtmdmult > 0) then
                    ntu = (mtmdsteps(2, i) - mtmdsteps(1, i))/mtmdsteps(3, i)
                    rntu = max(1, ntu)
                    mtmdsteps(3, i) = -mtmdsteps(3, i)
                    mtmdfrc(1, i) = (mtmdforce2/mtmdforce)**(1.0d0/rntu)
                    mtmdfrc(2, i) = mtmdforce
                else
                    mtmdfrc(1, i) = (mtmdforce2 - mtmdforce)/(step2 - step1)
                    mtmdfrc(2, i) = mtmdforce2 - mtmdfrc(1, i)*step2
                end if
            end if

        end do ! close loop over all targets

        return

    end subroutine mtmdread

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine mtmdenergy(enemtmd, refcrds, mtmdrmsdarr, mtmdfrc, mtmdmaskarr, &
        mtmdsteps, mtmdtgtatms, x, f, tgtrmsdcurr, rmsdcurr, name, irsnam, &
        mass, natom, nres, tgtnum, maxtgt, nstep)

        ! Subroutine Multiple Targeted Molecular Dynamics Energy.

        ! This subroutine       1) Calculates current target RMSD, and force constant
        !                       2) Calculates current RMSD to each reference
        !                       3) Calculates associated energies and forces for each reference

        ! Author: Matthew Seetin
        ! Date: 6/2008

        ! INPUT VARIABLES:
        ! ---------------

        !          X(*) : Coordinate Array
        !       NAME(*) : The name corresponding to atom I.
        !         NATOM : The number of atoms.
        !          NRES : The number of residues.
        !        MAXTGT : The maximum number of restraints for which array space is avail.
        !         TGTIN : The file name of the MTMD file given in the command line
        !  MTMDRMSDARR(2,I) : The slope and intercept, respectively, defining the dependency
        !                 the target RMSD relative to reference I.
        !   MTMDFRC(2,I): The slope and intercept for the force constant to reference I.
        !     REFCRDS(*): The coordinates of all reference structures.  The coordinates
        !                 for each structure are listed one after another in this array
        ! MTMDMASKARR(*):  An integer array identifying which atoms to use for fitting
        !                 the RMSD calculation, and where forces will be applied.  The
        !                 mask for each reference structure are listed one after another.
        ! MTMDSTEPS(3,I): The beginning (1) and ending (2) steps over which this
        !                 restraint is to be applied. MTMDSTEPS(3,I) gives the increment
        !                 between changes in the target values. If MTMDSTEPS(3,I)=0,
        !                 the target values are varied continuously.
        ! MTMDTGTATMS(I): The number of atoms in the mask for each reference.
        !        TGTNUM : The total number of reference structures defined.
        !         NSTEP : The current step of the calculation

        ! OUTPUT VARIABLES:
        ! ---------------
        !
        !       ENEMTMD : The summed energy of the RMSD restrain from each reference
        !          F(*) : Force Array
        !   RMSDCURR(I) : The current RMSD between the molecule and each reference
        ! TGTRMSDCURR(I): The current target RMSD to each reference

        implicit none

        _REAL_ enemtmd
        _REAL_ refcrds, mtmdrmsdarr, mtmdfrc, x, f, tgtrmsdcurr, rmsdcurr, mass
        dimension refcrds(*), mtmdrmsdarr(2, *), mtmdfrc(2, *), x(*), f(*), &
            tgtrmsdcurr(*), rmsdcurr(*), mass(*)

        character(len=4) name, irsnam
        dimension name(*), irsnam(*)

        integer natom, nres, tgtnum, maxtgt, nstep

        integer mtmdsteps, mtmdtgtatms, mtmdmaskarr
        dimension mtmdsteps(3, *), mtmdtgtatms(*), mtmdmaskarr(*)

        _REAL_ currforce
        integer crd_ptr, mask_ptr
        integer i, nstepu
        _REAL_ enemtmdcurr
        logical rmsok

        enemtmd = 0.0d0

        do i = 1, tgtnum

            rmsdcurr(i) = 0.0d0
            if (nstep < mtmdsteps(1, i) .or. (nstep > mtmdsteps(2, i) .and. &
                mtmdsteps(2, i) > 0)) cycle

            ! Calculate the values of tgtrmsdcurr and currforce for this step & restraint:

            ! Vary the step used in calculating the values only every ABS(MTMDSTEPS(3,I)) steps.
            ! If the user did not specify MTMDNINC in the input, ABS(MTMDSTEPS(3,I))=1.

            nstepu = nstep - mod(nstep - mtmdsteps(1, i), abs(mtmdsteps(3, i)))
            tgtrmsdcurr(i) = mtmdrmsdarr(1, i)*nstepu + mtmdrmsdarr(2, i)

            ! If mtmdsteps(3,I) > 0, then the weights are linearly interpolated.
            ! If mtmdsteps(3,I) < 0, then the weights are modified by a multplicative factor.

            if (mtmdsteps(3, i) > 0) then
                currforce = mtmdfrc(1, i)*nstepu + mtmdfrc(2, i)
            else
                nstepu = (nstep - mtmdsteps(1, i))/abs(mtmdsteps(3, i))
                currforce = mtmdfrc(2, i)*mtmdfrc(1, i)**nstepu
            end if

            crd_ptr = (i - 1)*3*natom + 1
            mask_ptr = (i - 1)*natom + 1

            ! Calculate the current RMSD to reference i.

            call rmsfit(refcrds(crd_ptr), x, mass, mtmdmaskarr(mask_ptr), &
                mtmdmaskarr(mask_ptr), rmsdcurr(i), mtmdtgtatms(i), mtmdtgtatms(i), rmsok)

            if (.not. rmsok) then
                write (6, *) 'Fatal Error calculating RMSD !'
                call mexit(6, 1)
            end if

            ! Calculate t he forces on each atom in the mask and the energy of the restraint.

            call xtgtmd(enemtmdcurr, mtmdmaskarr(mask_ptr), x, f, refcrds(crd_ptr), mass, tgtrmsdcurr(i), currforce, rmsdcurr(i), mtmdtgtatms(i))

            enemtmd = enemtmd + enemtmdcurr
        end do

    end subroutine mtmdenergy

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine mtmdprint(tgtnum, rmsdcurr, tgtrmsdcurr, iout)

        ! Subroutine Multiple Targeted Molecular Dynamics Print.

        ! This subroutine prints the current RMSD and the target RMSD for each
        ! reference in the out file.

        ! Author: Matthew Seetin
        ! Date: 6/2008

        implicit none

        _REAL_ rmsdcurr, tgtrmsdcurr
        dimension rmsdcurr(*), tgtrmsdcurr(*)

        integer tgtnum, iout
        integer i

        do i = 1, tgtnum
            write (iout, '(a,i3,a,f8.3)') 'Current RMSD from reference ', i, ':      ', rmsdcurr(i)
            write (iout, '(a,i3,a,f8.3)') 'Current target RMSD to reference ', i, ': ', tgtrmsdcurr(i)
        end do

        write (6, 8088)

8088    format(t2, 78('-'),/)

    end subroutine mtmdprint

end module multitmd
