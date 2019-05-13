! <compile=optimized>
!  ParameterReader handles user-defined parameters for MNDO-type calculations
!  Taisung Lee (Rutgers, 2011)

#include "copyright.h"
#include "dprec.fh"

module ParameterReader
  
  implicit none  
  
  type ParameterEntry
     character(8) :: name='Empty'
     character(2) :: elementName
     integer :: atomicNumber=-1
     _REAL_ :: value=0.0D0
  end type ParameterEntry

  public ParameterEntry, ReadParameterFile, GetNumberParameterEntries, GetParameterEntry
  


    
contains
  
  subroutine ReadParameterFile(fileName, parameterEntries,ParameterFileExisting)

    use UtilitiesModule, only : Upcase
    use ElementOrbitalIndex, only : GetAtomicNumber

    implicit none

    type(ParameterEntry), allocatable, dimension(:), intent(inout) :: parameterEntries
    character(len=*), intent(in) :: fileName
    
    integer,parameter :: readUnit=21
    
    integer :: ios, counter, i, j, atomicNumber
    character(len=80) :: buffer
    type(ParameterEntry), allocatable, dimension(:) :: tempEntries
    integer :: fsize                    ! to estimate file size
    
    logical, intent(inout) :: ParameterFileExisting
    
    ParameterFileExisting=.false.
    open(unit=readUnit,file=fileName,iostat=ios)

    if (ios==0) then
    
       rewind(readUnit)

       ! Esimate the largest possible number of entries
       ! inquire(unit=readUnit, size=fsize)  ! This is F2003 syntax
       ! fstatus=FSTAT(readUnit, tempbuf)    ! This is a non-standard extension
       ! fsize=min(2000,tempbuf(8))/20
       fsize = 1000
       allocate(tempEntries(fsize))
 
       counter=0       
       do while (ios==0) 
          counter=counter+1
            
          ! read a whole line and take out the first three data
          read(unit=readUnit,fmt='(A80)',iostat=ios) buffer  ! read a whole line
          read(buffer,*,iostat=ios) tempEntries(counter)%name, tempEntries(counter)%elementName, tempEntries(counter)%value

          tempEntries(counter)%name=upcase(tempEntries(counter)%name)
          atomicNumber=GetAtomicNumber(tempEntries(counter)%elementName)            

          if ( atomicNumber > 0 ) then
             tempEntries(counter)%atomicNumber=atomicNumber
          else
             counter=counter-1
          end if
          if (tempEntries(counter)%name .eq. 'END') then
             counter=counter-1
             exit
          end if
       end do
        
       if (counter>0) then          
          allocate(parameterEntries(counter))
          parameterEntries(1:counter)=tempEntries(1:counter)
          ParameterFileExisting=.true.
       end if
        
       deallocate(tempEntries)
        
    end if
    close(readUnit) 
   
  end subroutine ReadParameterFile

  function GetNumberParameterEntries(parameterEntries,ParameterFileExisting) result (numberOfEntries)

    implicit none
  
    type(ParameterEntry), dimension(:),intent(in) :: parameterEntries
    logical, intent(in) :: ParameterFileExisting
    integer :: numberOfEntries
    
    numberOfEntries=0

    if(ParameterFileExisting) numberOfEntries=size(parameterEntries)

  end function GetNumberParameterEntries


  function GetParameterEntry(parameterEntries,ParameterFileExisting, index)  result (entry)
    ! this function is to prevent user modifying 
    ! the stored parameterEntries--i.e., user
    ! can only get a copy of any entry

    implicit none
    
    type(ParameterEntry), dimension(:) :: parameterEntries
    logical, intent(in) :: ParameterFileExisting
    integer,intent(in) :: index
    type(ParameterEntry) :: entry

    if (index<=GetNumberParameterEntries(parameterEntries,ParameterFileExisting) .and. index>0) then
       entry=parameterEntries(index)
    end if
    
  end function GetParameterEntry

end module ParameterReader
    
