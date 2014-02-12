! These are parameters that allow more than 7 residue ranges to be specified
! on a single line in a card in the input file.  This limit of 7 was hardcoded
! in, though all references to that have been replaced with either MAX_RES or
! MAX_RES_FIELDS as needed. Each field is a range, so it requires storage for 2
! numbers (hence, MAX_RES_FIELDS = MAX_RES * 2), and the first element is a string
! that reads RES or RRES or LRES or something. The last element generally had to
! be a 0, or rfree would get confused while parsing, so I added an extra 2 values
! onto MAX_RES_FIELDS.  However, if MAX_RES is used as the end index in the do
! loop in rfree.f, this should no longer be necessary, but I will keep it just in
! case. Note that there are still some hard limits on the size of items, such as the
! 137-byte size of the line buffer in rfree.f.  Also, not all hard limits in rgroup
! were changed, just the ones I thought would be used the most. -- JMS, 11/2010

integer, parameter :: MAX_RES = 20
integer, parameter :: MAX_RES_FIELDS = 2 * MAX_RES + 2
