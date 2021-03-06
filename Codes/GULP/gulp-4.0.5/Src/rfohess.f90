  subroutine rfohess(hessian,nvar,nmin)
!
!  Form hessian from symmetry reduced second derivative matrix
!
!   6/05 Intent added
!
!  Conditions of use:
!
!  GULP is available free of charge to academic institutions
!  and non-commerical establishments only. Copies should be
!  obtained from the author only and should not be distributed
!  in any form by the user to a third party without the express
!  permission of the author. This notice applies to all parts
!  of the program, except any library routines which are
!  distributed with the code for completeness. All rights for
!  such routines remain with the original distributor.
!
!  No claim is made that this program is free from errors and
!  no liability will be accepted for any loss or damage that
!  may result. The user is responsible for checking the validity
!  of their results.
!
!  Copyright Curtin University 2005
!
!  Julian Gale, NRI, Curtin University, June 2005
!
  use derivatives
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)  :: nmin
  integer(i4), intent(in)  :: nvar
  real(dp),    intent(out) :: hessian(*)
!
!  Local variables
!
  integer(i4) :: i
  integer(i4) :: ind
  integer(i4) :: j
!
!  Transfer hessian to correct matrix
!
  ind = 0
  do i = nmin,nvar
    do j = nmin,i
      ind = ind + 1
      hessian(ind) = derv2(j,i)
    enddo
  enddo
!
  return
  end
