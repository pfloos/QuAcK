subroutine RGTpp_excitation_density(ispin,nBas,nC,nO,nV,nR,nOO,nVV,ERI,X1,Y1,rho1,X2,Y2,rho2)

! Compute excitation densities for T-matrix self-energy

  implicit none

! Input variables

  integer,intent(in)            :: ispin
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  integer,intent(in)            :: nOO
  integer,intent(in)            :: nVV 
  double precision,intent(in)   :: X1(nVV,nVV)
  double precision,intent(in)   :: Y1(nOO,nVV)
  double precision,intent(in)   :: X2(nVV,nOO)
  double precision,intent(in)   :: Y2(nOO,nOO)

! Local variables

  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  integer                       :: p,q
  integer                       :: ab,cd,ij,kl
  double precision,external     :: Kronecker_delta

! Output variables

  double precision,intent(out)  :: rho1(nBas,nBas,nVV)
  double precision,intent(out)  :: rho2(nBas,nBas,nOO)

  integer                       :: dim_1, dim_2
  double precision, allocatable :: ERI_1(:,:,:)
  double precision, allocatable :: ERI_2(:,:,:)

! Initialization

  rho1(:,:,:) = 0d0
  rho2(:,:,:) = 0d0

!----------------------------------------------
! Singlet manifold
!----------------------------------------------

  if(ispin == 1) then

    !$OMP PARALLEL DEFAULT(NONE)                                         &
    !$OMP          PRIVATE(p, q, a, b, ab, c, d, cd, i, j, ij, k, l, kl) &
    !$OMP          SHARED(nC, nBas, nR, nO, rho1, rho2, ERI, X1, Y1, X2, Y2)
    !$OMP DO COLLAPSE(2)

     do q=nC+1,nBas-nR
        do p=nC+1,nBas-nR
           
           ab = 0
           do a=nO+1,nBas-nR
           do b=a,nBas-nR
              ab = ab + 1
              cd = 0
              do c=nO+1,nBas-nR
                 do d=c,nBas-nR
                    cd = cd + 1
                    rho1(p,q,ab) = rho1(p,q,ab) & 
                                 + (ERI(p,q,c,d) + ERI(p,q,d,c))*X1(cd,ab)/ & 
                                   sqrt(1d0 + Kronecker_delta(c,d))
                 end do
              end do
          
              kl = 0
              do k=nC+1,nO
                 do l=k,nO
                    kl = kl + 1
                    rho1(p,q,ab) = rho1(p,q,ab) & 
                                 + (ERI(p,q,k,l) + ERI(p,q,l,k))*Y1(kl,ab)/ & 
                                   sqrt(1d0 + Kronecker_delta(k,l))
                 end do
              end do
              
           end do
           end do

           ij = 0
           do i=nC+1,nO
           do j=i,nO
              ij = ij + 1
              cd = 0
              do c=nO+1,nBas-nR
                 do d=c,nBas-nR
                    cd = cd + 1
                    rho2(p,q,ij) = rho2(p,q,ij) &
                                 + (ERI(p,q,c,d) + ERI(p,q,d,c))*X2(cd,ij)/ & 
                                   sqrt(1d0 + Kronecker_delta(c,d))
                 end do
              end do
              
              kl = 0
              do k=nC+1,nO
                 do l=k,nO
                    kl = kl + 1
                    rho2(p,q,ij) = rho2(p,q,ij) &
                                 + (ERI(p,q,k,l) + ERI(p,q,l,k))*Y2(kl,ij)/ & 
                                   sqrt(1d0 + Kronecker_delta(k,l))
                 end do
              end do
 
           end do
           end do

        end do
     end do
    !$OMP END DO
    !$OMP END PARALLEL
  end if

!----------------------------------------------
! Triplet manifold
!----------------------------------------------

  if(ispin == 2 .or. ispin == 4) then

    dim_1 = (nBas - nO) * (nBas - nO - 1) / 2
    dim_2 = nO * (nO - 1) / 2

    if((dim_1 .eq. 0) .or. (dim_2 .eq. 0)) then

      !$OMP PARALLEL DEFAULT(NONE)                                         &
      !$OMP          PRIVATE(p, q, a, b, ab, c, d, cd, i, j, ij, k, l, kl) &
      !$OMP          SHARED(nC, nBas, nR, nO, rho1, rho2, ERI, X1, Y1, X2, Y2)
      !$OMP DO COLLAPSE(2)
      do q = nC+1, nBas-nR
        do p = nC+1, nBas-nR
   
          ab = 0
  
          do a = nO+1, nBas-nR
            do b = a+1, nBas-nR
  
            ab = ab + 1
   
            cd = 0
            do c = nO+1, nBas-nR
              do d = c+1, nBas-nR
  
                cd = cd + 1
  
                rho1(p,q,ab) = rho1(p,q,ab) & 
                             + (ERI(p,q,c,d) - ERI(p,q,d,c))*X1(cd,ab) 
              end do ! d
            end do ! c
   
            kl = 0
            do k = nC+1, nO
              do l = k+1, nO
  
                kl = kl + 1
  
                rho1(p,q,ab) = rho1(p,q,ab) & 
                             + (ERI(p,q,k,l) - ERI(p,q,l,k))*Y1(kl,ab) 
              end do ! l
            end do ! k 
          end do ! b
        end do ! a
   
        ij = 0
        do i = nC+1, nO
          do j = i+1, nO
  
            ij = ij + 1
            
            cd = 0
            
            do c = nO+1, nBas-nR
              do d = c+1, nBas-nR
            
                cd = cd + 1
            
                rho2(p,q,ij) = rho2(p,q,ij) & 
                             + (ERI(p,q,c,d) - ERI(p,q,d,c))*X2(cd,ij) 
              end do ! d
            end do ! c
            
            kl = 0
            do k = nC+1, nO
              do l = k+1, nO
            
                kl = kl + 1
            
                rho2(p,q,ij) = rho2(p,q,ij) & 
                             + (ERI(p,q,k,l) - ERI(p,q,l,k))*Y2(kl,ij) 
              end do ! l
            end do ! k
          end do ! j
        end do ! i
      end do ! p
    end do ! q
    !$OMP END DO
    !$OMP END PARALLEL

    else

      allocate(ERI_1(nBas,nBas,dim_1), ERI_2(nBas,nBas,dim_2))
      ERI_1 = 0.d0
      ERI_2 = 0.d0
  
      !$OMP PARALLEL DEFAULT(NONE)                     &
      !$OMP          PRIVATE(p, q, c, d, cd, k, l, kl) &
      !$OMP          SHARED(nC, nBas, nR, nO, ERI_1, ERI_2, ERI)
      !$OMP DO COLLAPSE(2)
      do q = nC+1, nBas-nR
        do p = nC+1, nBas-nR
          cd = 0
          do c = nO+1, nBas-nR
            do d = c+1, nBas-nR
              cd = cd + 1
              ERI_1(p,q,cd) = ERI(p,q,c,d) - ERI(p,q,d,c)
            enddo
          enddo
          kl = 0
          do k = nC+1, nO
            do l = k+1, nO
              kl = kl + 1
              ERI_2(p,q,kl) = ERI(p,q,k,l) - ERI(p,q,l,k)
            end do
          end do
        enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL

      call dgemm("N", "N", nBas*nBas, dim_1, dim_1, 1.d0, &
                 ERI_1(1,1,1), nBas*nBas, X1(1,1), dim_1, &
                 0.d0, rho1(1,1,1), nBas*nBas)
  
      call dgemm("N", "N", nBas*nBas, dim_1, dim_2, 1.d0, &
                 ERI_2(1,1,1), nBas*nBas, Y1(1,1), dim_2, &
                 1.d0, rho1(1,1,1), nBas*nBas)
  
      call dgemm("N", "N", nBas*nBas, dim_2, dim_1, 1.d0, &
                 ERI_1(1,1,1), nBas*nBas, X2(1,1), dim_1, &
                 0.d0, rho2(1,1,1), nBas*nBas)
  
      call dgemm("N", "N", nBas*nBas, dim_2, dim_2, 1.d0, &
                 ERI_2(1,1,1), nBas*nBas, Y2(1,1), dim_2, &
                 1.d0, rho2(1,1,1), nBas*nBas)
  
      deallocate(ERI_1, ERI_2)
  
      rho1 = rho1
      rho2 = rho2

    endif
  endif

!----------------------------------------------
! alpha-beta block
!----------------------------------------------

  if(ispin == 3) then

    dim_1 = (nBas - nO) * (nBas - nO)
    dim_2 = nO * nO

    if((dim_1 .eq. 0) .or. (dim_2 .eq. 0)) then

      !$OMP PARALLEL DEFAULT(NONE)                                         &
      !$OMP          PRIVATE(p, q, a, b, ab, c, d, cd, i, j, ij, k, l, kl) &
      !$OMP          SHARED(nC, nBas, nR, nO, rho1, rho2, ERI, X1, Y1, X2, Y2)
      !$OMP DO COLLAPSE(2)
       
      do q = nC+1, nBas-nR
        do p = nC+1, nBas-nR
   
          ab = 0
          do a = nO+1, nBas-nR
            do b = nO+1, nBas-nR
  
              ab = ab + 1
   
              cd = 0
              do c = nO+1, nBas-nR
                do d = nO+1, nBas-nR
  
                  cd = cd + 1
  
                  rho1(p,q,ab) = rho1(p,q,ab) + ERI(p,q,c,d)*X1(cd,ab) 
                end do
              end do
   
              kl = 0
              do k = nC+1, nO
                do l = nC+1, nO
  
                  kl = kl + 1
  
                  rho1(p,q,ab) = rho1(p,q,ab) + ERI(p,q,k,l)*Y1(kl,ab) 
                end do
              end do
            end do
          end do
   
          ij = 0
          do i = nC+1, nO
            do j = nC+1, nO
  
              ij = ij + 1
               
              cd = 0
              do c = nO+1, nBas-nR
                do d = nO+1, nBas-nR
  
                  cd = cd + 1
  
                  rho2(p,q,ij) = rho2(p,q,ij) + ERI(p,q,c,d)*X2(cd,ij) 
                end do
              end do
   
              kl = 0
              do k = nC+1, nO
                do l = nC+1, nO
  
                  kl = kl + 1
  
                  rho2(p,q,ij) = rho2(p,q,ij) + ERI(p,q,k,l)*Y2(kl,ij) 
                end do
              end do
            end do
          end do
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

    else

      allocate(ERI_1(nBas,nBas,dim_1), ERI_2(nBas,nBas,dim_2))
      ERI_1 = 0.d0
      ERI_2 = 0.d0
  
      !$OMP PARALLEL DEFAULT(NONE)                     &
      !$OMP          PRIVATE(p, q, c, d, cd, k, l, kl) &
      !$OMP          SHARED(nC, nBas, nR, nO, ERI_1, ERI_2, ERI)
      !$OMP DO COLLAPSE(2)
      do q = nC+1, nBas-nR
        do p = nC+1, nBas-nR
          cd = 0
          do c = nO+1, nBas-nR
            do d = nO+1, nBas-nR
              cd = cd + 1
              ERI_1(p,q,cd) = ERI(p,q,c,d)
            enddo
          enddo
          kl = 0
          do k = nC+1, nO
            do l = nC+1, nO
              kl = kl + 1
              ERI_2(p,q,kl) = ERI(p,q,k,l)
            end do
          end do
        enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
  
      call dgemm("N", "N", nBas*nBas, dim_1, dim_1, 1.d0, &
                 ERI_1(1,1,1), nBas*nBas, X1(1,1), dim_1, &
                 0.d0, rho1(1,1,1), nBas*nBas)
  
      call dgemm("N", "N", nBas*nBas, dim_1, dim_2, 1.d0, &
                 ERI_2(1,1,1), nBas*nBas, Y1(1,1), dim_2, &
                 1.d0, rho1(1,1,1), nBas*nBas)
  
      call dgemm("N", "N", nBas*nBas, dim_2, dim_1, 1.d0, &
                 ERI_1(1,1,1), nBas*nBas, X2(1,1), dim_1, &
                 0.d0, rho2(1,1,1), nBas*nBas)
  
      call dgemm("N", "N", nBas*nBas, dim_2, dim_2, 1.d0, &
                 ERI_2(1,1,1), nBas*nBas, Y2(1,1), dim_2, &
                 1.d0, rho2(1,1,1), nBas*nBas)
  
      deallocate(ERI_1, ERI_2)

    endif
  endif

end subroutine 
