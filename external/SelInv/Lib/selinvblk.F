C
C  SELINVBLK computes the selected elements of inv(A) which has the same
C  sparsity pattern of the Cholesky factor of A.
C
C
C   INPUT PARAMETERS:
C       NSUPER          -   NUMBER OF SUPERNODES.
C       XSUPER          -   SUPERNODE PARTITION.
C       (XLINDX,LINDX)  -   ROW INDICES FOR EACH SUPERNODE.
C       (XLNZ,LNZ)      -   CHOLESKY FACTOR.
C   OUTPUT PARAMETERS:
C       COLPTR          -   Column partition of the selected components
C                           of invA
C       ROWIND          -   Row indicies of the selected components of 
C                           invA
C       LNZ             -   Values of the selected components of invA
C  
C  NOTE: 
C    CY: This is a pointwise implementation 
C    LL: The difference between selinvblk.f and exdiagblk.f is that
C        selinvblk.f returns both the diagonal as well as the off-diagonal
C        elements evaluated in the selected inversion process.
C  
C  Authors:
C  Chao Yang and Lin Lin
C  Computer Research Division, Lawrence Berkeley National Lab
C  Last modified:  09-12-2011
C
C***********************************************************************
C
      SUBROUTINE  SELINVBLK ( NSUPER, XSUPER, XLINDX, LINDX, XLNZ  ,
     &                        LNZ   , SNODES, PERM , NEQNS ,
     &                        colptr, rowind, dumpL )
      implicit none
C
C***********************************************************************
C
        INTEGER             NSUPER, NEQNS
        INTEGER             LINDX(*)      , XSUPER(*), PERM(*)
        INTEGER             XLINDX(*)     , XLNZ(*), SNODES(*)
        COMPLEX*16          LNZ(*)        
        INTEGER             dumpL
C LL
        INTEGER             rowind(*)     , colptr(*)
C***********************************************************************
C
        INTEGER             FJCOL , I     , IPNT  , IX    , IXSTOP,
     &                      IXSTRT, JCOL  , JPNT  , JSUP  , LJCOL
        COMPLEX*16    T, dinv
chao
        integer irow, snsize, idxrel, jxstrt, jxstop, jx, isup, 
     &          idest
        integer imat, ivec, imatch, vecrow, matrow, idel, nzvec, lx
        REAL                GTIMER, T0, T1
chaoadd
        complex*16,  allocatable, dimension(:)   :: dwork, ywork
        integer, allocatable, dimension(:)   :: newxlnz, ipiv
        integer, allocatable, dimension(:,:) :: blkmap
        integer  supsize, nnzlplus, colnnz, nrows, ncols, ix0, isup0
        integer  ierr, maxsup, ldsupi, ldsupj, nvrows, maxwork,
     &           maxldsup, matsize, nblks, irow0, ivblk, ib, ldsupk, kx,
     &           kxstrt, kxstop, iy, iy0, ldy
        complex*16   zone, zzero
        parameter (zone = (1.0d0,0.0d0), zzero = (0.0d0,0.0d0) )
        complex*16 zdotu
        real     tcopy0, tcopy1
        integer  l1, l2, ptr
C LL    
C       filnz: first nonzero index of L. 
C       lilnz: last nonzero index of L.
C       colord: reordered column indices
C       roword: reordered row indices
C       (iindx, jindx, inva): unsorted (i,j,val) triplet for invA.

        integer filnz, lilnz, cnt, pirow, pjcol, nnzl
        integer, allocatable, dimension(:):: colord, roword, rowtmp
        integer, allocatable, dimension(:):: iindx, jindx
        COMPLEX*16, allocatable, dimension(:):: inva, invatmp


 
C
C***********************************************************************
C
        t0 = gtimer()
        IF  ( NSUPER .LE. 0 )  RETURN
        if (dumpL .eq. 1) then
           call dumpLmat(nsuper, xsuper, xlindx, lindx, xlnz, lnz,
     &                  'debugL.m')  
        endif
        ! 
        ! find out how much extra storage needed
        !
        nnzlplus = xlnz(neqns+1)-1
#ifdef __DEBUG
        write(6,*) 'number of nonzeros = ', nnzlplus
#endif

        maxsup = 0
        do jsup = 1, nsuper
           supsize = xsuper(jsup+1)-xsuper(jsup)
           if (supsize .gt. maxsup) maxsup = supsize
           nnzlplus = nnzlplus + supsize*(supsize-1)/2
        end do
#ifdef __DEBUG
        write(6,*) 'increase to = ', nnzlplus
        write(6,*) 'maxsup =', maxsup
#endif
        !
        allocate(newxlnz(neqns+1))
        allocate(ipiv(maxsup))
        allocate(dwork(2*maxsup))
        do i = 1, maxsup
           ipiv(i) = i 
        end do
        !
        ! copy L and setup the pointer
        !  
        tcopy0 = gtimer()
        newxlnz(neqns+1) = nnzlplus+1
        maxwork = 0
        maxldsup = 0
        do jsup = nsuper,1,-1
           fjcol = xsuper(jsup)
           ljcol = xsuper(jsup+1)-1
           colnnz = xlnz(fjcol+1) - xlnz(fjcol)
           do jcol = ljcol,fjcol,-1
              ixstrt = xlnz(jcol)
              ixstop = xlnz(jcol+1)-1
              do ix = ixstop, ixstrt, -1
                  !
                  ! destination must be below the diagonal
                  !
                  idest = newxlnz(jcol+1)-1
                  lnz(idest+(ix-ixstop)) = lnz(ix)
              end do
              newxlnz(jcol) = newxlnz(jcol+1) - colnnz
           end do
           if (jsup .lt. nsuper) then
              matsize = colnnz*(ljcol-fjcol+1)
              if (colnnz .gt. maxldsup) maxldsup = colnnz
              if (matsize .gt. maxwork) maxwork = matsize
           endif
        end do
#ifdef __DEBUG
        write(6,*) 'newnnzl = ', newxlnz(neqns+1)-1
        write(6,*) 'maxwork = ', maxwork
#endif
        

        !
        !  modify all supernodes
        !
        do jsup = 1, nsuper
           ierr = 0
           fjcol = xsuper(jsup)
           ljcol = xsuper(jsup+1)-1
           supsize = ljcol-fjcol+1
           colnnz = newxlnz(fjcol+1) - newxlnz(fjcol)
           nrows = colnnz-supsize
           ncols = supsize
           ixstrt = newxlnz(fjcol) 
           if (supsize .gt. 1) then 
              ! 
              !  perform a triangular solve below the diagonal block
              ! 
              call ztrsm('Right','Lower','Notranspose','Unit', 
     &                   nrows, ncols, zone, lnz(ixstrt), colnnz, 
     &                   lnz(ixstrt+supsize), colnnz)
              !
              !  invert the diagonal block
              !  

C              ptr = ixstrt
C              do l2 = 1, supsize
C                do l1 = 1, supsize
C                  if( l1 == l2 ) then
C                    write(6, *) '(', l1, ', ', l2, ') = ', lnz(ptr)
C                  endif
C                  ptr = ptr + 1
C                enddo
C                ptr = ptr + colnnz - supsize
C              enddo

              call zsytri('Lower',supsize,lnz(ixstrt),colnnz,
     &                    ipiv, dwork, ierr)

                  
              if (ierr .ne. 0) then
                 write(6,*) 'triangular solve failed, ierr = ', ierr
                 write(6,*) 'jsup = ', jsup, 'supsize =', supsize
              endif 
              !
              ! need to store the upper triangular part of the diagonal
              ! block for for dgemm
              !
              do irow = 1,supsize-1
                 do jcol = irow+1,supsize
                    lnz(ixstrt+(jcol-1)*colnnz+irow-1) 
     &              = lnz(ixstrt+(irow-1)*colnnz + jcol-1)
                 end do
              end do
           else
C              write(6, *) '(', ixstrt, ', ', ixstrt, ') = ', lnz(ixstrt)
              lnz(ixstrt) = 1/lnz(ixstrt)
           endif
        end do
        tcopy1 = gtimer()
#ifdef __DEBUG
        write(6,123) tcopy1 - tcopy0
 123    format(1x,' Time copy = ', 1pe11.3)
#endif
C        do i = 1, supsize
C           write(6,*) 'diag = ', 
C     &                lnz(ixstrt+colnnz*(i-1)+i-1)
C        end do    

        allocate(ywork(maxwork))
        allocate(blkmap(maxldsup,3))

C       LL: Initialize the variable, EXTREMELY IMPORTANT! 
        do i = 1, maxwork
          ywork(i) = 0
        enddo
        do jsup = nsuper-1, 1, -1
           !write(6,*) 'jsup = ', jsup
           fjcol = xsuper(jsup)
           ljcol = xsuper(jsup+1)-1
           supsize = ljcol-fjcol+1
           !write(6,*) 'supernode size = ', supsize
           jpnt = xlindx(jsup)+supsize ! point to the beginning of
                                       ! the off-diagonal blocks  
                                       ! column
           ixstrt = newxlnz(fjcol)+supsize 
           ixstop = newxlnz(fjcol+1)-1
           ldsupj = newxlnz(fjcol+1)-newxlnz(fjcol)
           ldy    = ixstop - ixstrt + 1
           ! LL: do not do anything if there is no off diagonal elements
           if(ldy .gt. 0 ) then
           !
           ! construct a blockmap list for the jsup-th supernode
           ! a block row is defined to be a set of rows within 
           ! the same supernode that have consequtive row indices
           !
           ! blkmap(ib,1) --- row index of the first row in the block
           ! blkmap(ib,2) --- the supernode the block belongs to
           ! blkmap(ib,3) --- number of rows in the block
           !
           nblks = 0
           irow0 = -1
           isup0 = -1
           nrows = 0
           ipnt  = jpnt
           do ix = ixstrt, ixstop
              irow = lindx(ipnt)
              isup = snodes(irow)
              if (irow .eq. irow0 + 1 .and. isup .eq. isup0 ) then
                 !
                 ! within the same block
                 !
                 nrows = nrows + 1
              else
                 !
                 ! start a new block
                 !
                 isup0 = isup
                 nblks = nblks + 1
                 blkmap(nblks,1) = irow
                 blkmap(nblks,2) = isup
                 if (nblks .gt. 1) then
                    blkmap(nblks-1,3) = nrows
                 endif 
                 nrows = 1
              endif
              irow0 = irow
              ipnt = ipnt + 1               
           enddo
           if (nblks .ge. 1) then
              blkmap(nblks,3) = nrows 
           endif 
           !write(6,*) 'nblks = ', nblks
           !do i = 1, nblks
           !   write(6,*) blkmap(i,1),blkmap(i,2),blkmap(i,3)
           !end do 
           !
           !  do Y = -S*X, 
           !
           !  ix   pointer to the nonzero value array
           !  kx   pointer to the nonzero value array
           !  iy   pointer to the nonzero value work array y
           !  imat pointer to row index array
           !
           !  step thru blocks of X, each block belongs to a supernode
           !
           if (nblks .gt. 0) then
              irow0 = blkmap(1,1)
              iy0 = 1 
              ix0 = ixstrt
              do ib = 1, nblks
                 irow  = blkmap(ib,1)
                 isup  = blkmap(ib,2)
                 ncols = blkmap(ib,3) ! it is the numer of rows in a nonzero
                                      ! block of X. It corresponds to the number
                                      ! of columns in the trailing Schur 
                                      ! complement that will be used
                 ivblk = ib
                 !
                 kxstrt = newxlnz(irow) 
                 kxstop = newxlnz(irow+1)-1
                 imat   = xlindx(isup)  ! pointer to the row index array that 
                                        ! corresponds to the first nonzero
                                        ! entry of the isup-th supernode
                 ldsupk = kxstop-kxstrt+1
                 imatch = 0
                 iy = iy0
                 ix = ix0
                 kx = kxstrt
                 !
                 ! wall down the irow-th column (which belongs to isup-th
                 ! supernode) of the trailing Schur complement. Perform 
                 ! matrix operations when the row index of the Schur complement
                 ! matches the first row of a non-zero block of X
                 !
                 do while (kx .le. kxstop)
                    matrow = lindx(imat)
                    vecrow = blkmap(ivblk,1)
                    nvrows = blkmap(ivblk,3) 
                    if (matrow .lt. vecrow) then
                       ! 
                       ! skip this row
                       !
                       imat = imat + 1
                       kx   = kx + 1 
                    else 
                       ! matrow must match vecrow
                       !
                       ! the lower triangular contribution (axpy)
                       !
                       imatch = imatch + 1
                       !
                       ! separate the following two cases for performance
                       ! reason
                       !
                       if (supsize .eq. 1) then
                          if ( ncols .eq. 1) then
                             call zaxpy(nvrows,-lnz(ix), lnz(kx),
     &                                   1, ywork(iy), 1)
                          else
                             call zgemv('N',nvrows,ncols,-zone,
     &                                  lnz(kx),ldsupk, 
     &                                  lnz(ix), 1, zone,
     &                                  ywork(iy),1)
                          endif
                       else
                          call zgemm('N','N',nvrows,supsize,ncols,
     &                               -zone,lnz(kx),ldsupk,
     &                               lnz(ix),ldsupj,zone, 
     &                               ywork(iy), ldy)
                       endif
                       !
                       ! the upper triangular contribution (dot)
                       !
                       if (imatch .gt. 1) then
                          !
                          ! separate the following two cases
                          ! for performance reason
                          ! 
                          if (supsize .eq. 1) then
                             if (ncols .eq. 1) then
                                ywork(iy0) = ywork(iy0)
     &                                   - zdotu(nvrows,lnz(kx),1,
     &                                           lnz(ix+iy-iy0),1)
                             else
                                call zgemv('T',nvrows,ncols,-zone,
     &                                     lnz(kx),ldsupk,
     &                                     lnz(ix+iy-iy0),1,
     &                                     zone, ywork(iy0),1)
                             endif
                          else
                             call zgemm('T','N',ncols,supsize,nvrows,
     &                                  -zone,lnz(kx),ldsupk,
     &                                  lnz(ix+iy-iy0),ldsupj,
     &                                  zone, ywork(iy0),ldy)
                          endif
                       endif
                       !
                       ! move on to the next block in the Schur complement
                       !
                       ivblk = ivblk + 1
                       imat = imat + nvrows
                       kx   = kx + nvrows
                       iy   = iy + nvrows
                       if (imatch .ge. nblks-ib+1) goto 20
                    endif
                 enddo !while kx
 20              continue
                 !
                 ! remember that ncols is actually the number of rows
                 ! in a nonzero block of X that have just been processed.
                 !
                 iy0 = iy0 + ncols
                 ix0 = ix0 + ncols
              end do ! ib
           endif ! nblks > 0
           !
           ! do G = G + X'*SX
           !
           ix = ixstrt
           jx = newxlnz(fjcol)
           iy = 1
           do ib = 1, nblks
              nrows = blkmap(ib,3)
              if (supsize .eq. 1) then
                lnz(jx) = lnz(jx) 
     &                      - zdotu(nrows,lnz(ix),1,ywork(iy),1)  
              else
                 call zgemm('T','N',supsize,supsize,nrows,-zone,
     &                      lnz(ix), ldsupj, ywork(iy), ldy,
     &                      zone, lnz(jx), ldsupj)
              endif
              ix = ix + nrows
              iy = iy + nrows
           end do  
           !do jcol = fjcol,ljcol
           !   write(6,*) 'diag = ', lnz(newxlnz(jcol)+jcol-fjcol)
           !end do
           !
           ! copy y to lnz
           !
           call zlacpy('All',(ixstop-ixstrt+1),supsize,ywork,ldy,
     &                 lnz(ixstrt),ldsupj)
           ! 
           ! must zero out the workspace
           !
           do i = 1, ldy*supsize
              ywork(i) = zzero
           enddo
         endif ! if (ldy .gt. 0)
        enddo

        !
        ! LL: Reshuffle LNZ into the output (IINDX, JINDX, INVA)
        !
        nnzl = xlnz(neqns+1)-1
        allocate(iindx(nnzl))
        allocate(jindx(nnzl))
        allocate(inva(nnzl))
        allocate(colord(nnzl))
        allocate(roword(neqns))
        allocate(rowtmp(neqns))
        allocate(invatmp(neqns))

        cnt = 1
        do jsup = 1, nsuper
           fjcol = xsuper(jsup)
           ljcol = xsuper(jsup+1)-1
           do jcol = fjcol,ljcol
              filnz = newxlnz(jcol) + jcol-fjcol
              lilnz = newxlnz(jcol+1) - 1
              jpnt  = xlindx(jsup) + jcol - fjcol
              do jx = filnz, lilnz
                pirow = perm(lindx(jpnt))
                pjcol = perm(jcol)
                if( pirow > pjcol ) then
                  iindx(cnt) = pirow
                  jindx(cnt) = pjcol
                else
                  iindx(cnt) = pjcol
                  jindx(cnt) = pirow
                endif 
                inva(cnt)  = lnz(jx)
C                write(*,*) jx, lindx(jpnt), jcol, lnz(jx)
                cnt  = cnt + 1
                jpnt = jpnt + 1
              enddo
           enddo
        enddo


        if(nnzl  .ne. cnt-1 ) then
          write(6,*) 'The number of nonzeros do not match'
        endif 
        ! 
        ! LL: Write back (IINDX, JINDX, INVA) to (COLPTR, ROWIND, LNZ) in
        ! Column-major order
        ! NOTE: The size of LNZ is shrinked back to nnzl compared to
        ! nnzlplus in the computation phase.
        ! 
        ! Sort column first.
        call sortix(nnzl, jindx, colord)
C        call qsorti(colord,nnzl, jindx)

        colptr(1) = 1
        do jcol = 1, neqns
          cnt = 1
          jx  = colord(colptr(jcol))
          do while(jindx(jx) == jcol)
            rowtmp(cnt)  = iindx(jx)
            invatmp(cnt) = inva(jx)
            cnt = cnt + 1
            lx  = colptr(jcol)+cnt-1
            if( lx > nnzl ) then
              exit
            endif
            jx  = colord(lx)
          enddo 
          colptr(jcol+1) = colptr(jcol) + cnt - 1
          ! Sort row next
          call sortix(cnt-1, rowtmp, roword)
C        call qsorti(roword, cnt-1, rowtmp)

        do irow = 1, cnt-1
            jx = colptr(jcol)+irow-1
            rowind(jx)    =  rowtmp(roword(irow))
            ! Save data back to lnz
            lnz(jx)       =  invatmp(roword(irow))
              
          end do
        enddo


        t1 = gtimer()
#ifdef __DEBUG
        write(6,333) t1-t0
  333   format(1x,'Extraction time (including copy) = ',1pe11.3)
#endif

        deallocate(newxlnz)
        deallocate(ipiv)
        deallocate(dwork)
        deallocate(ywork)
        deallocate(blkmap)
        deallocate(iindx)
        deallocate(jindx)
        deallocate(inva)
        deallocate(roword)
        deallocate(colord)
        deallocate(rowtmp)
        deallocate(invatmp)
        RETURN
      END
