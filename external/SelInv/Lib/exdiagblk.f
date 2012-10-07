C
C  EXTRACT  THE DIAGONAL OF THE INVERSE OF A
C  This is a pointwise implementation (CY)
C
C   INPUT PARAMETERS:
C       NSUPER          -   NUMBER OF SUPERNODES.
C       XSUPER          -   SUPERNODE PARTITION.
C       (XLINDX,LINDX)  -   ROW INDICES FOR EACH SUPERNODE.
C       (XLNZ,LNZ)      -   CHOLESKY FACTOR.
C
C   UPDATED PARAMETERS:
C       SCHURC             -   ON INPUT, CONTAINS THE RIGHT HAND SIDE.  ON
C                           OUTPUT, CONTAINS THE SOLUTION.
C
C***********************************************************************
C
      SUBROUTINE  EXDIAGBLK ( NSUPER, XSUPER, XLINDX, LINDX, XLNZ  ,
     &                        LNZ   , SNODES, DIAG  , PERM , NEQNS ,
     &                        dumpL)
      implicit none
C
C***********************************************************************
C
        INTEGER             NSUPER, NEQNS
        INTEGER             LINDX(*)      , XSUPER(*), PERM(*)
        INTEGER             XLINDX(*)     , XLNZ(*), SNODES(*)
        COMPLEX*16          LNZ(*)        , DIAG(*)
        INTEGER             dumpL
C
C***********************************************************************
C
        INTEGER             FJCOL , I     , IPNT  , IX    , IXSTOP,
     &                      IXSTRT, JCOL  , JPNT  , JSUP  , LJCOL
        COMPLEX*16    T, dinv
chao
        integer irow, snsize, idxrel, jxstrt, jxstop, jx, isup, idest
        integer imat, ivec, imatch, vecrow, matrow, idel, nzvec
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
        write(6,*) 'number of nonzeros = ', nnzlplus
        maxsup = 0
        do jsup = 1, nsuper
           supsize = xsuper(jsup+1)-xsuper(jsup)
           if (supsize .gt. maxsup) maxsup = supsize
           nnzlplus = nnzlplus + supsize*(supsize-1)/2
        end do
        !write(6,*) 'increase to = ', nnzlplus
        !
        allocate(newxlnz(neqns+1))
        allocate(ipiv(maxsup))
        allocate(dwork(maxsup))
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
        write(6,*) 'newnnzl = ', newxlnz(neqns+1)-1
        write(6,*) 'maxwork = ', maxwork
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
              lnz(ixstrt) = 1/lnz(ixstrt)
           endif
        end do
        tcopy1 = gtimer()
        write(6,123) tcopy1 - tcopy0
 123    format(1x,' Time copy = ', 1pe11.3)
        !do i = 1, supsize
        !   write(6,*) 'diag = ', lnz(ixstrt+colnnz*(i-1)+i-1)
        !end do    

        allocate(ywork(maxwork))
        allocate(blkmap(maxldsup,3))
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
        enddo
        !
        ! permute the diagonal back to its original order
        !
        do jsup = 1, nsuper
           fjcol = xsuper(jsup)
           ljcol = xsuper(jsup+1)-1
           do jcol = fjcol,ljcol
              jx = newxlnz(jcol)+jcol-fjcol
              diag(perm(jcol)) = lnz(jx)
           enddo
        enddo
        t1 = gtimer()
        write(6,333) t1-t0
  333   format(1x,'Extraction time (including copy) = ',1pe11.3)

        deallocate(newxlnz)
        deallocate(ipiv)
        deallocate(dwork)
        deallocate(ywork)
        deallocate(blkmap)
        RETURN
      END
