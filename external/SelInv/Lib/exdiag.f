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
      SUBROUTINE  EXDIAG (  NSUPER, XSUPER, XLINDX, LINDX, XLNZ  ,
     &                      LNZ   , SNODES, DIAG  , Y    , PERM  ,
     &                      NEQNS )
      implicit none
C
C***********************************************************************
C
        INTEGER             NSUPER, NEQNS
        INTEGER             LINDX(*)      , XSUPER(*), PERM(*)
        INTEGER             XLINDX(*)     , XLNZ(*), SNODES(*)
        DOUBLE PRECISION    LNZ(*)        , DIAG(*), Y(*)
C
C***********************************************************************
C
        INTEGER             FJCOL , I     , IPNT  , IX    , IXSTOP,
     &                      IXSTRT, JCOL  , JPNT  , JSUP  , LJCOL
        DOUBLE PRECISION    T, dinv
chao
        integer irow, snsize, idxrel, jxstrt, jxstop, jx, isup
        integer imat, ivec, imatch, vecrow, matrow, idel, nzvec
        REAL                GTIMER, T0, T1
C
C***********************************************************************
C
        t0 = gtimer()
        IF  ( NSUPER .LE. 0 )  RETURN

!        call dumpLmat(nsuper, xsuper, xlindx, lindx, xlnz, lnz,
!     &                'debugL.m')  

!
!       The last supernode is dense
!
        jsup = nsuper
        !write(6,*) 'jsup = ', jsup
        FJCOL = XSUPER(JSUP)
        LJCOL = XSUPER(jsup+1)-1
        !write(6,*) 'supernode size = ', LJCOL-FJCOL+1
        jpnt = xlindx(jsup)+ljcol-fjcol-1 ! pointer to row index array
                                          ! that corresponds to the second 
                                          ! to last column
        diag(ljcol) = 1/lnz(xlnz(ljcol)) 
        lnz(xlnz(ljcol)) = diag(ljcol)
        do jcol = ljcol-1,fjcol,-1
           ixstrt = xlnz(jcol) 
           ixstop = xlnz(jcol+1)-1
           ipnt = jpnt + 1
           !
           ! perform a matvec
           ! going through the vector             
           do ix = ixstrt+1, ixstop
              irow = lindx(ipnt)
              !
              ! pick out the irow-th column
              !
              jxstrt = xlnz(irow) 
              jxstop = xlnz(irow+1)-1
              !
              do jx = jxstrt,jxstop
                 idxrel = jx-jxstrt+(ix-ixstrt)  ! relative index
                 !
                 ! the lower triangular contribution (axpy)
                 !
                 y(idxrel) = y(idxrel)+lnz(jx)*lnz(ix)
                 !
                 ! the upper triangular contribution (dot)
                 !
                 if (jx .gt. jxstrt) then
                    y(ix-ixstrt) = y(ix-ixstrt) 
     &                           + lnz(jx)*lnz(ix+jx-jxstrt) 
                 endif
              end do
              !
              ipnt = ipnt + 1
           enddo
           jpnt = jpnt - 1
           !
           ! dot product for updating the inverse
           !
           diag(jcol) = 1/lnz(xlnz(jcol)) 
           do ix = ixstrt+1,ixstop
              idxrel = ix - ixstrt
              diag(jcol) = diag(jcol) + y(idxrel)*lnz(ix)
           end do  
           lnz(xlnz(jcol)) = diag(jcol)
           !
           ! copy y to lnz
           !
           do ix = ixstrt+1,ixstop
              lnz(ix) = -y(ix-ixstrt)
              y(ix-ixstrt) = 0.0
           enddo
        end do
        !do jcol = fjcol,ljcol
        !   write(6,111) jcol, diag(jcol)  
        !end do 
 111    format(1x,'diag(',I10,') = ',1pe11.4)
        write(6,*)
!
        do jsup = nsuper-1, 1, -1
           !write(6,*) 'jsup = ', jsup
           fjcol = xsuper(jsup)
           ljcol = xsuper(jsup+1)-1
           !write(6,*) 'supernode size = ', ljcol-fjcol+1
           jpnt = xlindx(jsup)+ljcol-fjcol ! pointer to row index array
                                           ! that corresponds to the last 
                                           ! column
           !
           ! update diag(jcol) using computed Schur complement 
           ! stored in the trailing LNZ
           !
           do jcol = ljcol,fjcol,-1
              !write(6,*) 'jcol = ', jcol 
              ixstrt = xlnz(jcol) 
              ixstop = xlnz(jcol+1)-1
              ipnt = jpnt + 1  ! start from the first non-diagonal element
              !
              ! perform a matvec
              ! step thru the jcol-th column of L, do axpy and dot simultaneously 
              !
              nzvec = ixstop - ixstrt  ! number of nonzeros in the jcol-th column
                                       ! (below the diagonal). This is also the
                                       ! dimension of the matrix to be multiplied
              !write(6,*) 'nzvec = ', nzvec 
              do ix = ixstrt+1, ixstop
                 irow = lindx(ipnt)
                 !write(6,*) '    irow = ', irow 
                 !
                 ! pick out the irow-th column of the trailing Schur complement
                 !
                 jxstrt = xlnz(irow) 
                 jxstop = xlnz(irow+1)-1
                 !
                 ! find out which supernode it belongs to
                 ! 
                 isup = snodes(irow)
                 idel = irow - xsuper(isup)! offset from the first column
                                           ! of this supernode 
                 imat = xlindx(isup)+idel  ! pointer to the first useful index
                                           ! in the row index array that provides
                                           ! row index associated with the first
                                           ! entry of the matrix to be multiplied 
                 !write(6,*) 'isup = ', isup, 'idel = ', idel
                 imatch = 0
                 ivec = ipnt
                 !
                 ! Now step thru the irow-th column of the trailing
                 ! Schur complement, and extract only those elements
                 ! that need to be multiplied.
                 !
                 ! The nonzero pattern of the jcol-th column of L
                 ! should be a subset of that associated with each 
                 ! column of the trailing Schur complement
                 !
                 do jx = jxstrt,jxstop
                    idxrel = jx-jxstrt+(ix-ixstrt)  ! relative index
                    matrow = lindx(imat)
 10                 continue
                    vecrow = lindx(ivec) 
                    if (matrow < vecrow) then
                       ! 
                       ! skip this row
                       !
                       imat = imat + 1
                       cycle
                    else 
                       ! matrow must match vecrow
                       !
                       ! the lower triangular contribution (axpy)
                       !
                       imatch = imatch + 1
                       y(imatch+ix-ixstrt-1) 
     &                     = y(imatch+ix-ixstrt-1)+lnz(jx)*lnz(ix)
                       !
                       ! the upper triangular contribution (dot)
                       !
                       if (imatch .gt. 1) then
                          y(ix-ixstrt) = y(ix-ixstrt) 
     &                           + lnz(jx)*lnz(ix+imatch-1) 
                       endif
                       ivec = ivec + 1
                       imat = imat + 1 
                       if (imatch .ge. nzvec) goto 20
                    endif
                 enddo
 20              continue
                 !
                 ipnt = ipnt + 1
                 nzvec = nzvec - 1
              enddo
              jpnt = jpnt - 1
              !
              ! dot product for updating the inverse
              !
              dinv = 1/lnz(xlnz(jcol)) 
              do ix = ixstrt+1,ixstop
                 idxrel = ix - ixstrt
                 dinv = dinv + y(idxrel)*lnz(ix)
              end do  
              lnz(xlnz(jcol)) = dinv
              !
              ! copy y to lnz
              !
              do ix = ixstrt+1,ixstop
                 lnz(ix) = -y(ix-ixstrt)
                 y(ix-ixstrt) = 0.0
              enddo
           enddo
        enddo
        !
        ! permute the diagonal back to its original order
        !
        do jcol = 1, neqns
           diag(perm(jcol)) = lnz(xlnz(jcol))
        enddo
        t1 = gtimer()
        write(6,333) t1-t0
  333   format(1x,'Extraction time = ',1pe11.3)
        RETURN
      END
