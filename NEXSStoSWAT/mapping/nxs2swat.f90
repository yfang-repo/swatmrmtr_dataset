program tree_sort
! Sorts a file of integers by building a
! tree, sorted in infix order.
! This sort has expected behavior n log n,
! but worst case (input is sorted) n ** 2.
! combine swat and nexss which has the same decreasing or increasing trend and
! grid number with +/- 4 grid cells.  -- not so good
! consider if there is a node landing in the swat reach, add the rest too.
! sum up nexss nhd segment flux

   implicit none
   type nxs_node
      type (nxs_node_fld), pointer :: t
      real*8 :: rate(2)
      real*8 :: area
      real*8 :: xsarea !cross sectional area
      real*8 :: width
      real*8 :: tlength
      integer :: fnd = 0
   end type nxs_node

   type nxs_node_fld
      type (nxs_node_fld), pointer :: left, right
      integer :: FID_reach
      real*8 :: qhef_lat_m
      real*8 :: qhef_ver_m
      real*8 :: qhef_total
      real*8 :: w_ma 
      real*8 :: d_ma 
      real*8 :: Lreach 
      real*8 :: t50_lat_ma
      real*8 :: t50_ver_ma
      real*8 :: leng_int
      real*8 :: frac
      logical :: fnd = .false.
   end type nxs_node_fld
   
   type (nxs_node), dimension(:), pointer :: rch  ! A tree

   type swat_node
      type (swat_node_fld), pointer :: t
      integer :: nxs_rch
      real*8 :: tlength
      integer :: skip = 1
   end type swat_node
   
   type swat_node_fld
      type (swat_node_fld), pointer :: left, right
      integer :: FID_reach
      real*8 :: leng_int
   end type swat_node_fld

   type (swat_node), dimension(:), pointer :: swat_rch  ! A tree

   integer :: number, ios
   integer :: fid_nhdplu
   real*8 :: w_ma
   real*8 :: d_ma
   real*8 :: wbkf
   real*8 :: t50_ver_ma
   real*8 :: t50_lat_ma
   real*8 :: qhef_ver_m
   real*8 :: qhef_lat_m
   real*8 :: Lreach
   integer :: FID_reach ! structured grid
   real*8 :: leng_int ! length of line segment intersect with grid
   integer :: i
   integer :: nxs_rch
   integer :: swat_rchs
   real*8 :: length_tmp ! total length from geometry segment, to calculate proportion of rate on each seg

! swat variables
   integer :: FID_riv1
   real*8 :: length_int
   real*8 :: qsum(2)
   real*8 :: tqsumv
   real*8 :: tqsuml
   integer :: ich
   integer :: overlap ! on the same if overlap by two nodes
   integer :: j !nexss rch counter
   real*8 :: xsarea

   open(10,file='nexss_gridn.txt',status='old')
   read(10,*) !skip names
! from the first number in the txt file + 1
   nxs_rch = 677
   swat_rchs = 101

   allocate(rch(nxs_rch))
!FID_riv1	FID_reach_	length_int
   
   ! Start with empty tree
   ! for nexss
   do i=1,nxs_rch
     nullify (rch(i)%t)  
   enddo
   do
      read (10, *, iostat = ios) fid_nhdplu,w_ma,wbkf,d_ma,t50_ver_ma,t50_lat_ma, &
          qhef_ver_m, qhef_lat_m, Lreach, FID_reach, leng_int
      if (ios < 0) exit
      rch(fid_nhdplu+1)%rate(1) = qhef_ver_m
      rch(fid_nhdplu+1)%rate(2) = qhef_lat_m
!area
      rch(fid_nhdplu+1)%area = Lreach*w_ma
      rch(fid_nhdplu+1)%width = w_ma
      rch(fid_nhdplu+1)%xsarea = w_ma*d_ma
!      rch(fid_nhdplu+1)%depth = d_ma
      rch(fid_nhdplu+1)%fnd = 0
!residence time

      call insert (rch(fid_nhdplu+1)%t, w_ma,t50_ver_ma,t50_lat_ma,qhef_ver_m, &
           qhef_lat_m, Lreach, FID_reach, leng_int) ! Put next number in tree
   end do
   ! Print nodes of tree in infix order
   do i=1,nxs_rch
     length_tmp = 0.d0
     call print_length_tree (rch(i)%t)
     rch(i)%tlength  = length_tmp
!     call print_frac_tree (rch(i)%t, rch(i)%tlength)
   enddo

   open(11,file='reach_grid.txt',status='old')
   read(11,*) !skip names
   allocate(swat_rch(swat_rchs))
   
   ! Start with empty tree
   do i=1,swat_rchs
     nullify (swat_rch(i)%t)  
   enddo
   do
!FID_riv1	FID_reach_	length_int
      read (11, *, iostat = ios) FID_riv1, FID_reach, length_int
      if (ios < 0) exit
      call swat_insert (swat_rch(FID_riv1+1)%t, FID_reach, length_int) ! Put next number in tree
   end do
!
   do i=1,swat_rchs
     length_tmp = 0.d0
     call print_length_swat_tree (swat_rch(i)%t)
     swat_rch(i)%tlength  = length_tmp
   enddo
!
   do i=1,swat_rchs
     qsum(:) = 0.d0
     tqsumv = 0.d0
     tqsuml = 0.d0
     xsarea = 0.d0
       ich = 0
     do j=1,nxs_rch
       overlap = 0
       call nxs_to_swat (swat_rch(i)%t,j,swat_rch(i)%tlength,qsum,tqsumv,tqsuml,ich,xsarea)
!if(overlap>2) print *,'overlap1-',overlap,i,j
     enddo
     if(qsum(1) == 0.d0) then
        tqsumv = 0.d0
     else
        tqsumv = tqsumv/qsum(1)
     endif
     if(qsum(2) == 0.d0) then
        tqsuml = 0.d0
     else
       tqsuml = tqsuml/qsum(2)
     endif
     if(qsum(1) + qsum(2) == 0.d0) then
       xsarea = 100.d0
     else
       xsarea = xsarea/(qsum(1)+qsum(2))
     endif
     write(*,100) qsum(1:2),tqsumv,tqsuml,xsarea,ich
   enddo
   100 format(5e12.3,i5)
!
contains

   recursive subroutine insert (t,w_ma,t50_ver_ma,t50_lat_ma, &
      qhef_ver_m, qhef_lat_m, Lreach, FID_reach, leng_int)

      type (nxs_node_fld), pointer :: t  ! A tree
      real*8, intent(in) :: w_ma
      real*8, intent(in) :: t50_ver_ma
      real*8, intent(in) :: t50_lat_ma
      real*8, intent(in) :: qhef_ver_m
      real*8, intent(in) :: qhef_lat_m
      real*8, intent(in) :: Lreach
      real*8, intent(in) :: leng_int
      integer, intent(in) :: FID_reach
      ! If (sub)tree is empty, put number at root
      if (.not. associated (t)) then
         allocate (t)
         t % FID_reach = FID_reach
           t % w_ma = w_ma
           t % t50_ver_ma = t50_ver_ma
           t % t50_lat_ma = t50_lat_ma
           t % qhef_ver_m = qhef_ver_m
           t % qhef_lat_m = qhef_lat_m
           t % Lreach = Lreach
           t % leng_int = leng_int
         nullify (t % left)
         nullify (t % right)
      ! Otherwise, insert into correct subtree
      else if (FID_reach < t % FID_reach) then
         call insert (t % left,w_ma,t50_ver_ma,t50_lat_ma, &
      qhef_ver_m, qhef_lat_m, Lreach, FID_reach, leng_int)
      else
         call insert (t % right,w_ma,t50_ver_ma,t50_lat_ma, &
      qhef_ver_m, qhef_lat_m, Lreach, FID_reach, leng_int)
      end if

   end subroutine insert

   recursive subroutine swat_insert (t,FID_reach, leng_int)

      type (swat_node_fld), pointer :: t  ! A tree
      integer, intent(in) :: FID_reach
      real*8, intent(in) :: leng_int
      ! If (sub)tree is empty, put number at root
      if (.not. associated (t)) then
         allocate (t)
         t % FID_reach = FID_reach
         t % leng_int = leng_int
         nullify (t % left)
         nullify (t % right)
      ! Otherwise, insert into correct subtree
      else if (FID_reach < t % FID_reach) then
         call swat_insert (t % left, FID_reach, leng_int)
      else
         call swat_insert (t % right, FID_reach, leng_int)
      end if

   end subroutine swat_insert

   recursive subroutine print_tree (t)
   ! Print tree in infix order

      type (nxs_node_fld), pointer :: t  ! A tree

      if (associated (t)) then
         call print_tree (t % left )
         print *, t % FID_reach
         call print_tree (t % right )
      end if

   end subroutine print_tree

   recursive subroutine print_length_tree (t)
   ! Print tree in infix order
      type (nxs_node_fld), pointer :: t  ! A tree

      if (associated (t)) then
         call print_length_tree (t % left )
!         rate_tmp = rate_tmp + t%qhef_ver_m + t%qhef_lat_m
         length_tmp = length_tmp + t% leng_int 
         call print_length_tree (t % right )
      end if

   end subroutine print_length_tree

   recursive subroutine print_frac_tree (t,tlength)
   ! Print tree in infix order
      type (nxs_node_fld), pointer :: t  ! A tree
      real*8, intent(in) :: tlength
      if (associated (t)) then
         call print_frac_tree (t % left, tlength )
         t%frac =  t%leng_int / tlength 
!         print *,'tfrac',t%frac
         call print_frac_tree (t % right, tlength )
      end if

   end subroutine print_frac_tree
!
   recursive subroutine print_swat_tree (t)
   ! Print tree in infix order

      type (swat_node_fld), pointer :: t  ! A tree

      if (associated (t)) then
         call print_swat_tree (t % left )
         print *, t % FID_reach
         call print_swat_tree (t % right )
      end if

   end subroutine print_swat_tree
!
   recursive subroutine print_length_swat_tree (t)
   ! Print tree in infix order
      type (swat_node_fld), pointer :: t  ! A tree

      if (associated (t)) then
         call print_length_swat_tree (t % left )
         length_tmp = length_tmp + t% leng_int 
         call print_length_swat_tree (t % right )
      end if

   end subroutine print_length_swat_tree
!
   recursive subroutine nxs_to_swat (t,i,s_length,q,tqv,tql,ich,xsarea)
   ! Print tree in infix order
      real*8, intent(inout) :: q(2)
      real*8, intent(inout) :: tqv
      real*8, intent(inout) :: tql
      integer, intent(inout) :: ich
      real*8, intent(inout) :: xsarea

      type (swat_node_fld), pointer :: t  ! A tree
      real*8 :: ratio ! total length ratio between swat and nexss 
      real*8 :: s_length ! total swat reach length
      real*8 :: add(2)
      integer :: grid_id
      integer :: i
      integer :: swat_grid_id
      if (associated (t)) then
!         do i=1,nxs_rch
         call nxs_to_swat (t % left, i, s_length, q ,tqv, tql, ich,xsarea)
           if(rch(i)%fnd == 0) then
             call get_nex_grid_id(rch(i)%t,t%FID_reach,grid_id)
             ratio = rch(i)%tlength / s_length
!             if(s_length < rch(i)%tlength) ratio = 1.d0 /ratio
!             if(ratio > 10000.d0) return ! swat reach not in nhd
             if(t%FID_reach == grid_id) then  !nxs id = swat id
                 overlap = overlap + 1
!print *,'nxrch,swrc-',i,ich,overlap,grid_id
             endif
!             if(t%FID_reach == grid_id .and. overlap > 2 ) then
             if(t%FID_reach == grid_id .and. overlap > 2 .and. ratio > 1.d-8 ) then
               rch(i)%rate(:) = max(0.d0,rch(i)%rate(:))
               if(ratio > 1.d0) then
                  add(1:2) = rch(i)%rate(1:2)*rch(i)%area*s_length/rch(i)%tlength
               else
                  add(1:2) = rch(i)%rate(1:2)*rch(i)%area*rch(i)%tlength/s_length
               endif
!               if(ratio > 1.d0) then
!                 q(1:2) = q(1:2) + add(1:2)
!               else
                 q(1:2) = q(1:2) + add(1:2)
!               endif
!               q(1:2) = rch(i)%rate(1:2)
               tqv = tqv + add(1)*rch(i)%t%t50_ver_ma
               tql = tql + add(2)*rch(i)%t%t50_lat_ma
               xsarea = xsarea + (add(1)+add(2))*rch(i)%xsarea
!               rch(i)%fnd = 1
               ich = i 
               return
             endif
           endif
         call nxs_to_swat (t % right, i, s_length, q, tqv, tql ,ich,xsarea)
!         enddo
!         print *, t % FID_reach
      end if

   end subroutine nxs_to_swat
!
   recursive subroutine get_nex_grid_id (t,swat_grid_id,grid_id)
   ! get grid id of nex segment

      type (nxs_node_fld), pointer :: t  ! A tree
      integer :: grid_id
      integer :: swat_grid_id

      if (associated (t)) then
         call get_nex_grid_id (t % left, swat_grid_id, grid_id )
         grid_id = t % FID_reach
         if(grid_id == swat_grid_id) return
         call get_nex_grid_id (t % right, swat_grid_id, grid_id )
      end if

   end subroutine get_nex_grid_id
!
!   recursive subroutine get_nex_flux (t,q)
   ! get grid id of nex segment

!      type (nxs_node_fld), pointer :: t  ! A tree
!      integer :: grid_id

!      if (associated (t)) then
!         call get_nex_grid_id (t % left )
!         q = t % FID_reach
!         call get_nex_grid_id (t % right )
!      end if

!   end subroutine get_nex_flux
!
end program tree_sort
