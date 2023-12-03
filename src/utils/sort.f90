 subroutine quick_sort(x,iorder,isize)                        

  implicit none                                                  

  integer,intent(in)             :: isize                        
  double precision,intent(inout) :: x(isize)          
  integer,intent(inout)          :: iorder(isize)                

  call rec_quicksort(x,iorder,isize,1,isize,1)           

end subroutine

recursive subroutine rec_quicksort(x,iorder,isize,first,last,level) 

  implicit none                                                  

  integer, intent(in)            :: isize, first, last, level    
  integer,intent(inout)          :: iorder(isize)                
  double precision, intent(inout):: x(isize)          
  double precision               :: c, tmp            
  integer                        :: itmp                         
  integer                        :: i, j                         

  c = x( shiftr(first+last,1) )                                  
  i = first                                                      
  j = last                                                       
  do                                                             
    do while (x(i) < c)                                          
      i=i+1                                                      
    end do                                                       
    do while (c < x(j))                                          
      j=j-1                                                      
    end do                                                       
    if (i >= j) then                                             
 exit                                                            
  end if                                                          
    tmp  = x(i)                                                  
    x(i) = x(j)                                                  
    x(j) = tmp                                                   
    itmp      = iorder(i)                                        
    iorder(i) = iorder(j)                                        
    iorder(j) = itmp                                             
    i=i+1                                                        
    j=j-1                                                        
  end do
  if ( ((i-first <= 10000).and.(last-j <= 10000)).or.(level<=0) ) then
    if (first < i-1) then                                        
      call rec_quicksort(x, iorder, isize, first, i-1,level/2) 
    end if                                                        
    if (j+1 < last) then                                         
      call rec_quicksort(x, iorder, isize, j+1, last,level/2)  
    end if                                                        
  else                                                           
    if (first < i-1) then                                        
      call rec_quicksort(x, iorder, isize, first, i-1,level/2) 
    end if                                                        
    if (j+1 < last) then                                         
      call rec_quicksort(x, iorder, isize, j+1, last,level/2)  
    end if                                                        
  end if                                                          
end subroutine

subroutine set_order(x,iorder,isize,jsize)

  implicit none                                                  
  
  integer                        :: isize,jsize
  double precision               :: x(isize,jsize)              
  double precision,allocatable   :: xtmp(:,:)           
  integer                        :: iorder(*)                    
  integer                        :: i,j

  allocate(xtmp(isize,jsize)) 

  do i=1,isize                                                   
    do j=1,jsize                                                   
      xtmp(i,j) = x(i,iorder(j)) 
    end do                                                          
  end do

  do i=1,isize                                                   
    do j=1,jsize                                                   
      x(i,j) = xtmp(i,j)
    end do                                                          
  end do                                                          

  deallocate(xtmp)                                               

end subroutine 
