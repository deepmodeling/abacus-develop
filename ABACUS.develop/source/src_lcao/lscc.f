      subroutine sphbsl(n,r,A,val) 
        integer :: n
        real*8 :: r,A
        real*8 :: x,val
        x = r*A
        if (n .eq. 0) then
        
          if ( x .lt. 1.d-3 ) then
            val = 1 + x**2/6
          else
            val = dsinh(x)/x
          end if
        else if (n .eq. 2) then
        
          if ( x .lt. 1.d-2 ) then
            val = -x**2/15 -x**4/210 - x**6/7560
          else
            val = 3*dcosh(x)/x**2 + (-3-x**2)*dsinh(x)/x**3
          end if
        
        else if (n .eq. 4) then

          if( x .lt. 5.d-1)then
            val = x**4/945 + x**6/20790 + x**8/1081080 + x**10/97297200
          else
            val = -5*(21+2*x**2)*dcosh(x)/x**4+(105+45*x**2+x**4)*
     &       dsinh(x)/x**5
          end if
        
        else if (n .eq. 6) then
        
          if ( x .lt. 9.d-1) then
            val = -x**6/135135-x**8/4054050-x**10/275675400
          else
            val = 21*(495+60*x**2+x**4)*dcosh(x)/x**6 +
     &       (-10395-4725*x**2-210*x**4-x**6)*dsinh(x)/x**7
          end if
        
        else
        end if
      END subroutine sphbsl

      subroutine sphhnk(n,r,A,val)
        integer :: n
        real*8 :: r,A
        real*8 :: x,val
        x = r*A
        if (n .eq. 0) then
        
          if ( x .lt. 1.d-3 ) then
            val = -1/x + 1 -x/2 + x**2/6
          else
            val = -dexp(-x)/x
          endif
        
        else if (n .eq. 2) then
        
          if ( x .lt. 1.d-2) then
            val = 3/x**3-1/(2*x)+x/8-x**2/15+x**3/48
          else
            val = dexp(-x)*(3+3*x+x**2)/x**3
          endif
        
        else if (n .eq. 4) then
        
          if (x .lt. 5.d-1) then
            val = -105/x**5 + 15/(2*x**3) - 3/(8*x) + x/48 - x**3/384 
     &        +x**4/945
          else
            val = -dexp(-x)*(105+105*x+45*x**2+10*x**3+x**4)/x**5
          endif
        
        else if (n .eq. 6) then

          if (x .lt. 9.d-1) then
            val = 10395/x**7 - 945/(2*x**5) + 105/(8*x**3) -5/(16*x) + 
     &            x/128-x**3/3840 + x**5/46080 - x**6/135135
          else
            val = dexp(-x)*(10395+10395*x+4725*x**2+1260*x**3+210*x**
     &       4+21*x**5+x**6)/x**7
          endif
        
        else
        endif
      END SUBROUTINE sphhnk


