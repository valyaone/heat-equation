function C = progonka(a,b,c,F)
    N = size(F,2)-1;
    alpha=zeros(N+1);
    betta=zeros(N+1);
    alpha(1)=b(1)/c(1);
    betta(1)=F(1)/c(1);
    
    for i=2:N
        alpha(i)=b(i)/(c(i)-a(i)*alpha(i-1)); % матрица значение alphai
        betta(i)=(F(i)+a(i)*betta(i-1))/(c(i)-a(i)*alpha(i-1));
    end
    
    betta(N+1)=(F(N+1)+a(N+1)*betta(N))/(c(N+1)-a(N+1)*alpha(N));
    y(N+1)=betta(N+1);
    
    for i=1:N
        y(N-i+1)=alpha(N-i+1)*y(N-i+2)+betta(N-i+1); 
    end
    
    for i=1:N+1   
        C(i)=y(i);
    end