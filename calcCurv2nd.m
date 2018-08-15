%Estimate y(x) and calculate Curvature from a 2D point Cloud. We do a new 
% estimation in every step using the Least Square Method and 2nd order
% polynomials

%****** INPUT PARAMETERS ******

%vectors x=1xN and y=1xN: x and y coordinates of the cloud

%ds: the length in pixels that we move on the curve to define the next 
%point to test. Minimum value should be 2

%dsError: the % error that we allow

%lseN: the number of the points prior and after the currert point that we
%use for the LSE 

%maxSteps: the maximum increase of the index of the vector x that is
%allowed. You have to adjust it to the problem or disable it setting it inf

%masDist: the maximum Distance that we allow the next estimated point to
%have. You can disable it by setting it to inf


clearvars -except x y;

ds=2;
dserror=ds*0.001; 

lseN=10; 
maxSteps=6; 
            
maxDist=3*ds;
pSteps=int16(ds);


i=1;
miki_S(i)=0;

p=lseN+1;  
p_old=p;     
xnew=x(p);
ynew=y(p);


while p<length(x)-lseN+1  
    pis(i)=p;

    samples=[ (x(p-lseN:p+lseN))', (y(p-lseN:p+lseN))' ]; 
    samples2=samples(2:end,1);
    samples2=samples2-samples(1:end-1,1);
    
    if ( sum(samples2>0) && sum(samples2<0) )|| sum(samples2==0)>8
       
        xRot=samples(:,2);
        yRot=-samples(:,1);
        samples(:,1)=xRot;
        samples(:,2)=yRot;   
       
        %CS2
        lseA=[ samples(:,1).^2 , samples(:,1) , ones( length(samples(:,1)),1 ) ];
        lseY=samples(:,2);
        lseB=lseA \ lseY;
        %lseB=(lseA'*lseA)^-1*(lseA'*lseY);
        if any(isnan(lseB)) || any(isinf(lseB))
            lseB=pinv(lseA)*lseY;   
        end 
        
        %rotation 90 degrees
        xnewCS2=ynew;
        %ynewCS2=-xnew;
        yEktCS2=lseB(1)*(xnewCS2^2)+lseB(2)*xnewCS2+lseB(3);
        
        %rot -90      
        xEkt(i)=-yEktCS2;
        yEkt(i)=ynew;      


        kamp(i)=abs( 2*lseB(1) )/(  1 + ( 2*lseB(1)*xnewCS2+lseB(2) )^2  )^(3/2) ;


        A=4*lseB(1)^2;
        B=4*lseB(1)*lseB(2);
        C=1+lseB(2)^2;

        if A<10^-10
            A=10^-10;
        end

        x1=xnewCS2;
        x2=xnewCS2+ds;
        Sx1=(2*A*x1+B) * sqrt(A*x1^2+B*x1+C)/(4*A) +((4*A*C-B^2)/(8*A*sqrt(A))) *log( 2*A*x1+B+2*sqrt(A)*sqrt(A*x1^2+B*x1+C) );
        Sx2=(2*A*x2+B) * sqrt(A*x2^2+B*x2+C)/(4*A) + ((4*A*C-B^2)/(8*A*sqrt(A))) *log( 2*A*x2+B+2*sqrt(A)*sqrt(A*x2^2+B*x2+C) );

        S=Sx2-Sx1;

        while S>ds+dserror || S<ds-dserror
            xtemp=(x2+x1)/2;
            Sx2=(2*A*xtemp+B) * sqrt(A*xtemp^2+B*xtemp+C)/(4*A) +((4*A*C-B^2)/(8*A*sqrt(A))) *log( 2*A*xtemp+B+2*sqrt(A)*sqrt(A*xtemp^2+B*xtemp+C) );
            S=Sx2-Sx1;
            if (S<ds-dserror)
                x1=xtemp;
            else 
                x2=xtemp;
            end

        end

        xnewTemp=x2;
        ynewTemp=lseB(1)*(xnewTemp^2)+lseB(2)*xnewTemp+lseB(3);

        xRot=-ynewTemp;
        yRot=xnewTemp;
        xnewTemp=xRot;
        ynewTemp=yRot;

        apostaseis=(pdist2( [xnewTemp ynewTemp], [x' y'])) ;
        [apos , p]=min(apostaseis);


      
        if p>p_old+maxSteps  

            while p<=p_old || p>p_old+maxSteps
                apostaseis(p)=inf;
               [apos ,p]=min(apostaseis);
            end
            if apos>maxDist
                    p=p_old+pSteps;
                    S=pdist2([x(p) y(p)],[xnew ynew] );
                    xnew=x(p);
                    ynew=y(p);
                else
                    xnew=xnewTemp;
                    ynew=ynewTemp;
            end

        elseif p<=p_old
            x1=xnewCS2;
            x2=xnewCS2-ds;
            Sx2=(2*A*x2+B) * sqrt(A*x2^2+B*x2+C)/(4*A) + ((4*A*C-B^2)/(8*A*sqrt(A))) *log( 2*A*x2+B+2*sqrt(A)*sqrt(A*x2^2+B*x2+C) );
            S=abs(Sx1-Sx2);
            while S>ds+dserror || S<ds-dserror
                    xtemp=(x2+x1)/2;
                    Sx2=(2*A*xtemp+B) * sqrt(A*xtemp^2+B*xtemp+C)/(4*A) +((4*A*C-B^2)/(8*A*sqrt(A))) *log( 2*A*xtemp+B+2*sqrt(A)*sqrt(A*xtemp^2+B*xtemp+C) );
                    S=abs(Sx1-Sx2);
                    if (S<ds-dserror)
                        x1=xtemp;
                    else 
                        x2=xtemp;
                    end
            end
            xnewTemp=x2;
            ynewTemp=lseB(1)*(xnewTemp^2)+lseB(2)*xnewTemp+lseB(3);
            
            xRot=-ynewTemp;
            yRot=xnewTemp;
            xnewTemp=xRot;
            ynewTemp=yRot;      
            
            apostaseis=(pdist2( [xnewTemp ynewTemp], [x' y'])) ;
            [apos , p]=min(apostaseis);

            if p<=p_old
                p=p_old+pSteps; 
                S=pdist2([x(p) y(p)],[xnew ynew] );
                xnew=x(p);
                ynew=y(p);
            else
                while p<=p_old || p>p_old+maxSteps
                    apostaseis(p)=inf;
                    [apos ,p]=min(apostaseis);
                end
                if apos>maxDist
                    p=p_old+pSteps;
                    S=pdist2([x(p) y(p)],[xnew ynew] );                
                    xnew=x(p);
                    ynew=y(p);
                else
                    xnew=xnewTemp;
                    ynew=ynewTemp;
                end
            end

        else
            if apos>maxDist
                p=p_old+pSteps;
                S=pdist2([x(p) y(p)],[xnew ynew] );
                xnew=x(p);
                ynew=y(p);
            else
                xnew=xnewTemp;
                ynew=ynewTemp;
            end

        end
        
        
    
        
    else
        xEkt(i)=xnew;
        lseA=[ samples(:,1).^2 , samples(:,1) , ones( length(samples(:,1)),1 ) ];
        lseY=samples(:,2);
        lseB=lseA \ lseY;
        %lseB=(lseA'*lseA)^-1*(lseA'*lseY);
        if any(isnan(lseB)) || any(isinf(lseB))
            lseB=pinv(lseA)*lseY;   
        end 

        %y=a*x^2+b*x+c------a=lseB(1), b=lseB(2), c=lseB(3) 
        yEkt(i)=lseB(1)*(xnew^2)+lseB(2)*xnew+lseB(3);

        kamp(i)=abs( 2*lseB(1) )/(  1 + ( 2*lseB(1)*xnew+lseB(2) )^2  )^(3/2) ;


        A=4*lseB(1)^2;
        B=4*lseB(1)*lseB(2);
        C=1+lseB(2)^2;

        if A<10^-10
            A=10^-10;
        end

        x1=xnew;
        x2=xnew+ds;
        Sx1=(2*A*x1+B) * sqrt(A*x1^2+B*x1+C)/(4*A) +((4*A*C-B^2)/(8*A*sqrt(A))) *log( 2*A*x1+B+2*sqrt(A)*sqrt(A*x1^2+B*x1+C) );
        Sx2=(2*A*x2+B) * sqrt(A*x2^2+B*x2+C)/(4*A) + ((4*A*C-B^2)/(8*A*sqrt(A))) *log( 2*A*x2+B+2*sqrt(A)*sqrt(A*x2^2+B*x2+C) );

        S=Sx2-Sx1;

        while S>ds+dserror || S<ds-dserror
            xtemp=(x2+x1)/2;
            Sx2=(2*A*xtemp+B) * sqrt(A*xtemp^2+B*xtemp+C)/(4*A) +((4*A*C-B^2)/(8*A*sqrt(A))) *log( 2*A*xtemp+B+2*sqrt(A)*sqrt(A*xtemp^2+B*xtemp+C) );
            S=Sx2-Sx1;
            if (S<ds-dserror)
                x1=xtemp;
            else 
                x2=xtemp;
            end

        end

        xnewTemp=x2;
        ynewTemp=lseB(1)*(xnewTemp^2)+lseB(2)*xnewTemp+lseB(3);

        apostaseis=(pdist2( [xnewTemp ynewTemp], [x' y'])) ;
        [apos , p]=min(apostaseis);

        if p>p_old+maxSteps  

            while p<=p_old || p>p_old+maxSteps
                apostaseis(p)=inf;
               [apos ,p]=min(apostaseis);
            end
            if apos>maxDist
                    p=p_old+pSteps;
                    S=pdist2([x(p) y(p)],[xnew ynew] );
                    xnew=x(p);
                    ynew=y(p);
                else
                    xnew=xnewTemp;
                    ynew=ynewTemp;
            end

        elseif p<=p_old
            x1=xnew;
            x2=xnew-ds;
            Sx2=(2*A*x2+B) * sqrt(A*x2^2+B*x2+C)/(4*A) + ((4*A*C-B^2)/(8*A*sqrt(A))) *log( 2*A*x2+B+2*sqrt(A)*sqrt(A*x2^2+B*x2+C) );
            S=abs(Sx1-Sx2);
            while S>ds+dserror || S<ds-dserror
                    xtemp=(x2+x1)/2;
                    Sx2=(2*A*xtemp+B) * sqrt(A*xtemp^2+B*xtemp+C)/(4*A) +((4*A*C-B^2)/(8*A*sqrt(A))) *log( 2*A*xtemp+B+2*sqrt(A)*sqrt(A*xtemp^2+B*xtemp+C) );
                    S=abs(Sx1-Sx2);
                    if (S<ds-dserror)
                        x1=xtemp;
                    else 
                        x2=xtemp;
                    end
            end
            xnewTemp=x2;
            ynewTemp=lseB(1)*(xnewTemp^2)+lseB(2)*xnewTemp+lseB(3);
            apostaseis=(pdist2( [xnewTemp ynewTemp], [x' y'])) ;
            [apos , p]=min(apostaseis);

            if p<=p_old
                p=p_old+pSteps; 
                S=pdist2([x(p) y(p)],[xnew ynew] );
                xnew=x(p);
                ynew=y(p);
            else
                while p<=p_old || p>p_old+maxSteps
                    apostaseis(p)=inf;
                    [apos ,p]=min(apostaseis);
                end
                if apos>maxDist
                    p=p_old+pSteps;
                    S=pdist2([x(p) y(p)],[xnew ynew] );                
                    xnew=x(p);
                    ynew=y(p);
                else
                    xnew=xnewTemp;
                    ynew=ynewTemp;
                end
            end

        else
            if apos>maxDist
                p=p_old+pSteps;
                S=pdist2([x(p) y(p)],[xnew ynew] );
                xnew=x(p);
                ynew=y(p);
            else
                xnew=xnewTemp;
                ynew=ynewTemp;
            end

        end
        
    end 

    p_old=p;
    i=i+1;
    miki_S(i)=miki_S(i-1)+S;
end
miki_S=miki_S(1:end-1);


%PLOT
figure
plot(x,y,'o')
hold on
plot(xEkt,yEkt,'r','Linewidth',2)
daspect([1 1 1])
legend('points','ektimisi')
plot(xEkt(1),yEkt(1),'x g','MarkerSize',10,'Linewidth',2)
hold off
daspect([1 1 1])
legend('yEkt', 'y')
figure
plot(miki_S,kamp,'r')


%filter
windowSize = 11;
b = (1/windowSize)*ones(1,windowSize);
a = 1;
kamp2=filtfilt(b,a,kamp);
hold on
plot(miki_S,kamp2,'o-g')
legend('ektimomeni kampilotita','filtrarismeni kampilotita')
hold off
