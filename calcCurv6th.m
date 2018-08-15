%Estimate y(x) and calculate Curvature from a 2D point Cloud. We cut the
%cloud in proper regions and then estimate each of them using a 6th order
%polynomial.

%****** INPUT PARAMETERS ******

%vectors x=1xN and y=1xN: x and y coordinates of the cloud

%ds: the length that we move on the curve to define the next point to test 

%dsError: the % error that we allow

%masDist: the maximum Distance that we allow the next estimated point to
%have. You can diable it by setting it to inf

%maxPixelEfthia: the maximum pixel that we allow to have the same coordinate
%Y without beiing considered as a vertical line

%maxRegion: the maximum regions that we allow to exist (in pixel number)


clearvars -except x y;

%filter points 
samples=x;
samples2=x(2:end);
samples2=samples2-samples(1:end-1);
samples2=samples2(2:end)-samples2(1:end-1);
while find( abs(samples2)==2 )
    x(find( abs(samples2)==2,1 )+1)=[];
    y(find( abs(samples2)==2,1 )+1)=[];
    samples=x;
    samples2=x(2:end);
    samples2=samples2-samples(1:end-1);
    samples2=samples2(2:end)-samples2(1:end-1);
end

clear samples samples2


ds=0.1;
dserror=ds*0.001; 
if ds<1
    pSteps=1;
else 
    pSteps=int16(ds);
end 
maxDist=2;
maxPixelEfthia=15;
maxRegion=15;



%cut the curve 
pCutPoint=0;
samples=x(1);
pCut=1;
i=2;
p=2;
zerosInRow=0;
caseNum=0;


while p<=length(x)
    
    samples(p-pCutPoint)=x(p);
    samples2=samples(2:end);
    samples2=samples2-samples(1:end-1);
    
    if ~isempty(samples2) && (samples2(end)~=0)
        
        if zerosInRow>=maxPixelEfthia-1
            caseNum=1;
        else
            zerosInRow=0;
        end
        
    elseif ~isempty(samples2)
        zerosInRow=zerosInRow+1;
    end
    
    if caseNum==1
        caseNum=0;
        pCut(i)=p-zerosInRow-1;
        zerosInRow=0;
        i=i+1;
        pCut(i)=p; 
        p=p-1;
        pCutPoint=p;
        clear samples samples2       
        i=i+1;     
    elseif sum(samples2>0) && sum(samples2<0)         
        pCut(i)=p;
        p=p-1;
        pCutPoint=p;
        clear samples samples2       
        i=i+1;        
    end
    
    p=p+1;
    
end
if zerosInRow>=maxPixelEfthia
    pCut(end+1)=p-zerosInRow-1;
end    
%pCut( [false pCut(1:end-1)-pCut(2:end)==0] )=[];
%pCut( find( (pCut(1:end-1)-pCut(2:end))==0) )=[];

pCut2=pCut(2:end);
pCut2=pCut2-pCut(1:end-1);
x(pCut(pCut2==1))=[];
y(pCut(pCut2==1))=[];
pCut3=find(pCut2==1);
while ~isempty(pCut3)   
    pCut(pCut3(1)+1:end)=pCut(pCut3(1)+1:end)-1;
    pCut3(1)=[];
end
pCut(pCut2==1)=[];
pCut(end+1)=length(x);
pCut=unique(pCut);

%find max Regions
pCut2=pCut(2:end);
pCut2=pCut2-pCut(1:end-1);
pCut3=find(pCut2>maxRegion,1);
while pCut3
    pnew=pCut(pCut3)+floor( (pCut(pCut3+1)-pCut(pCut3))/2 );
    pCut(end+1)=pnew;
    pCut=sort(pCut);
    pCut2=pCut(2:end);
    pCut2=pCut2-pCut(1:end-1);
    pCut3=find(pCut2>maxRegion,1);
end



%Finish Cut

%clear i p samples samples2 pCutPoint
i=1;
miki_S(i)=0;
caseNum=1;

p=pCut(1);
p_old=p;     
xnew=x(p);
ynew=y(p);

k=1;
while p<length(x)-2
    
    if k+1<=length(pCut)
        if p>=pCut(k) 
            k=sum(p>=pCut);
            samples=[ (x (pCut(k):pCut(k+1)-1) )', ( y(pCut(k):pCut(k+1)-1) )' ]; 
            
            xnew=x(p); 
            ynew=y(p);           
            sampleSum=sum(samples(2:end,1)-samples(1:end-1,1));
            

            if sampleSum~=0 
                
                if sampleSum>0
                    caseNum=1; %x++
                else
                    caseNum=2; %x--
                end
                
                lseA=[ samples(:,1).^6 ,samples(:,1).^5 ,samples(:,1).^4 ,samples(:,1).^3 ,samples(:,1).^2 , samples(:,1) , ones( length(samples(:,1)),1 ) ];
                lseY=samples(:,2);
                lseB=lseA \ lseY;
                %lseB=(lseA'*lseA)\lseA'*lseY;
                if any(isnan(lseB)) || any(isinf(lseB))
                     lseB=pinv(lseA)*lseY;   
                end
                if any( abs(lseB>10^6))
                    if sum(samples(2:end,2)-samples(1:end-1,2))<0
                        caseNum=3;
                    else 
                        caseNum=4;
                    end                       
                    
                end              
                                
            else
                if samples(1,2)>samples(2,2)
                    caseNum=3; %straight line up
                else
                    caseNum=4;%straight line down
                end
            end
            
            k=k+1;
        end
    end
    
    xEkt(i)=xnew;
    
   
    if caseNum==1        
        yEkt(i)=lseB(1)*(xnew^6)+lseB(2)*(xnew^5)+lseB(3)*(xnew^4)+lseB(4)*(xnew^3)+lseB(5)*(xnew^2)+lseB(6)*xnew+lseB(7);
        kamp(i)=abs( 30*lseB(1)*xnew^4+20*lseB(2)*xnew^3+12*lseB(3)*xnew^2+6*lseB(4)*xnew+2*lseB(5) )/(  1 + ( 6*lseB(1)*xnew^5+5*lseB(2)*xnew^4+4*lseB(3)*xnew^3+3*lseB(4)*xnew^2+2*lseB(5)*xnew+lseB(6) )^2  )^(3/2) ;

        x1=xnew;
        x2=xnew+ds;
        y1=sqrt( 1+( 6*lseB(1)*x1^5+5*lseB(2)*x1^4+4*lseB(3)*x1^3+3*lseB(4)*x1^2+2*lseB(5)*x1+lseB(6) )^2 );
        y2=sqrt( 1+( 6*lseB(1)*x2^5+5*lseB(2)*x2^4+4*lseB(3)*x2^3+3*lseB(4)*x2^2+2*lseB(5)*x2+lseB(6) )^2 );

        S=(y1+y2)*(x2-x1)/2;

        while S>ds+dserror || S<ds-dserror
            xtemp=(x2+x1)/2;
            y2=sqrt( 1+( 6*lseB(1)*xtemp^5+5*lseB(2)*xtemp^4+4*lseB(3)*xtemp^3+3*lseB(4)*xtemp^2+2*lseB(5)*xtemp+lseB(6) )^2 );
            S=(y1+y2)*(xtemp-xnew)/2;
            if (S<ds-dserror)
                x1=xtemp;
            else 
                x2=xtemp;
            end
        end
        
        xnew=x2;
        ynew=lseB(1)*(x2^6)+lseB(2)*(x2^5)+lseB(3)*(x2^4)+lseB(4)*(x2^3)+lseB(5)*(x2^2)+lseB(6)*x2+lseB(7);
   
    elseif caseNum==2
        yEkt(i)=lseB(1)*(xnew^6)+lseB(2)*(xnew^5)+lseB(3)*(xnew^4)+lseB(4)*(xnew^3)+lseB(5)*(xnew^2)+lseB(6)*xnew+lseB(7);
        kamp(i)=abs( 30*lseB(1)*xnew^4+20*lseB(2)*xnew^3+12*lseB(3)*xnew^2+6*lseB(4)*xnew+2*lseB(5) )/(  1 + ( 6*lseB(1)*xnew^5+5*lseB(2)*xnew^4+4*lseB(3)*xnew^3+3*lseB(4)*xnew^2+2*lseB(5)*xnew+lseB(6) )^2  )^(3/2) ;
        
        
        x1=xnew;
        x2=xnew-ds;
        y1=sqrt( 1+( 6*lseB(1)*x1^5+5*lseB(2)*x1^4+4*lseB(3)*x1^3+3*lseB(4)*x1^2+2*lseB(5)*x1+lseB(6) )^2 );
        y2=sqrt( 1+( 6*lseB(1)*x2^5+5*lseB(2)*x2^4+4*lseB(3)*x2^3+3*lseB(4)*x2^2+2*lseB(5)*x2+lseB(6) )^2 );
        S=(y1+y2)*(x1-x2)/2;
        while S>ds+dserror || S<ds-dserror
            xtemp=(x2+x1)/2;
            y2=sqrt( 1+( 6*lseB(1)*xtemp^5+5*lseB(2)*xtemp^4+4*lseB(3)*xtemp^3+3*lseB(4)*xtemp^2+2*lseB(5)*xtemp+lseB(6) )^2 );
            S=(y1+y2)*(xnew-xtemp)/2;
            if (S<ds-dserror)
                x1=xtemp;
            else
                x2=xtemp;
            end
        end
        
        xnew=x2;
        ynew=lseB(1)*(x2^6)+lseB(2)*(x2^5)+lseB(3)*(x2^4)+lseB(4)*(x2^3)+lseB(5)*(x2^2)+lseB(6)*x2+lseB(7);
        
    elseif caseNum==3
        yEkt(i)=ynew;
        ynew=ynew-ds;
        kamp(i)=10^-5;
        S=ds;
   
    else %caseNum==4
        yEkt(i)=ynew;
        ynew=ynew+ds;
        kamp(i)=10^-5;
        S=ds;    
    end
    
    apostaseis=(pdist2( [xnew ynew], [x' y'])) ;
    [apos , p]=min(apostaseis);
    

      
      if p<p_old
          p=p_old+pSteps;
          if apos>maxDist
              xnew=x(p);ynew=y(p);
          end
      elseif p==p_old
          if apos>maxDist
              p=p_old+pSteps;
              xnew=x(p);ynew=y(p);
          end
      end
      
      p_old=p;
      i=i+1;
      miki_S(i)=miki_S(i-1)+S;
end
miki_S=miki_S(1:end-1);


kamp(kamp>5*mean(kamp))=mean(kamp);


%Plots
figure
plot(xEkt,yEkt,'r','Linewidth',2)
hold on
plot(x,y,'o')
plot(x(pCut),y(pCut),'x r','MarkerSize',10,'Linewidth',2)
plot(xEkt(1),yEkt(1),'x g','MarkerSize',10,'Linewidth',2)
daspect([1 1 1])
hold off
legend('yEkt', 'y')
figure
plot(miki_S,kamp,'o-')
legend('ektimomeni kampilotita')
%}
