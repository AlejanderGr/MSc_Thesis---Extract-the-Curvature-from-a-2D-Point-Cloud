%Given the vectors x and y we extract the variable S which is the 
%length-distance of every pixel starting from the beginning. Then we 
%estimate the parametrcical functions x(s) and y(s) and calculate the 
%curvature of the point-cloud. We make new estimation in every step using
%LSE and 2nd order polynomials

%****** INPUT PARAMETERS ******

%vectors x=Nx1 and y=Nx1 : x and y coordinates of the cloud

%ds: the length in pixels that we move on the curve to define the next 
%point to test 

%lseN: the number of the points prior and after the currert point that we
%use for the LSE 

clearvars -except x y ;

ds=1;
lseN=10;


x=[x(end-lseN:end-1); x; x(2:1+lseN)];
y=[y(end-lseN:end-1); y; y(2:1+lseN)];


S=zeros(length(x),1);
for i=1:length(x)-1
    S(i+1)=S(i)+pdist2( [x(i),y(i)], [x(i+1) y(i+1)]);
end

p=lseN+1;
s=S(p);
miki_s=s;
i=1;
while p<=length(S)-lseN+1
    
    sampleS=S(p-lseN:p+lseN-1);
    sampleX=x(p-lseN:p+lseN-1);
    sampleY=y(p-lseN:p+lseN-1); 
    
    matA=[sampleS.^2 , sampleS, ones( length(sampleS),1 ) ];
    matBx=matA \ sampleX;
    matBy=matA \ sampleY; 
    
    
    xEkt(i)=matBx(1)*(s^2) + matBx(2)*s + matBx(3); 
    yEkt(i)=matBy(1)*(s^2) + matBy(2)*s + matBy(3);
        
    dx=2*matBx(1)*s + matBx(2);
    dy=2*matBy(1)*s + matBy(2);
    
    dx2=2*matBx(1);
    dy2=2*matBy(1); 
    
    %kamp(i)= abs( dx.*dy2-dx2.*dy ) ./ ( dx.^2 + dy.^2).^(3/2);
    kamp(i)=sqrt( (dx2)^2+(dy2)^2 );
    
    i=i+1;
    s=s+ds;
    p=find(S>s,1); 
    miki_s(i)=s;
    
    
end
miki_s=miki_s(1:end-1)-S(lseN+1);

S=S(lseN+1:end-lseN)-S(lseN+1);
x=x(lseN+1:end-lseN);
y=y(lseN+1:end-lseN);

%PLOT
figure
plot(x,y,'o')
hold on
plot(xEkt,yEkt,'r','Linewidth',2)
daspect([1 1 1])
legend('points','ektimisi')
plot(xEkt(1),yEkt(1),'x g','MarkerSize',10,'Linewidth',2)

figure
plot(S,x,'o')
hold on
plot(miki_s,xEkt,'r')
legend('x(s)','xEkt(s)')

figure
plot(S,y,'o')
hold on
plot(miki_s,yEkt,'r')
legend('y(s)','yEkt(s)')
hold off
 

figure
plot(kamp,'-')
legend('ektimomeni kampilotita')

%FILTER
windowSize = 9;
b = (1/windowSize)*ones(1,windowSize);
a = 1;
kamp2=filtfilt(b,a,kamp);
hold on
plot(kamp2,'o-g')
legend('ektimomeni kampilotita','filtrarismeni kampilotita')
hold off


