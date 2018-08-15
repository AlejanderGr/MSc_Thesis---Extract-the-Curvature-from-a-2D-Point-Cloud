%Given the vectors x and y we extract the variable S which is the 
%length-distance of every pixel starting from the beginning. Then we 
%estimate the parametrcical functions x(s) and y(s) and calculate the 
%curvature of the point-cloud using a single N order polynomial.

%****** INPUT PARAMETERS ******

%vectors x=Nx1 and y=Nx1 : x and y coordinates of the cloud

%ds: the length in pixels that we move on the curve to define the next 
%point to test 



clearvars -except x y ;
ds=1;

S=zeros(length(x),1);
for i=1:length(x)-1
    S(i+1)=S(i)+pdist2( [x(i),y(i)], [x(i+1) y(i+1)]);
end


%{
%13th
 matA=[S.^13, S.^12, S.^11, S.^10, S.^9, S.^8, S.^7, S.^6 ,S.^5 ,S.^4 ,S.^3 ,S.^2 , S, ones( length(S),1 ) ];
 matB=matA\x;
 xEkt=matB(1)*(s.^13)+matB(2)*(s.^12)+matB(3)*(s.^11)+matB(4)*(s.^10)+matB(5)*(s.^9)+matB(6)*(s.^8)+matB(7)*(s.^7)+matB(8)*(s.^6)+matB(9)*(s.^5)++matB(10)*(s.^4)+matB(11)*(s.^3)+matB(12)*(s.^2)+matB(13)*s+matB(14);
 plot(s,xEkt,'y')
 error13=mean((x-xEkt).^2);

%10th
 matA=[S.^10, S.^9, S.^8, S.^7, S.^6 ,S.^5 ,S.^4 ,S.^3 ,S.^2 , S, ones( length(S),1 ) ];
 matB=matA\x;
 s=S;
 xEkt=matB(1)*(s.^10)+matB(2)*(s.^9)+matB(3)*(s.^8)+matB(4)*(s.^7)+matB(5)*(s.^6)+matB(6)*(s.^5)+matB(7)*(s.^4)+matB(8)*(s.^3)+matB(9)*(s.^2)++matB(10)*s+matB(11);
 plot(s,xEkt,'r')
 error10=mean((x-xEkt).^2);
 
 %7th
 matA=[S.^7, S.^6 ,S.^5 ,S.^4 ,S.^3 ,S.^2 , S, ones( length(S),1 ) ];
 matB=matA\x;
 s=S;
 xEkt=matB(1)*(s.^7)+matB(2)*(s.^6)+matB(3)*(s.^5)+matB(4)*(s.^4)+matB(5)*(s.^3)+matB(6)*(s.^2)+matB(7)*s+matB(8);
 plot(s,xEkt,'g')
 error7=mean((x-xEkt).^2);

 %3rd
 matA=[S.^3 ,S.^2 , S, ones( length(S),1 ) ];
 matB=matA \ x;
 xEkt=matB(1)*(s.^3)+matB(2)*(s.^2)+matB(3)*(s.^1)+matB(4);
 plot(s,xEkt,'x g')
 error3=mean( (x-xEkt).^2 );

%6th
 matA=[ S.^6 ,S.^5 ,S.^4 ,S.^3 ,S.^2 , S, ones( length(S),1 ) ];
 matBx=matA \ x;
 xEkt=matBx(1)*(s.^6)+matBx(2)*(s.^5)+matBx(3)*(s.^4)+matBx(4)*(s.^3)+matBx(5)*(s.^2)+matBx(6)*s+matBx(7);
 plot(s,xEkt,'g') 
 error6=mean( (x-xEkt).^2 );

  %4th
 matA=[S.^4 ,S.^3 ,S.^2 , S, ones( length(S),1 ) ];
 matBx=matA \ x;
 xEkt=matBx(1)*(s.^4)+matBx(2)*(s.^3)+matBx(3)*(s.^2)+matBx(4)*s+matBx(5);
 plot(s,xEkt,'y')
 error4=mean( (x-xEkt).^2 );
%}

 %5th
 matA=[S.^5 ,S.^4 ,S.^3 ,S.^2 , S, ones( length(S),1 ) ];
 
 matBx=matA \ x;
 matBy=matA \ y;


s=S;
 xEkt=matBx(1)*(s.^5)+matBx(2)*(s.^4)+matBx(3)*(s.^3)+matBx(4)*(s.^2)+matBx(5)*s+matBx(6); 
 yEkt=matBy(1)*(s.^5)+matBy(2)*(s.^4)+matBy(3)*(s.^3)+matBy(4)*(s.^2)+matBy(5)*s+matBy(6);
 
figure
plot(x,y,'o')
daspect([1 1 1])
legend('object ( y(x) )')
 figure
 plot(S,x,'o')
 hold on
 plot(s,xEkt,'r')
 legend('x(s)','xEkt(s)')
 %daspect([1 1 1])
 figure
 plot(S,y,'o')
 hold on
 plot(s,yEkt,'r')
 %daspect([1 1 1])
  legend('y(s)','yEkt(s)')
 error5x=mean( (x-xEkt).^2 );
 error5y=mean( (y-yEkt).^2 );

 
 s=0:ds:S(end);
 xEkt=matBx(1)*(s.^5) + matBx(2)*(s.^4) + matBx(3)*(s.^3) + matBx(4)*(s.^2) + matBx(5)*s + matBx(6);
 yEkt=matBy(1)*(s.^5) + matBy(2)*(s.^4) + matBy(3)*(s.^3) + matBy(4)*(s.^2) + matBy(5)*s + matBy(6);
 
 dx=5*matBx(1)*s.^4 + 4*matBx(2)*s.^3 + 3*matBx(3)*s.^2 + 2*matBx(4)*s + matBx(5);
 dy=5*matBy(1)*s.^4 + 4*matBy(2)*s.^3 + 3*matBy(3)*s.^2 + 2*matBy(4)*s + matBy(5);

dx2=20*matBx(1)*s.^3 + 12*matBx(2)*s.^2 + 6*matBx(3)*s + 2*matBx(4);
dy2=20*matBy(1)*s.^3 + 12*matBy(2)*s.^2 + 6*matBy(3)*s + 2*matBy(4); 
kamp= abs( dx.*dy2-dx2.*dy ) ./ ( dx.^2 + dy.^2).^(3/2);
 

figure
plot(kamp,'o')
legend('ektimomeni kampilotita')
figure
plot(x,y,'o')
hold on
plot(xEkt,yEkt,'r','Linewidth',2)
legend('y(x), yEkt(x)')



