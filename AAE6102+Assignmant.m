
%AAE6102 Assignment
%Zhang Hengwei
tic
clc
clear
GM=3.986005e+14;     %万有引力常数与地球质量乘积 m3/s2
we=7.2921151467e-5;  %地球自转角速度 rad/s
c=299792458;         %光速
F = -4.442807633e-10; % Constant, [sec/(meter)^(1/2)]
PRN=[6,17,30,10,23,22,26,5];
load('C:\Users\18717\Desktop\Assignment\Data\rcvr.dat')  
for i=1:length(PRN)
    Prange(1,i)=abs(rcvr(i,3)); 
end
load('C:\Users\18717\Desktop\Assignment\Data\eph.dat')
eph=eph([2:end 1],:);            
for i=1:length(PRN)   
     e=eph(i,9);      idt=eph(i,17);
     af0=eph(i,5);    w=eph(i,13);     cic=eph(i,21);
     af1=eph(i,6);    cuc=eph(i,19);   cis=eph(i,20);
     af2=eph(i,7);    cus=eph(i,18);   om0=eph(i,14);
     M0=eph(i,12);    crc=eph(i,23);   om=eph(i,16);
     A=eph(i,10);     crs=eph(i,22);   toe=eph(i,4);
     dn=eph(i,11);    i0=eph(i,15);
         
      DT=440992-toe;
      tk=440992-toe-Prange(1,i)/c;
      %计算时刻平近点角M，偏近点角E，真近点角f
      M=M0+(sqrt(GM/(A^6))+dn)*tk;
      E0=0;
      E=M;
%      E = kepOrb2E(M,e);
      while (abs(E-E0)>1e-12)
             E0=E;
             E=M+e*sin(E0); 
      end
    %Compute relativistic correction term
    dtr = F * e * A * sin(E);
    RelDelay(i) = dtr*c;
    C_NA(i)=af0+af1*DT+af2*DT*DT+RelDelay(i)/c;
      v1=sqrt((1-e^2))*sin(E);
      v2=cos(E)-e;
      f=atan2(v1,v2);
      %计算卫星升交角，卫星地心距离和轨道面倾角
      u0=w+f;                                              
      u=u0+cuc*cos(2*u0)+cus*sin(2*u0);                   
      r=A^2*(1-e*cos(E))+crc*cos(2*u0)+crs*sin(2*u0);        
      is=i0+idt*tk+cic*cos(2*u0)+cis*sin(2*u0);             

      xr=r*cos(u);
      yr=r*sin(u);
      L=om0+(om-we)*tk-we*toe;  

      X(i)=xr*cos(L)-yr*cos(is)*sin(L);
      Y(i)=xr*sin(L)+yr*cos(is)*cos(L);
      Z(i)=yr*sin(is); 
end

UserPosition=[-2694685.473,-4293642.366,3857878.924];              
UserDeltaT=0;                       
V_N=length(PRN);                    
R=ones(1,V_N);
Wxyznew=ones(3,V_N);
RN=ones(1,V_N);
Error=0;
ComputeTime=0;
Wxyz=ones(V_N,3);
for i=1:V_N
Wxyz(i,1)=X(i);
Wxyz(i,2)=Y(i);
Wxyz(i,3)=Z(i);
end

while ComputeTime<10
       ComputeTime=ComputeTime+1;       
       for n=1:V_N
           R(1,n)=sqrt((Wxyz(n,:)-UserPosition)*(Wxyz(n,:)-UserPosition)');
       end 
       travelTime=zeros(1,V_N);                   
       for n=1:V_N   
           Wxyznew(:,n)=e_r_corr(travelTime(n), Wxyz(n,:)');
       end                             

       for n=1:V_N
           RN(1,n)=norm(Wxyznew(:,n)-UserPosition','fro');  
       end                
       DeltaP=RN-Prange-c*C_NA+c*UserDeltaT;     
       A=ones(V_N,3);
       for n=1:V_N
           A(n,:)=(Wxyznew(:,n)'-UserPosition)/RN(1,n);
       end  
       h=[A ones(V_N,1)];              
       DeltaX=(h'*h)\(h'*DeltaP');                
       TempDeltaX=DeltaX(1:3,:); 
       Error=max(abs(TempDeltaX));
       UserPosition=UserPosition+TempDeltaX';    
       UserDeltaT=UserDeltaT+DeltaX(4,1)/(-c);
end

