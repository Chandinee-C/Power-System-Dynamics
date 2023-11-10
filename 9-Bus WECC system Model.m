%----------------------------------------------------------------------%
%             Power System Dynamics                                    % 
%                                                                      %       
%                   7th Order full-netwrok Model                       %
%                                                                      %
%  -Chandinee C                                                        %
%                                                                      %
% Refernce: POWER SYSTEM DYNAMICS                                      % 
%               AND STABILITY                                          % 
%           Peter W. Sauer and M. A. Pa                                %                                     %
%----------------------------------------------------------------------%
clear all;
n=9; % No. of nodes
g=3; %no. of generator nodes
l=7; % no of non-generator or load nodes

% Z matrix
z= [0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i 0.0000 + 0.0576i   0.0000 + 0.0000i   0.0000 + 0.0000i 0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i 0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i  0.0000 + 0.0625i   0.0000 + 0.0000i   0.0000 + 0.0000i
   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i 0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i 0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0586i
   0.0000 + 0.0576i   0.0000 + 0.0000i   0.0000 + 0.0000i 0.0000 + 0.0000i   0.0100 + 0.0850i   0.0170 + 0.0920i 0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i 0.0100 + 0.0850i   0.0000 + 0.0000i   0.0000 + 0.0000i 0.0320 + 0.1610i   0.0000 + 0.0000i   0.0000 + 0.0000i
   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i 0.0170 + 0.0920i   0.0000 + 0.0000i   0.0000 + 0.0000i 0.0000 + 0.0000i   0.0000 + 0.0000i   0.0390 + 0.1700i
   0.0000 + 0.0000i   0.0000 + 0.0625i   0.0000 + 0.0000i 0.0000 + 0.0000i   0.0320 + 0.1610i   0.0000 + 0.0000i 0.0000 + 0.0000i   0.0085 + 0.0720i   0.0000 + 0.0000i
   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i 0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i 0.0085 + 0.0720i   0.0000 + 0.0000i   0.0119 + 0.1008i
   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0586i  0.0000 + 0.0000i   0.0000 + 0.0000i   0.0390 + 0.1700i 0.0000 + 0.0000i   0.0119 + 0.1008i   0.0000 + 0.0000i];

% Calulate Y Admittance matrix
for i=1:1:n
    for j=1:1:n
        if(i~=j)
            if z(i,j)==0
                y(i,j)=0; % If Z= 0 then Y = 0 as it can't be infinity
            else
                y(i,j)=-(1/z(i,j)); 
            end
        end
        
    end
   
end

x=y;

% Yii or the diagonal element of Y matrix is sum of all rows
for i=1:1:n
    temp=0;
    for j=1:1:n
            %sprintf('%i',j);
            temp= temp+x(i,j);
            y(i,i)=temp;
            y(i,i)=-y(i,i);
    end

end
alpha=zeros(n,n);
y_abs=zeros(n,n);

%COnvert the Y Bus Matrix into Polar form to get Y magnitute and alpha
%Angle
for i=1:1:n
    for j=1:1:n
        
        
        y_abs(i,j)=abs(y(i,j));   % Magnitude

        alpha(i,j)=rad2deg(angle(y(i,j)));  % In degree

        
        
    end
end



% Initialization of Values

V=[1.04,1.025,1.025,1.026,0.996,1.013,1.026,1.016,1.032]; % Power flow solution V
Theta=[0,9.3,4.7,-2.2,-4,-3.7,3.7,0.7,2]; % Corresponding del value of V 

%Machine Values
H_sec=[23.64, 6.4,3.01]; 
Xd=[0.146, 0.8958,1.3125];
Xd_dash=[0.0608,0.1198,0.1813];
Xq=[0.0969,0.8645,1.2578];
Xq_dash=[0.0969,0.1969,0.25];
Tdo_dash=[8.96,6.0,5.89];
Tqo_dash=[0.31,0.535,0.6];


Ka=[20,20,20];              %
Ta=[0.2,0.2,0.2];           %
Ke=[1,1,1];                 %
Te=[0.314,0.314,0.314];     %
Kf=[0.063,0.063,0.063];     %
Tf=[0.35,0.35,0.35];      %

D_divided_M=[0,0,0];  %pu
One_divided_M=[8,29.5,62.6];

%Initial Conditions
del=[3.58, 61.1, 54.2,0,0,0,0,0,0]; %degree
Id=[0.302, 1.29, 0.562];
Iq=[0.671, 0.931, 0.619];
Vd=[0.065,0.805, 0.779];
Vq=[1.038,0.634,0.666];
Ed_dash=[0,0.622,0.624];
Eq_dash=[1.056,0.788,0.768];
Efd=[1.082,1.789,1.403];
Rf=[0.195,0.322,0.252];
Vr=[1.105,1.902,1.453];
Vref=[1.095,1.12,1.09];
Tm=[0.716,1.63,0.85];  %Mechanical Input


% Se(Efd) Function
for i=1:1:g
    %Efd(i)
    sprintf('%i',Efd(i));
    Se_Efd(i)=0.0039*exp(1.555*Efd(i));
   
end

% fs(Efd(i)) Function
for i=1:1:g
    temp=0;
    par_diff_SeEfd=0.0039*1.555*(exp(1.555*Efd(i)));
    %par_diff_SeEfd
    temp=Ke(i)+(Efd(i)*par_diff_SeEfd)+Se_Efd(i);
    fs(i)=-temp/Te(i);
    
end



x=l*g;

A1=zeros(x,l);
B1=zeros(x,6);
B2=zeros(x,2);
E1=zeros(x,2);

%-------------------A1-------------------%
%-------------------B1-------------------%
%-------------------B2-------------------%
%-------------------E1-------------------%
for i=1:1:g
 
    % A1
     t=(i-1)*7;
     c=(i-1)*2;

     A1(t+1,t+2)=1;
     
     A1(t+2,t+2)=-D_divided_M(i);
     A1(t+2,t+3)=-Iq(i)*One_divided_M(i);
     A1(t+2,t+4)=-Id(i)*One_divided_M(i);
     
     A1(t+3,t+3)=-1/Tdo_dash(i);
     A1(t+3,t+5)=1/Tdo_dash(i);
     
     A1(t+4,t+4)=-1/Tqo_dash(i);
     
     A1(t+5,t+5)=fs(i);
     A1(t+5,t+6)=1/Te(i);
     
     A1(t+6,t+5)=-(Ka(i)*Kf(i))/(Ta(i)*Tf(i));
     A1(t+6,t+6)=-1/Ta(i);
     A1(t+6,t+7)=Ka(i)/Ta(i);

     A1(t+7,t+5)=Kf(i)/(Tf(i)*Tf(i));
     A1(t+7,t+7)=-1/Tf(i);
     
    % B1
     temp=Iq(i)*(Xd_dash(i)-Xq_dash(i));

     B1(t+2,c+1)=(temp-Ed_dash(i))*One_divided_M(i);
     
     temp=Id(i)*(Xd_dash(i)-Xq_dash(i));

     B1(t+2,c+2)=(temp-Eq_dash(i))*One_divided_M(i);
     
     B1(t+3,c+1)=-(Xd(i)-Xd_dash(i))/Tdo_dash(i);

     B1(t+4,c+2)=-(Xq(i)-Xq_dash(i))/Tqo_dash(i);


     
     % B2
     B2(t+6,c+2)=-Ka(i)/Ta(i);
     
     % E1
     E1(t+2,c+1)=One_divided_M(i);
     E1(t+6,c+2)=Ka(i)/Ta(i);
     



end


A1
B1
B2
E1

%-------------------C1-------------------%
%-------------------D1-------------------%
%-------------------D2-------------------%
C1=zeros(6,x);
D1=zeros(6,6);
D2=zeros(6,6);


for i=1:1:g

    t=(i-1)*2;
    c=(i-1)*7;

    C1(t+1,c+1)=-V(i)*cosd(del(i)-Theta(i));
    C1(t+1,c+4)=1;
    
    C1(t+2,c+1)=V(i)*sind(del(i)-Theta(i));
    C1(t+2,c+3)=1;

    D1(t+1,t+2)=Xq_dash(i);
    D1(t+2,t+1)=-Xd_dash(i);

    D2(t+1,t+1)=V(i)*cosd(del(i)-Theta(i));
    D2(t+1,t+2)=-sind(del(i)-Theta(i));
    
    D2(t+2,t+1)=-V(i)*sind(del(i)-Theta(i));
    D2(t+2,t+2)=-cosd(del(i)-Theta(i));
    
    

end



C1
D1
D2


C2=zeros(6,x);
D3=zeros(6,6);



%-------------------C2-------------------%
%-------------------D3-------------------%
for i=1:1:g
    t=(i-1)*2;
    c=(i-1)*7;
    a=del(i)-Theta(i);
    C2(t+1,c+1)=(Id(i)*V(i)*cosd(a))-(Iq(i)*V(i)*sind(a));
    C2(t+2,c+1)=-(Id(i)*V(i)*sind(a))-(Iq(i)*V(i)*cosd(a));
    
    D3(t+1,t+1)=(V(i)*sind(a));
    D3(t+1,t+2)=(V(i)*cosd(a));
    
    D3(t+2,t+1)=(V(i)*cosd(a));
    D3(t+2,t+2)=-(V(i)*sind(a));
   

end

C2
D3

%-------------------D4-------------------%

%-------------------D4(i,k) Where i=k-------------------%

i=1;
k=1;
angle2=Theta(i)-Theta(k)-alpha(i,k);
a=del(i)-Theta(i);

D411=[-(Id(i)*V(i)*cosd(a))+Iq(i)*V(i)*sind(a)   (Id(i)*sind(a))+Iq(i)*cosd(a)
    (Id(i)*V(i)*sind(a))+Iq(i)*V(i)*cosd(a)     (Id(i)*cosd(a))-(Iq(i)*sind(a))];

D411(1,2)=D411(1,2)-V(k)*y_abs(i,k)*cosd(angle2);
D411(2,2)=D411(2,2)-V(k)*y_abs(i,k)*sind(angle2);


i=2;
k=2;
angle2=Theta(i)-Theta(k)-alpha(i,k);
a=del(i)-Theta(i);

D422=[-(Id(i)*V(i)*cosd(a))+Iq(i)*V(i)*sind(a)   (Id(i)*sind(a))+Iq(i)*cosd(a)
    (Id(i)*V(i)*sind(a))+Iq(i)*V(i)*cosd(a)     (Id(i)*cosd(a))-(Iq(i)*sind(a))];

D422(1,2)=D422(1,2)+(-V(k)*y_abs(i,k)*cosd(angle2));
D422(2,2)=D422(2,2)+(-V(k)*y_abs(i,k)*sind(angle2));


i=3;
k=3;
angle2=Theta(i)-Theta(k)-alpha(i,k);
a=del(i)-Theta(i);


D433=[-(Id(i)*V(i)*cosd(a))+Iq(i)*V(i)*sind(a)   (Id(i)*sind(a))+Iq(i)*cosd(a)
    (Id(i)*V(i)*sind(a))+Iq(i)*V(i)*cosd(a)     (Id(i)*cosd(a))-(Iq(i)*sind(a))];

D433(1,2)=D433(1,2)+(-V(k)*y_abs(i,k)*cosd(angle2));
D433(2,2)=D433(2,2)+(-V(k)*y_abs(i,k)*sind(angle2));

sum_V_real=0;
sum_V_react=0;
sum_Theta_real=0;
sum_Theta_react=0;

for k=1:1:3
    i=1;
    angle2=Theta(i)-Theta(k)-alpha(i,k);
    sum_V_real=sum_V_real-(V(k)*y_abs(i,k)*cosd(angle2));
    sum_V_react=sum_V_react-(V(k)*y_abs(i,k)*sind(angle2));

    if(i~=k)
        
        angle2=Theta(i)-Theta(k)-alpha(i,k);
        sum_Theta_real=sum_Theta_real-(V(i)*V(k)*y_abs(i,k)*sind(angle2));
        sum_Theta_react=sum_Theta_react+(V(i)*V(k)*y_abs(i,k)*cosd(angle2));
        
    end

end

D411(1,2)=D411(1,2)+sum_V_real;
D411(2,2)=D411(2,2)+sum_V_react;

D411(1,1)=D411(1,1)+sum_Theta_real;
D411(2,1)=D411(2,1)+sum_Theta_react;

sum_V_real=0;
sum_V_react=0;
sum_Theta_real=0;
sum_Theta_react=0;

for k=1:1:3
    i=2;
    angle2=Theta(i)-Theta(k)-alpha(i,k);
    sum_V_real=sum_V_real+(-V(k)*y_abs(i,k)*cosd(angle2));
    sum_V_react=sum_V_react+(-V(k)*y_abs(i,k)*sind(angle2));

    if(i~=k)
        angle2=Theta(i)-Theta(k)-alpha(i,k);
        sum_Theta_real=sum_Theta_real-(V(i)*V(k)*y_abs(i,k)*sind(angle2));
        sum_Theta_react=sum_Theta_react+(V(i)*V(k)*y_abs(i,k)*cosd(angle2));
        
    end

end

D422(1,2)=D422(1,2)+sum_V_real;
D422(2,2)=D422(2,2)+sum_V_react;

D422(1,1)=D422(1,1)+sum_Theta_real;
D422(2,1)=D422(2,1)+sum_Theta_react;

D422;

sum_V_real=0;
sum_V_react=0;
sum_Theta_real=0;
sum_Theta_react=0;

for k=1:1:3
    i=3;
    angle2=Theta(i)-Theta(k)-alpha(i,k);
    sum_V_real=sum_V_real+(-V(k)*y_abs(i,k)*cosd(angle2));
    sum_V_react=sum_V_react+(-V(k)*y_abs(i,k)*sind(angle2));

    if(i~=k)
        angle2=Theta(i)-Theta(k)-alpha(i,k);
        sum_Theta_real=sum_Theta_real-(V(i)*V(k)*y_abs(i,k)*sind(angle2));
        sum_Theta_react=sum_Theta_react+(V(i)*V(k)*y_abs(i,k)*cosd(angle2));
        
    end

end
D433(1,2)=D433(1,2)+sum_V_real;
D433(2,2)=D433(2,2)+sum_V_react;

D433(1,1)=D433(1,1)+sum_Theta_real;
D433(2,1)=D433(2,1)+sum_Theta_react;

D433;


%-------------------D4(i,k) Where i!=k-------------------%
%---------------Calulating only the values where i node is connected to k-------------------%

i=1;
k=2;
angle2=Theta(i)-Theta(k)-alpha(i,k);

D412 = [-V(i)*V(k)*y_abs(i,k)*sind(angle2)  -V(i)*y_abs(i,k)*cosd(angle2)
       V(i)*V(k)*y_abs(i,k)*cosd(angle2)   -V(i)*y_abs(i,k)*sind(angle2)
];

i=1;
k=3;
angle2=Theta(i)-Theta(k)-alpha(i,k);

D413 = [-V(i)*V(k)*y_abs(i,k)*sind(angle2)  -V(i)*y_abs(i,k)*cosd(angle2)
       V(i)*V(k)*y_abs(i,k)*cosd(angle2)   -V(i)*y_abs(i,k)*sind(angle2)
];

i=2;
k=1;
angle2=Theta(i)-Theta(k)-alpha(i,k);

D421 = [-V(i)*V(k)*y_abs(i,k)*sind(angle2)  -V(i)*y_abs(i,k)*cosd(angle2)
       V(i)*V(k)*y_abs(i,k)*cosd(angle2)   -V(i)*y_abs(i,k)*sind(angle2)
];

i=2;
k=3;
angle2=Theta(i)-Theta(k)-alpha(i,k);

D423 = [-V(i)*V(k)*y_abs(i,k)*sind(angle2)  -V(i)*y_abs(i,k)*cosd(angle2)
       V(i)*V(k)*y_abs(i,k)*cosd(angle2)   -V(i)*y_abs(i,k)*sind(angle2)
];


i=3;
k=1;
angle2=Theta(i)-Theta(k)-alpha(i,k);

D431 = [-V(i)*V(k)*y_abs(i,k)*sind(angle2)  -V(i)*y_abs(i,k)*cosd(angle2)
       V(i)*V(k)*y_abs(i,k)*cosd(angle2)   -V(i)*y_abs(i,k)*sind(angle2)
];

i=3;
k=2;
angle2=Theta(i)-Theta(k)-alpha(i,k);

D432 = [-V(i)*V(k)*y_abs(i,k)*sind(angle2)  -V(i)*y_abs(i,k)*cosd(angle2)
       V(i)*V(k)*y_abs(i,k)*cosd(angle2)   -V(i)*y_abs(i,k)*sind(angle2)
];

D432;

D4=[D411 D412 D413
    D421 D422 D423
    D431 D432 D433];

D4

%---------------------D5------------------------------%
%---------------Calulating only the values where i node is connected to k-------------------%


i=1;
k=4;
angle2=Theta(i)-Theta(k)-alpha(i,k);

D514 = [-V(i)*V(k)*y_abs(i,k)*sind(angle2)  -V(i)*y_abs(i,k)*cosd(angle2)
       V(i)*V(k)*y_abs(i,k)*cosd(angle2)   -V(i)*y_abs(i,k)*sind(angle2)
];

D514;

i=2;
k=7;
angle2=Theta(i)-Theta(k)-alpha(i,k);

D527 = [-V(i)*V(k)*y_abs(i,k)*sind(angle2)  -V(i)*y_abs(i,k)*cosd(angle2)
       V(i)*V(k)*y_abs(i,k)*cosd(angle2)   -V(i)*y_abs(i,k)*sind(angle2)
];


D527;


i=3;
k=9;
angle2=Theta(i)-Theta(k)-alpha(i,k);

D539 = [-V(i)*V(k)*y_abs(i,k)*sind(angle2)  -V(i)*y_abs(i,k)*cosd(angle2)
       V(i)*V(k)*y_abs(i,k)*cosd(angle2)   -V(i)*y_abs(i,k)*sind(angle2)
];


D539;

empty_D=zeros(2,2);
D5=[D514 empty_D empty_D empty_D empty_D empty_D
    empty_D empty_D empty_D D527 empty_D empty_D
    empty_D empty_D empty_D empty_D empty_D D539];

D5

%---------------------D6------------------------------%
%---------------Calulating only the values where i node is connected to k-------------------%


i=4;
k=1;
angle2=Theta(i)-Theta(k)-alpha(i,k);

D641 = [-V(i)*V(k)*y_abs(i,k)*sind(angle2)  -V(i)*y_abs(i,k)*cosd(angle2)
       V(i)*V(k)*y_abs(i,k)*cosd(angle2)   -V(i)*y_abs(i,k)*sind(angle2)
];

D641;

i=7;
k=2;
angle2=Theta(i)-Theta(k)-alpha(i,k);

D672 =[-V(i)*V(k)*y_abs(i,k)*sind(angle2)  -V(i)*y_abs(i,k)*cosd(angle2)
       V(i)*V(k)*y_abs(i,k)*cosd(angle2)   -V(i)*y_abs(i,k)*sind(angle2)
];


D672;

i=9;
k=3;
angle2=Theta(i)-Theta(k)-alpha(i,k);

D693 = [-V(i)*V(k)*y_abs(i,k)*sind(angle2)  -V(i)*y_abs(i,k)*cosd(angle2)
       V(i)*V(k)*y_abs(i,k)*cosd(angle2)   -V(i)*y_abs(i,k)*sind(angle2)
];


D693;
empty_D=zeros(2,2);

D6=[D641 empty_D empty_D
    empty_D empty_D empty_D
    empty_D empty_D empty_D
    empty_D D672 empty_D
    empty_D empty_D empty_D
    empty_D empty_D D693];

D6

%---------------------D7------------------------------%
sum_Theta_real=0;
sum_Theta_react=0;
sum_V_real=0;
sum_V_react=0;

for k=1:1:9
    i=4;
    angle2=Theta(i)-Theta(k)-alpha(i,k);
    sum_V_real=sum_V_real+(V(k)*y_abs(i,k)*cosd(angle2));
    sum_V_react=sum_V_react+(V(k)*y_abs(i,k)*sind(angle2));

    if(i~=k)
        angle2=Theta(i)-Theta(k)-alpha(i,k);
        sum_Theta_real=sum_Theta_real+(V(i)*V(k)*y_abs(i,k)*sind(angle2));
        sum_Theta_react=sum_Theta_react+(V(i)*V(k)*y_abs(i,k)*cosd(angle2));
        
    end

end
i=4;
k=4;
angle2=Theta(i)-Theta(k)-alpha(i,k);

D744(1,1)=sum_Theta_real;
D744(2,1)=-sum_Theta_react;
D744(1,2)=-sum_V_real+(-V(k)*y_abs(i,k)*cosd(angle2));
D744(2,2)=-sum_V_react+((-V(k)*y_abs(i,k)*sind(angle2)));

D744;

i=5;
k=5;
sum_Theta_real=0;
sum_Theta_react=0;
sum_V_real=0;
sum_V_react=0;

for k=1:1:9
    i=5;
    angle2=Theta(i)-Theta(k)-alpha(i,k);
    sum_V_real=sum_V_real+(V(k)*y_abs(i,k)*cosd(angle2));
    sum_V_react=sum_V_react+(V(k)*y_abs(i,k)*sind(angle2));

    if(i~=k)
        angle2=Theta(i)-Theta(k)-alpha(i,k);
        sum_Theta_real=sum_Theta_real+(V(i)*V(k)*y_abs(i,k)*sind(angle2));
        sum_Theta_react=sum_Theta_react+(V(i)*V(k)*y_abs(i,k)*cosd(angle2));
        
    end

end

i=5;
k=5;
angle2=Theta(i)-Theta(k)-alpha(i,k);

D755(1,1)=sum_Theta_real;
D755(2,1)=-sum_Theta_react;
D755(1,2)=-sum_V_real+(-V(k)*y_abs(i,k)*cosd(angle2));
D755(2,2)=-sum_V_react+((-V(k)*y_abs(i,k)*sind(angle2)));

D755;


sum_Theta_real=0;
sum_Theta_react=0;
sum_V_real=0;
sum_V_react=0;

for k=1:1:9
    i=6;
   angle2=Theta(i)-Theta(k)-alpha(i,k);
    sum_V_real=sum_V_real+(V(k)*y_abs(i,k)*cosd(angle2));
    sum_V_react=sum_V_react+(V(k)*y_abs(i,k)*sind(angle2));

    if(i~=k)
        angle2=Theta(i)-Theta(k)-alpha(i,k);
        sum_Theta_real=sum_Theta_real+(V(i)*V(k)*y_abs(i,k)*sind(angle2));
        sum_Theta_react=sum_Theta_react+(V(i)*V(k)*y_abs(i,k)*cosd(angle2));
        
    end

end

i=6;
k=6;
angle2=Theta(i)-Theta(k)-alpha(i,k);

D766(1,1)=sum_Theta_real;
D766(2,1)=-sum_Theta_react;
D766(1,2)=-sum_V_real+(-V(k)*y_abs(i,k)*cosd(angle2));
D766(2,2)=-sum_V_react+((-V(k)*y_abs(i,k)*sind(angle2)));

D766;


sum_Theta_real=0;
sum_Theta_react=0;
sum_V_real=0;
sum_V_react=0;

for k=1:1:9
    i=7;
    angle2=Theta(i)-Theta(k)-alpha(i,k);
    sum_V_real=sum_V_real+(V(k)*y_abs(i,k)*cosd(angle2));
    sum_V_react=sum_V_react+(V(k)*y_abs(i,k)*sind(angle2));

    if(i~=k)
        angle2=Theta(i)-Theta(k)-alpha(i,k);
        sum_Theta_real=sum_Theta_real+(V(i)*V(k)*y_abs(i,k)*sind(angle2));
        sum_Theta_react=sum_Theta_react+(V(i)*V(k)*y_abs(i,k)*cosd(angle2));
        
    end

end

i=7;
k=7;
angle2=Theta(i)-Theta(k)-alpha(i,k);

D777(1,1)=sum_Theta_real;
D777(2,1)=-sum_Theta_react;
D777(1,2)=-sum_V_real+(-V(k)*y_abs(i,k)*cosd(angle2));
D777(2,2)=-sum_V_react+((-V(k)*y_abs(i,k)*sind(angle2)));

D777;


sum_Theta_real=0;
sum_Theta_react=0;
sum_V_real=0;
sum_V_react=0;

for k=1:1:9
    i=8;
    angle2=Theta(i)-Theta(k)-alpha(i,k);
    sum_V_real=sum_V_real+(V(k)*y_abs(i,k)*cosd(angle2));
    sum_V_react=sum_V_react+(V(k)*y_abs(i,k)*sind(angle2));

    if(i~=k)
        angle2=Theta(i)-Theta(k)-alpha(i,k);
        sum_Theta_real=sum_Theta_real+(V(i)*V(k)*y_abs(i,k)*sind(angle2));
        sum_Theta_react=sum_Theta_react+(V(i)*V(k)*y_abs(i,k)*cosd(angle2));
        
    end

end

i=8;
k=8;
angle2=Theta(i)-Theta(k)-alpha(i,k);

D788(1,1)=sum_Theta_real;
D788(2,1)=-sum_Theta_react;
D788(1,2)=-sum_V_real+(-V(k)*y_abs(i,k)*cosd(angle2));
D788(2,2)=-sum_V_react+((-V(k)*y_abs(i,k)*sind(angle2)));

D788;


sum_Theta_real=0;
sum_Theta_react=0;
sum_V_real=0;
sum_V_react=0;

for k=1:1:9
    i=9;
   angle2=Theta(i)-Theta(k)-alpha(i,k);
    sum_V_real=sum_V_real+(V(k)*y_abs(i,k)*cosd(angle2));
    sum_V_react=sum_V_react+(V(k)*y_abs(i,k)*sind(angle2));

    if(i~=k)
        angle2=Theta(i)-Theta(k)-alpha(i,k);
        sum_Theta_real=sum_Theta_real+(V(i)*V(k)*y_abs(i,k)*sind(angle2));
        sum_Theta_react=sum_Theta_react+(V(i)*V(k)*y_abs(i,k)*cosd(angle2));
        
    end

end

i=9;
k=9;
angle2=Theta(i)-Theta(k)-alpha(i,k);

D799(1,1)=sum_Theta_real;
D799(2,1)=-sum_Theta_react;
% D799(1,2)=-sum_V_real;
% D799(2,2)=-sum_V_react;

D799(1,2)=-sum_V_real+(-V(k)*y_abs(i,k)*cosd(angle2));
D799(2,2)=-sum_V_react+((-V(k)*y_abs(i,k)*sind(angle2)));

D799;

i=4;
k=5;
angle2=Theta(i)-Theta(k)-alpha(i,k);

D745 = [-V(i)*V(k)*y_abs(i,k)*sind(angle2)  -V(k)*y_abs(i,k)*cosd(angle2)
       V(i)*V(k)*y_abs(i,k)*cosd(angle2)   -V(k)*y_abs(i,k)*sind(angle2)
];

D745;


i=4;
k=6;
angle2=Theta(i)-Theta(k)-alpha(i,k);

D746 = [-V(i)*V(k)*y_abs(i,k)*sind(angle2)  -V(i)*y_abs(i,k)*cosd(angle2)
       V(i)*V(k)*y_abs(i,k)*cosd(angle2)   -V(i)*y_abs(i,k)*sind(angle2)
];


D746;


i=5;
k=4;
angle2=Theta(i)-Theta(k)-alpha(i,k);

D754 = [-V(i)*V(k)*y_abs(i,k)*sind(angle2)  -V(i)*y_abs(i,k)*cosd(angle2)
       V(i)*V(k)*y_abs(i,k)*cosd(angle2)   -V(i)*y_abs(i,k)*sind(angle2)
];


D754;

i=5;
k=7;
angle2=Theta(i)-Theta(k)-alpha(i,k);

D757 = [-V(i)*V(k)*y_abs(i,k)*sind(angle2)  -V(i)*y_abs(i,k)*cosd(angle2)
       V(i)*V(k)*y_abs(i,k)*cosd(angle2)   -V(i)*y_abs(i,k)*sind(angle2)
];

D757;

i=6;
k=4;
angle2=Theta(i)-Theta(k)-alpha(i,k);

D764 = [-V(i)*V(k)*y_abs(i,k)*sind(angle2)  -V(i)*y_abs(i,k)*cosd(angle2)
       V(i)*V(k)*y_abs(i,k)*cosd(angle2)   -V(i)*y_abs(i,k)*sind(angle2)
];


D764;

i=6;
k=9;
angle2=Theta(i)-Theta(k)-alpha(i,k);

D769 = [-V(i)*V(k)*y_abs(i,k)*sind(angle2)  -V(i)*y_abs(i,k)*cosd(angle2)
       V(i)*V(k)*y_abs(i,k)*cosd(angle2)   -V(i)*y_abs(i,k)*sind(angle2)
];

D769;

i=7;
k=5;
angle2=Theta(i)-Theta(k)-alpha(i,k);

D775 = [-V(i)*V(k)*y_abs(i,k)*sind(angle2)  -V(i)*y_abs(i,k)*cosd(angle2)
       V(i)*V(k)*y_abs(i,k)*cosd(angle2)   -V(i)*y_abs(i,k)*sind(angle2)
];


D775;

i=7;
k=8;
angle2=Theta(i)-Theta(k)-alpha(i,k);

D778 = [-V(i)*V(k)*y_abs(i,k)*sind(angle2)  -V(i)*y_abs(i,k)*cosd(angle2)
       V(i)*V(k)*y_abs(i,k)*cosd(angle2)   -V(i)*y_abs(i,k)*sind(angle2)
];

D778;


i=8;
k=7;
angle2=Theta(i)-Theta(k)-alpha(i,k);

D787 = [-V(i)*V(k)*y_abs(i,k)*sind(angle2)  -V(i)*y_abs(i,k)*cosd(angle2)
       V(i)*V(k)*y_abs(i,k)*cosd(angle2)   -V(i)*y_abs(i,k)*sind(angle2)
];


D787;

i=8;
k=9;
angle2=Theta(i)-Theta(k)-alpha(i,k);

D789 = [-V(i)*V(k)*y_abs(i,k)*sind(angle2)  -V(i)*y_abs(i,k)*cosd(angle2)
       V(i)*V(k)*y_abs(i,k)*cosd(angle2)   -V(i)*y_abs(i,k)*sind(angle2)
];


D789;

i=9;
k=6;
angle2=Theta(i)-Theta(k)-alpha(i,k);

D796 = [-V(i)*V(k)*y_abs(i,k)*sind(angle2)  -V(i)*y_abs(i,k)*cosd(angle2)
       V(i)*V(k)*y_abs(i,k)*cosd(angle2)   -V(i)*y_abs(i,k)*sind(angle2)
];

D796;




i=9;
k=8;
angle2=Theta(i)-Theta(k)-alpha(i,k);

D798 = [-V(i)*V(k)*y_abs(i,k)*sind(angle2)  -V(i)*y_abs(i,k)*cosd(angle2)
       V(i)*V(k)*y_abs(i,k)*cosd(angle2)   -V(i)*y_abs(i,k)*sind(angle2)
];


D798;

D7=[D744 D745 D746 empty_D empty_D empty_D
    D754 D755 empty_D D757 empty_D empty_D
    D764 empty_D D766 empty_D empty_D D769
    empty_D D775 empty_D D777 D778 empty_D
    empty_D empty_D empty_D D787 D788 D789
    empty_D empty_D D796 empty_D D798 D799];

D7

%---------------Solving Equation for ∆x in 8.51, 8.50, 8.49-------------------%

K1=D4-(D3*inv(D1)*D2);
K2=C2-(D3*inv(D1)*C1);

K1;
K2;
K3=inv(K1)*(-K2);

K4=D5*inv(D7)*D6;

K5=inv(-K4)*(-K2);
K5;

K6=K3+K5;

K7=(B2-(B1*inv(D1)*D2))*K6;
K7;

K8=A1-B1*inv(D1)*C1;

Asym=K8+K7;
Asym;

% 
% Asym(1,:)=[];  % Removing node w.r.t to generator 1 (Slack Bus) and making it del reference
% Asym(:,1)=[];
% %  
% Asym(1,1)=-1;  % Changing the del reference w.r.t to Generator 1 as 
% Asym(2,1)=-1;
% %  
 e=eig(Asym);

fprintf('--------------------Natural Mode of the System-------------------\n\n');

e
