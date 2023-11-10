clear all;

fprintf('--------------------Admittance Matrix------------------------------\n\n');


i=0;
j=0;
z=[]; % Impedence Matix
y=[]; % Admitance Matrix
H=[]; %∂ P/∂δ 
K=[]; %∂ P/∂V
M=[]; %∂ Q/∂V
K=[]; %∂ P/∂δ
Jac=[]; % Jacobian Matrix
n=9; % Length of bus



z= [0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i 0.0000 + 0.0576i   0.0000 + 0.0000i   0.0000 + 0.0000i 0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i 0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i  0.0000 + 0.0625i   0.0000 + 0.0000i   0.0000 + 0.0000i
   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i 0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i 0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0586i
   0.0000 + 0.0576i   0.0000 + 0.0000i   0.0000 + 0.0000i 0.0000 + 0.0000i   0.0100 + 0.0850i   0.0170 + 0.0920i 0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i
   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i 0.0100 + 0.0850i   0.0000 + 0.0000i   0.0000 + 0.0000i 0.0320 + 0.1610i   0.0000 + 0.0000i   0.0000 + 0.0000i
   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i 0.0170 + 0.0920i   0.0000 + 0.0000i   0.0000 + 0.0000i 0.0000 + 0.0000i   0.0000 + 0.0000i   0.0390 + 0.1700i
   0.0000 + 0.0000i   0.0000 + 0.0625i   0.0000 + 0.0000i 0.0000 + 0.0000i   0.0320 + 0.1610i   0.0000 + 0.0000i 0.0000 + 0.0000i   0.0085 + 0.0720i   0.0000 + 0.0000i
   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i 0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0000i 0.0085 + 0.0720i   0.0000 + 0.0000i   0.0119 + 0.1008i
   0.0000 + 0.0000i   0.0000 + 0.0000i   0.0000 + 0.0586i  0.0000 + 0.0000i   0.0000 + 0.0000i   0.0390 + 0.1700i 0.0000 + 0.0000i   0.0119 + 0.1008i   0.0000 + 0.0000i];

z

%Uncomment the following code if you want to enter Z manually
%n=input("Enter number of nodes:");
% for i=1:1:n
%     for j=1:1:n
% 
%         if(i<=j) % only taking the first half of matrix as we it is
%         symmetrica matrix
%             z(i,j,1)=input(sprintf("Enter Z%i%i :", i ,j));
%         else
%             z(i,j)=z(j,i);
%         end
%         
%     end
%    
% end
%z

% Calculate Y Admittance matrix

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


% Print Admittance Matrix
y

fprintf('--------------------Jacobian Matrix-------------------\n\n');

V=[1.04,1.025,1.025,1.026,0.996,1.013,1.026,1.016,1.032]; % Power flow solution V
del=[0,9.3,4.7,-2.2,-4,-3.7,3.7,0.7,2]; % Corresponding del value of V 


fprintf('--------------------H Matrix-------------------\n\n');

for i=1:1:n
    for j=1:1:n
           if(i~=j)
               B=imag(y(i,j)); % Getting imaginary value of yij
               G=real(y(i,j)); % Getting real value of yij
               a=del(i)-del(j);
               H(i,j)=-V(i)*V(j)*(B*cosh(a)-G*sinh(a));
           end

        
    end

end
x=H;
% Finding Diagonal value of Hii matrix which is sum of all element in the
% row
for i=1:1:n
    temp=0;
    for j=1:1:n
           temp= temp+x(i,j);
            H(i,i)=temp;
           
    end

end
H

fprintf('--------------------N Matrix-------------------\n\n');

for i=1:1:n
    for j=1:1:n
           if(i~=j)
               B=imag(y(i,j));
               G=real(y(i,j));
               a=del(i)-del(j);
               N(i,j)=-(V(i)*V(j)*(G*cosh(a)+B*sinh(a)));
           end

        
    end

end
x=N;

%Digonal Value
for i=1:1:n
    temp=0;
    for j=1:1:n
           temp= temp+x(i,j);
            N(i,i)=temp;
           
    end

end
N

fprintf('--------------------M Matrix-------------------\n\n');

for i=1:1:n
    for j=1:1:n
           if(i~=j)
               B=imag(y(i,j));
               G=real(y(i,j));
               a=del(i)-del(j);
               M(i,j)=(V(i)*V(j)*(G*cosh(a)+B*sinh(a)));
           end

        
    end

end
x=M;
for i=1:1:n
    temp=0;
    for j=1:1:n
           temp= temp+x(i,j);
            M(i,i)=temp;
           
    end
     G=real(y(i,i));
     %G
    M(i,i)=M(i,i)+(2*V(i)*V(i)*G);

end
M

fprintf('--------------------K Matrix-------------------\n\n');

for i=1:1:n
    for j=1:1:n
           if(i~=j)
               B=imag(y(i,j));
               G=real(y(i,j));
               a=del(i)-del(j);
               K(i,j)=V(i)*V(j)*(G*sinh(a)-B*cosh(a));
           end

        
    end

end
x=K;
x
for i=1:1:n
    temp=0;
    for j=1:1:n
           temp= temp+x(i,j);
            K(i,i)=temp;
          
    end
     B=imag(y(i,i));
     %B
    K(i,i)=K(i,i)-(2*V(i)*V(i)*B);
end
K
fprintf('--------------------Final Jacobian Matrix-------------------\n\n');

Jac= [H M; N K]; % concating H, N M and K
Jac

% Eigen value of Jacobian Matrix
E=eig(Jac);
E

