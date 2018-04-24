%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lecture: Digital Communication
% Gram-Schmidt procedure is used to obtain set of BASIS FUNCTIONS from a
% given set of signal Elements

% Author: Abdur Rahman Mohamed Ismail
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;
close all;
M=4;   %No of signals
%%%% The Given signals (Checking for DiCo exam WS 2018)

s1=[1 1 0];
s2=[1 -1 0]; 
s3=[1 1 -1];
s4=[0 0 -1];
%%%%%%%%%%%%%%%%%%%%


Smat=[s1;s2;s3;s4]'; %signal_matrix

Umat=zeros(length(s1),M); % U vector matrix
Umatc=zeros(1,M);         % Normalising constants of U vector matrix
Umat(:,1)=Smat(:,1);      % First elements are always same
Umatc(1,1)=1/(sqrt(sum(Smat(:,1).*Smat(:,1)))); 
temp=zeros(length(s1),1);
temp1=zeros(length(s1),1);


for i=2:M
     % Finding the span of vectors (projection)    
    for j=1:i-1
           temp = temp + ( (Umatc(1,j)^2)*sum(Smat(:,i).*Umat(:,j)) ).*Umat(:,j);
    end

    temp1=(Smat(:,i))-temp(:,1); 
    Umat(:,i)= temp1;
   
    Umatc(1,i)= 1/(sqrt(sum(round(Umat(:,i).^2,4)))); % Normalising formula
    temp=0;
end
 
disp("The basis function of given signals are given as a matrix:");
disp("Note:Inf and NaN are invalid basis function and can be ignored");
 Basis_fn=Umat.*Umatc

 

    