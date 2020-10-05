
clear all
close all
clc

% This syetm that you need to solve will be singular. Matlab gives you a
% warning at each time step. To switch of this warning, remove the comment
% in the next line

%warning on


Re = 1000;              % Reynolds number
N =64;                  % Number of volumes in the x- and y-direction
Delta = 1/N;            % uniform spacing to be used in the mapping to compute tx


% Determine a suitable time step and stopping criterion, tol
%dt =2.8e-4;             % time step
tol =10^(-6);             % tol determines when steady state is reached and the program terminates

% wall velocities
U_wall_top = -1;
U_wall_bot = 0;
U_wall_left = 0;
U_wall_right = 0;
V_wall_top = 0;
V_wall_bot = 0;
V_wall_left = 0;
V_wall_right = 0;

%
%   Generation of a non-uniform mesh
%

%
%   tx are the coordinates of the nodal points on the outer-oriented grid
%
tx = zeros(1,N+1);
for i=1:N+1
    xi = (i-1)*Delta;
    tx(i) = 0.5*(1. - cos(pi*xi));
end

% Local mesh size on outer oriented grid
th = zeros(N,1);
th = tx(2:N+1) - tx(1:N);

%
%  x are the coordinates of the nodal points on the inner-oriented grid (including
%  endpoints 0 and 1)
%  h contains the edge lengths on the inner-oriented grid
%
x = 0.5*(tx(1:N) + tx(2:N+1));
x = [0 x 1];

h = zeros(N+1,1);
h = x(2:N+2) - x(1:N+1);

%
%   Initial condition u=v=0
%
%   Both u and v will be stored in one big vector called 'u'
%
%   The vector u only contains the true unknowns, not the velocities that
%   are prescribed by the boundary conditions
%
%   The vector u contains the *inner-oriented* circulations as unknowns

u = zeros(2*N*(N-1),1);

%% Creating tE21
 %Preparing diagonals for tE21
B=zeros(N*N,N+3); %B is a matrix containing diagonals of sparse matrix tE21
B(:,N+2)=(-1)*ones(N*N,1); %Diagonal containing incoming fluxes in y-direction
B(:,N+3)=ones(N*N,1); %Diagonal containing outgoing fluxes in y-direction
B(1:N,1)=-1;
B(N*N-N+1:N*N,N+1)=1;

%Diagonal containing indexes related to fluxes in x-direction
temp1=(-1)*ones(2*N,1);
temp1(1:N)=temp1(1:N)*(-1);

for i=2:N
    for j=1:2*N
        B(j+(i-2)*N,i)=temp1(j);
    end
end

%d contains the placement of the diagonals in tE21 contained as collums in B

d=zeros(N+3,1);
for i=2:N+1
    d(i)=i-1;
end
d(N+2)=(N+1)*N;
d(N+3)=(N+2)*N;

tE21=spdiags(B,d,N*N,2*(N+1)*N);


%Removing colums of tE21 corresponding to predescribed fluxes and make
%u_normal
B=sparse(N*N,4*N); % Matrix to save columns removed from tE21
%Removinng columns related to flux on the upper boundary (y=1)
j=4*N; 
for i=2*N*(N+1):-1:2*N*(N+1)-N+1
    B(:,j)=tE21(:,i);
    tE21(:,i)=[];
    j=j-1;
end
%Removinng columns related to flux on the lower boundary (y=0)
for i=N*(N+2):-1:N*(N+1)+1
    B(:,j)=tE21(:,i);
    tE21(:,i)=[];
    j=j-1;
end

%Remowing coumns related to flux on the left boundary (x=0)
j=1;
for i=1:N
    B(:,j)=tE21(:,1+(i-1)*N);
    tE21(:,1+(i-1)*N)=[];
    j=j+2;
end
%Remowing columns relted to flux on the right boundary (x=1)
j=2;
for i=N:N-1:N*(N-1)+1
    B(:,j)=tE21(:,i);
    tE21(:,i)=[];
    j=j+2;
end
%disp(full(tE21));
%disp(full(B));

%Storing Wall Velocities in vector u_h and relate wall normal fluxes to the
%boundary condition conditions for the wall velocities. According to the
%definition of the wall fluxes, these should all be zero

u_norm=zeros(N*N,1);


%Find E10 and show its relation to tE21
E10=tE21';

%% Setting up E21

%Prepering diagonals for E21 
B=zeros((N+1)*(N+1),N+4); %B is a matrix containing diagonals of sparse matrix E21

d=zeros(N+4,1); %d is the location of the diagonals
d(2)=N+1;
d(N+4)=2*(N+2)*(N+1);
iter=3;
for i=(N+2)*(N+1):(N+2)*(N+1)+(N+1)
    d(iter)=i;
    iter=iter+1;
end

B(:,1)=1; %Circulations in positive x-direction in circulation direction
B(:,2)=-1; %Circulations in negative x-direction in circulation direction
B(1:N+1,3)=-1; %Circulations perpendicular to y=0 (coming out from boundary)
B((N+1)*N+1:(N+1)*(N+1),N+4)=1; %Circulations perpendicular to y=1 (going into boundary)
%Prepare other diagonals containing circulations in y direction
temp1=ones(2*(N+1),1);
temp1(N+2:(2*(N+1)),1)=-1;

iter=0;
for j=4:N+3
    for i=1:2*(N+1)
        B(iter+i,j)=temp1(i);
    end
   iter=iter+(N+1);
end

E21=spdiags(B,d,(N+1)*(N+1),2*(N+2)*(N+1));
%Remove columns that are prescribed


%removing prescribed along circulations at y=1 boundary
for i=(N+2)*(N+1):-1:(N+1)*(N+1)+1
    E21(:,i)=[];
end
%removing prescribed along circulations at y=0 boundary
for i=(N+1):-1:1
    E21(:,i)=[];
end
%removing prescribed circulations along  x=0 boundary
for i=0:N
    E21(:,(N+1)*N+1+(N+1)*i)=[];
end
%removing prescribed circulations along  x=1 boundary
for i=0:N
    E21(:,(N+1)*N+(N+1)+N*i)=[];
end
%removing prescribed circulations normal to x=0 
E21(:,1)=[];
for i=1:N-1
    E21(:,i*N+1)=[];
end
%removing prescribed circulations normal to x=1
for i=1:N
    E21(:,i*(N-1)+1)=[];
end
%removing prescribed circulations normal to y=0
for i=1:N
    E21(:,N*(N-1)+1)=[];
end
%removing prescribed circulations normal to y=0
temp=size(E21);
temp=temp(2);
for i=1:N
    E21(:,temp)=[];
    temp=temp-1;
end

%% Finding u_pres
%u_pres
u_pres=zeros((N+1)*(N+1),1);

temp=((N+1)*(N+1)-N);
counter=1;
for i=temp:(N+1)*(N+1)
    u_pres(i)=U_wall_top*h(counter)*(-1);
    counter=counter+1;
end

%% Set up Ht02

areaVec=zeros((N+1)^2,1);
for i=0:N
    for j=1:N+1
        areaVec(i*(N+1)+j)=h(i+1)*h(j);
    end
    
end
diag=zeros(((N+1)^2),1);

for i=0:N
    for j=1:N+1
        diag(i*(N+1)+j)=1/areaVec(i*(N+1)+j);
        
    end
end

Ht02=spdiags(diag,0,(N+1)*(N+1),(N+1)*(N+1));

%% Set up H1t1



H1t1 = sparse(N*N+N,N*N+N);

counter = 0;

for j = 1 : size(th,2)
    for g = 2 : size(h,2)-1
        counter = counter + 1;
        H1t1(counter,counter) = h(g)/th(j);
    end 
end
for g = 2 : size(h,2)-1
    for j = 1 : size(th,2)
        counter = counter +1;
        H1t1(counter,counter) = h(g)/th(j);
    end
end

% Hu_norm is the vector which will contain the Hodge of the prescribed
% normal fluxes. Calculate the vector 'Hu_norm'

%Since the vector u_norm is a zero vector, the vector operation done by
%applying the hodge to it is also a zero vector

%% Finding Ht11
%This Matrix should be the inverse of H1t1
Ht11=inv(H1t1);
%% Time steps
hmin=min(h);

%dt=min(hmin,((1/2)*Re*hmin^2));
%The time steps used differ a bit. When N gets high, it becomes more
%convinient to follow a timestep a little bit higher than what indicated by
%the smallest cell size.
dt=2.5e-4;


%% Solver code
Integrated_vor_pre=u_pres; %Save prescribed intergrated vorticity
u_pres_tvort=Ht02*u_pres; %U_pres to outer oriented 0 form representing contribution of boundary conditions to point wise vorticity
u_pres = H1t1*E21'*Ht02*u_pres; %U_pres to inner oriented 1 forms

diff = 1;
iter = 1;


% Set up the matrix for the Poisson equation    

A = -tE21*Ht11*tE21';


% Perform an LU-decomposition for the pressure matrix A

[L,U] = lu(A);

% Abbreviation for some matrix products which are constant in the time loop

VLaplace = H1t1*E21'*Ht02*E21;
DIV = tE21*Ht11;

t=0; %See how long time it takes before programme thinks steady state is reached

while diff>tol
        
    %Vector chi is obtained. It corresponds with the point-wise vorticity
    %at each cell
    chi=Ht02*E21*u+u_pres_tvort;  
    
    %Vectors uxchi and uychi correspond with the multiplications of
    %chi with the horizontal and vertical velocity components at each cell.
    %Only the cells required for vector convective are calculated. The
    %ordering of each vector with respect to the ordering of cells in the
    %grid is different (left to right for uxchi and bottom to top for
    %uychi)
    
    uxchi=zeros((N+1)*(N-1),1);
    uychi=zeros((N+1)*(N-1),1);
    
    for i=1:N-1
        for j=1:N+1
            k=j+(i-1)*(N+1); %Number of vector component
            if j==1                
                uxchi(k)=U_wall_left*chi(j+i*(N+1));                
                uychi(k)=V_wall_bot*chi((i+1)+(j-1)*(N+1));
            elseif j==N+1
                uxchi(k)=U_wall_right*chi(j+i*(N+1));
                uychi(k)=V_wall_top*chi((i+1)+(j-1)*(N+1));
            else
                uxchi(k)=(u(j-1+(i-1)*(N-1))+u(j-1+i*(N-1)))/(2*h(j))*chi(j+i*(N+1));
                uychi(k)=(u(N*(N-1)+i+(j-2)*N)+u(N*(N-1)+i+1+(j-2)*N))/(2*h(j))*chi((i+1)+(j-1)*(N+1));
            end
        end
    end
   
    %Given vectors uxchi and uychi, vector convective can be constructed
    convective=zeros(2*N*(N-1),1);
    for i=1:N
        for j=1:N-1
            convective(j+(i-1)*(N-1))=-h(j+1)/2*(uychi(i+(j-1)*(N+1))+uychi(i+1+(j-1)*(N+1))); %Components along horizontal cell lines
            convective(N*(N-1)+i+(j-1)*N)=h(j+1)/2*(uxchi(i+(j-1)*(N+1))+uxchi(i+1+(j-1)*(N+1))); %Components along vertical cell lines
        end
    end
    
    % Set up the right hand side for the Poisson equation for the pressure
 



    rhs_Poisson  =   DIV*(u/dt  - convective - VLaplace*u/Re - u_pres/Re) + u_norm/dt; 
    
    % Solve for the new pressure
    
    temp = L\rhs_Poisson;
    p = U\temp;
    
    % Store the velocity from the previous time step in the vector u_old
    
    uold = u;
    
    % Udate the velocity field
    
    u = u - dt* (convective - tE21'*p + VLaplace*u/Re + u_pres/Re); 
    
    %
    %  Every other 1000 iterations check whether you approach steady state
    %  and check whether yopu satisfy conservation of mass. The largest
    %  rate at which mass is destroyed or created is denoted by 'maxdiv'.
    %  This number should be very small, in the order of machine precision.
    
    if mod(iter,1000) == 0
    
        maxdiv = max(DIV*u + u_norm)
        
        diff = max(abs(u-uold))/dt
        
    end
    iter = iter + 1;
    t=t+dt;

end

%% Plots of results
%Counter is an index that indicates at wich index the velocity component of
%x=0.5 is located in the vector u, i.e midpoint of the vector. It is
%usefull when plotting.

chi=Ht02*E21*u+u_pres_tvort;

if N==3
    counter=1;
elseif mod(N,2)~=0
    counter=(N-1)/2;
else
    counter=N/2-1;
end

% Extracting usefull results. 
%ux contains velocities in x direction
ux=u(1:length(u)/2);
%uy containing velocity in y direction
uy=u(length(u)/2+1:length(u)); 
%Convert circulations to actual velocities by dividing by the length
%segments

for i=0:N-1 %for velocities in x direction
    for j=1:length(h)-2
    ux(i*(N-1)+j)=ux(i*(N-1)+j)/h(j+1);
    end
end

for i=0:N-2 %for velocities in y direction
    for j=1:N
        uy(i*N+j)=uy(i*N+j)/h(i+2);
    end
end


%Find the vector containing the pressure Pres. This vector corresponds to
%p-(1/2)*(u^2+v^2)
  pres=p;

%Substract x velocities from p

ctemp=1;
for i=1:length(p)
    if mod(i,N)~=0
        pres(i,1)=pres(i,1)-0.5*(ux(ctemp))^2;
        disp(i);
        ctemp=ctemp+1;
    end
    
end

%substract y velocities

for i=1:length(uy)
    pres(i,1)=pres(i,1)-0.5*((uy(i))^2); 
end

% Set the pressure at the middle point equal to zero, and substract the
% precomputed value from the pressure vector, just as was done in the
% article by Botella.

if mod(N/2,2)==0
    pres=pres-(1/4)*(pres(counter*N+N/2)+pres(counter*N+N/2+1)+pres((counter+1)*N+N/2)+pres((counter+1)*N+N/2+1));
else
    pres=pres-pres(counter*N+counter);
    
end



%Plot of the pressure along the line y=0.5;
presMat=reshape(pres,[N,N])';
figure(1)
if mod(N,2)==0
    plot(x(2:length(x)-1),0.5*(presMat(N/2,:)+presMat(N/2+1,:)));
else
    plot(x(2:length(x)-1),pres(counter*N+1:counter*N+N));
end
title('Pressure along midline (y=0.5)');
xlabel('x')
ylabel('p')
hold off

%Plot of the vertical velocity along the line y=0.5

figure(2)
uyly05=zeros(length(x),1);
for i=2:length(x)-1
    uyly05(i,1)=uy(counter*N+i-1);
end
plot(x,uyly05);
title('Velocity in y direction along y=0.5');
xlabel('x');
ylabel('v');
hold off

%Plot of horizontal velocity along y=0.5
uxly05=zeros(length(x),1);
for i=2:length(x)-2
    uxly05(i,1)=ux(counter*(N-1)+i-1);
end
figure(3)
plot(x,uxly05);
title('Velocity in x direction along y=0.5');
xlabel('x');
ylabel('u');
hold off

%Plot of vorticity along y=0.5
vortly05=zeros(N+1,1);
for i=1:N+1
    vortly05(i,1)=chi((N+1)*(N/2)+i);
end
figure(4)
plot(tx,vortly05);
title('Vorticity along y=0.5');
xlabel('x');
ylabel('vorticity')



%Plot of velocity in y direction along x=0.5
uylx05=zeros(length(x),1);
for i=2:length(x)-2
    uylx05(i,1)=ux(counter+(i-2)*N);
end
figure(5)
plot(x,uyly05);
title('Velocity in y direction along x=0.5');
xlabel('y');
ylabel('v');
hold off



%plot of velocity in x direction along x=0.5
uxlx05=zeros(length(x),1);
uxlx05(length(x),1)=U_wall_top;
for i=2:length(x)-1
    uxlx05(i,1)=ux(counter+(i-2)*(N-1));
end
figure(6)
plot(x,uxlx05);
title('Velocity in x direction along x=0.5')
xlabel('y')
ylabel('u')
hold off

%Plot of pressure along x=0.5

% Extracting pressure value along x=0.5
if mod(N,2)==0
    pres_lx05=0.5*(presMat(:,N/2)+presMat(:,N/2+1));
else
    pres_lx05=presMat(:,(N-1)/2+1);
end
figure(7)
plot(x(2:length(x)-1),pres_lx05);
title('Pressure along x=0.5');
xlabel('y');
ylabel('x');
hold off

%Plot of vorticity along x=0.5

vortMat=reshape(chi,[N+1,N+1])';
if mod(N+1,2)==0
    vort_lx05=0.5*(vortMat(:,(N+1)/2)+vortMat(:,(N+1)/2+1));
else
    vort_lx05=vortMat(:,(N)/2+1);
end
figure(8)
plot(tx,vort_lx05);
title('Vorticity along x=0.5');
xlabel('y');
ylabel('Vorticity');
hold off

%Plot of pressure contour

figure(9) %Pressure contour
presCon=[0.3,0.17,0.12,0.11,0.09,0.07,0.05,0.02,0,-0.002];
Cpres=contour(x(2:length(x)-1),x(2:length(x)-1),presMat,presCon);
clabel(Cpres,presCon);
%['a','b','c','d','e','f','h','i','j']
title('Isobaric lines of flow')
hold off


%Vorticity Contour 
figure(10) %Vorticity contour plot
vortCon=[5,4,3,2,1,0.5,0,-0.5,-1,-2,-3];
Cvort=contour(tx,tx,vortMat,vortCon);
clabel(Cvort,vortCon);
title('Iso vorticity lines')
colorbar
xlabel('x')
ylabel('y');
hold off

%Streamfunction Contour 
%Preparing
tu=Ht11*u;
tu_ux_matrix=reshape(tu(1:(N-1)*N),[N-1,N])'; %u fluxes in one matrix
tu_ux_zero=zeros(N,1);

full_tu_ux_matrix=[tu_ux_zero,tu_ux_matrix,tu_ux_zero]; 

%Since the difference in value of the stream function representing two
%streamlines is the volume flux between them (assuming incompresibility),
%one can assume the value zero for the streamfunction along the bottom of
%the cavity (stagnation streamline). The value of the streamfunction on the
% 0 cochain one row up from this this line, is hence the flux crossing the cell edge
% connecting the two points. The process of adding the flux between
% concexutice points can then be repeated up to the top of the cavity.
streamFun=zeros(1,N+1); 
for i=1:N
    streamFun=[streamFun;streamFun(i,:)+full_tu_ux_matrix(i,:)];
end

figure(11)
streamCon=[0.1175, 0.115 ,0.11, 0.1 ,9*10^(-2), 7*10^(-2),5*10^(-2),3*10^(-2),1*10^(-2),10^(-4),10^(-5),10^(-10),0,-10^(-6),...
    10^(-5),5*10^(-5),-10^(-4),-2.5*10^(-4),-5*10^(-4),-10^(-3),-1.5*10^(-3)];
Cstream=contour(tx,tx,streamFun,streamCon);
clabel(Cstream,streamCon);
title('Streamlines');
xlabel('x');
ylabel('y');



%% Summing up the intergrated vorticity
Total_vorticity=sum(E21*u+Integrated_vor_pre);
disp(Total_vorticity);

%% Additional tasks 

%Setting up E10 and proving that matrix is equal to -tE21 transpose

%coresponding to N in this special case
E10=sparse(zeros(2*N*(N-1),N*N));
%Setting up upper half corresponding to 
for i=1:N %for each upper half row 
    for j=1:N-1
        E10((i-1)*(N-1)+j,(i-1)*N+j)=-1;
        E10((i-1)*(N-1)+j,(i-1)*N+j+1)=1;
    end
end

%Setting up lower half corresponding to v fluxes
for i = N*(N-1)+1:2*N*(N-1)
    E10(i,i-N*(N-1))=-1;
    E10(i,i-N*(N-1)+N)=1;
end


%Show equality between E10 and -tE21 transpose
isEqualE10=E10==(-tE21');
if sum(sum(isEqualE10))==size(E10,1)*size(E10,2)
   disp("E10 is equal to tE21'");
end


%% Showing the symetry and singularity of pressure matrix A

% if A is symetric, then it is equal to its transpose 

if A==A'
    disp('A is symetric, i.e Aij=Aji for all entries in matrix A');
end

% Since row sum equal to zero implies singularity, singularity can be
% proven by summing up the rows of A, and then sum up the sum of the rows
% of A and show that this sum is equal to  zero (or extremly close to due to
% rounding errors) 

disp('Sum of rowsums of pressure matrix A:')
disp(sum(sum(A')));

