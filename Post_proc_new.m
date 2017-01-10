clc; clear;

N_loc=importdata('domain.txt');
N_loc=N_loc(:,2:4);


Nx_tot=sum(N_loc([1, 10, 19],1));
Ny_tot=sum(N_loc([1, 4, 7],2));
Nz_tot=sum(N_loc(1:3,3));

%i_min i_max j_min j_max k_min k_max
Dim=zeros(27,6);
Dim(:,1:2:5)=1;
Dim(:,2:2:6)= N_loc;

domL=1:9;
domR=19:27;
domD=[1,2,3,10,11,12,19,20,21];
domU=domD+6;
domB=1:3:25;
domF=domB+2;

% Dim(domL,1)=Dim(domL,1)+1;
% Dim(domR,2)=Dim(domR,2)-1;
% Dim(domD,3)=Dim(domD,3)+1;
% Dim(domU,4)=Dim(domU,4)-1;
% Dim(domB,5)=Dim(domB,5)+1;
% Dim(domF,6)=Dim(domF,6)-1;

ptot=NaN(Nx_tot, Ny_tot, Nz_tot);

fid = fopen('pressure.txt');
P = fscanf(fid, '%f');%, N_loc(1,1)*N_loc(1,2)*N_loc(1,3));
fclose(fid);

count=0;

for p=1:27
    for k=Dim(p,5):Dim(p,6)
        for j=Dim(p,3):Dim(p,4)
            for i=Dim(p,1):Dim(p,2)
                
                count=count+1;
                %%
                if(p/3>1 && mod(p,9)~=1)
                    if (p>18)% && mod(p,6)==1)
                        appo=p-18;
                    elseif (p>9)
                        appo=p-9;
                    else
                        appo=p;
                    end
                    
                    if (appo>6)
                        ind_j=j+Dim(p-3,4)+Dim(p-6,4);
                    elseif (appo>3)
                        ind_j=j+Dim(p-3,4);
                    else
                        ind_j=j;
                    end
                else
                    ind_j=j;
                end
                
                if (p>9)
                    if(p<=18)
                        ind_i=i+Dim(p-9,2);
                    else
                        ind_i=i+Dim(p-9,2)+Dim(p-18,2);
                    end
                else
                    ind_i=i;
                end
                
                if  (mod(p,3)==1)
                    ind_k=k;
                else
                    if( mod(p,3)==2)
                        ind_k=k+Dim(p-1,6);
                    end
                    if( mod(p,3)==0)
                        ind_k=k+Dim(p-1,6)+Dim(p-2,6);
                    end
                end
                %%
                ptot(ind_i,ind_j,ind_k)=P(count);
                
            end
        end
    end
end

[row,col,pl]=find(isnan(ptot));


%La condizione alla dirichlet è imposta sulla faccia BACK
% ptot(:,:,1)

% BcD=1;
% BcN=0;
%
% ptot(:,:,1)=BcD;
%
% ptot(:,:,end)=ptot(:,:,end-1);
% ptot(:,1,:)=ptot(:,2,:);
% ptot(:,end,:)=ptot(:,end-1,:);
% ptot(1,:,:)=ptot(2,:,:);
% ptot(end,:,:)=ptot(end-1,:,:);

%A=ptot(5,:,:);
A(:,:)=ptot(10,:,:);
B(:,:)=ptot(:,10,:);
C(:,:)=ptot(:,:,10);

[G1, G2]=meshgrid(linspace(0,1,length(A(1,:))),linspace(0,1,length(A(1,:))));

figure
subplot(2,2,1)
surfc(A), pbaspect([1 1 0.5]), xlabel('y'), ylabel('z'), zlabel('p'), title('Sez. x = cost.')
subplot(2,2,2)
surfc(B), pbaspect([1 1 0.5]), xlabel('x'), ylabel('z'), zlabel('p'), title('Sez. y = cost.')
subplot(2,2,[3,4])
surfc(C), pbaspect([1 1 0.5]), xlabel('x'), ylabel('y'), zlabel('p'), title('Sez. z = cost.')


