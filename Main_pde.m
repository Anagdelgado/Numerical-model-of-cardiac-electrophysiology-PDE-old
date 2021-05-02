clear all; close all; clc;
tic

constants

%% ---> Spatial discretization
elem = reference_element();
load('pde_electrophysiology.mat')
load('mesh5.mat')
nng = length(X(:,1)); %---> Number of global nodes (nng)
[ne, nne] = size(R); %---> Number of elements, Number of nodes per element(nne=8)(integration points)

%% ---> Time discretization
n_periodes = 1;
Time = 1;
dt = 0.1;
dtt = 0.1;

tEnd    = Time*n_periodes;     
nStep   = tEnd/dt;      %---> Number of global time steps
nnStep  = tEnd/dtt;     %---> Number of intern time steps

%% ---> Initial conditions
V0 = zeros(nng,1);
%leftside = find(abs(X(:,1)) < 1e-5);
V0(:,1)= -86;
c_na = 11.6; c_k = 138.3; c_ca = 0.08*10^(-3); c_srca = 0.56;
g_m = 0; g_h = 0.75; g_j = 0.75; g_d = 0; g_f = 1; g_fca = 1; g_r = 0; g_s = 1; g_xs = 0; g_xr1 = 0; g_xr2 = 0; g_k1 = 0.5; g_g = 1;

%% ---> Iterative process
% Potential vector
V = zeros(nng,nStep); V(:,1)= V0; V_ant = V(:,1);

% Gates vectors
gant_m = g_m.*ones(ne,nne);
gant_h = g_h.*ones(ne,nne);
gant_j = g_j.*ones(ne,nne);
gant_d  = g_d.*ones(ne,nne);
gant_f = g_f.*ones(ne,nne);
gant_fca = g_fca.*ones(ne,nne);
gant_r = g_r.*ones(ne,nne);
gant_s = g_s.*ones(ne,nne);
gant_xs = g_xs.*ones(ne,nne);
gant_xr1 = g_xr1.*ones(ne,nne);
gant_xr2 = g_xr2.*ones(ne,nne);
gant_g = g_g.*ones(ne,nne);

% Concentration vectors
cant_k = c_k.*ones(ne,nne);
cant_na = c_na.*ones(ne,nne);
cant_ca = c_ca.*ones(ne,nne);
cant_srca = c_srca.*ones(ne,nne);

% Save results
t = zeros (1 ,nStep/10 + 1);
V_vtk = zeros(nng, nStep/10 + 1);
cont = 2;

% time newton iterations 
for n = 2:nStep + 1
    n
    I_stim = stim(X,n,nng,dt);
    
    % loop elements
    [V_nou,cnou_k,cnou_na,cnou_ca,cnou_srca,...
    gnou_m,gnou_h,gnou_j,gnou_xr1,gnou_xr2,gnou_xs,...
    gnou_r,gnou_s,gnou_d,gnou_f,gnou_fca,gnou_g] = newton_iter_time(X,R,elem,V_ant,dt,dtt,nng,ne,nne,nnStep,I_stim,...
                                                                    cant_k,cant_na,cant_ca,cant_srca,...
                                                                    gant_m,gant_h,gant_j,gant_xr1,gant_xr2,gant_xs,...
                                                                    gant_r,gant_s,gant_d,gant_f,gant_fca,gant_g);              
    V(:,n) = V_nou
    if(mod(n-1,10) < 1e-3) 
        V_vtk(:,cont) = V_nou;
        t(cont) = n*dt;
        cont = cont+1;
    end
    
    cnou_k
    cnou_na
    cnou_ca 
    cnou_srca
        
  
    % Update old (ant) variables
    V_ant = V_nou;
    cant_k = cnou_k;
    cant_na = cnou_na;
    cant_ca = cnou_ca;
    cant_srca = cnou_srca;
    gant_m = gnou_m;
    gant_h = gnou_h;
    gant_j = gnou_j;
    gant_xr1 = gnou_xr1;
    gant_xr2 = gnou_xr2;
    gant_xs = gnou_xs;
    gant_r = gnou_r;
    gant_s = gnou_s;
    gant_d = gnou_d;
    gant_f = gnou_f;
    gant_fca = gnou_fca;
    gant_g = gnou_g;
    
    
end

save('pde_electrophysiology.mat', 'V_vtk', 't')
vtk_results

figure(2)
plot(V(1,:)), hold on;