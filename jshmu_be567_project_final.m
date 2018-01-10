%%AIM 1a: Arabinose Inducible Promoter
clear all
close all

%initialize values, initial conditions of molecule states
mol = [1000 0 0 0 0 0 0 0]; 
ara=mol(1);
a_m=mol(2);
b_m=mol(3);
cde_m=mol(4);
a_p=mol(5);
b_p=mol(6);
cde_p=mol(7);
luc=mol(8);

%defining reaction rate
rate = [3.054 0.1386 0.1386 0.1386 2.0016 2.0016 2.0016 0.000578 0.000578 0.000578 0.1695 .00385 .0693]; %rates all in molecules/minute
r.lux_mtran=rate(1); %luxABCDE transcription rate
r.A_mde=rate(2); %luxA mRNA degradation rate
r.B_mde=rate(3); %luxB mRNA degradation rate
r.CDE_mde=rate(4); %luxCDE mRNA degradation rate
r.A_ptran=rate(5); %luxA mRNA translation rate
r.B_ptran=rate(6); %luxB mRNA translation rate
r.CDE_ptran=rate(7); %luxCDE mRNA translation rate
r.A_pde=rate(8); %luxA protein degradation rate
r.B_pde=rate(9); %luxB protein degradation rate
r.CDE_pde=rate(10); %luxCDE protein degradation rate
r.AB_dim=rate(11); %luxA and luxB dimerization rate
r.AB_deg=rate(12); %luxAB dimer degradation/dissocation rate
r.ara_deg=rate(13); %arabinose degradation rate

%propensity functions
prop=@(m,r)[r.lux_mtran*m(1) r.A_mde*m(2) r.B_mde*m(3) r.CDE_mde*m(4) r.A_ptran*m(2) r.B_ptran*m(3) r.CDE_ptran*m(4) r.A_pde*m(5) r.B_pde*m(6) r.CDE_pde*m(7) r.AB_dim*(m(5)*m(6)*m(7)) r.AB_deg*m(8) r.ara_deg*m(1)];

%stoich matrix
stm=[0	0	0	0	0	0	0	0	0	0	0	0 -1; 1	-1	0	0	0	0	0	0	0	0	0	0 0; 1	0	-1	0	0	0	0	0	0	0	0	0 0; 1	0	0	-1	0	0	0	0	0	0	0	0 0; 0	0	0	0	1	0	0	-1	0	0	-1	0 0; 0	0	0	0	0	1	0	0	-1	0	-1	0 0; 0	0	0	0	0	0	1	0	0	-1	0	0 0; 0	0	0	0	0	0	0	0	0	0	1	-1 0];

runs=5; %put how many runs of simulation you want
step=50; %put how many timesteps you want each run to have
end_time=360; %how many seconds you want sim to run for
%720 minutes chosen because want to run simulation for 6 h

peak_luc=zeros(runs,1);
peak_luc_time=zeros(runs,1);
end_luc=zeros(runs,1);
cum_luc=zeros(runs,1);

t_runs=zeros(runs,step);

out=linspace(0,end_time,step); %output trajectory
t=0;

for k=1:length(end_luc)
    t=0; %initial time
    j=1; %index
    times=zeros(1,length(out)); %matrix to collect time snapshots
    state=mol(:); %system state
    results=zeros(length(mol),length(out)); %setup output matrix
    results(:,1)=mol(:); % initialize reactants
    while t<out(end)
           %prop=[mol(1)*lambda mol(2)*gamma mol(2)*mu mol(3)*delta mol(3)*mu_p mol(4)*delta_p];
           alpha=feval(prop,state,r); %calculates propensity
           alpha_o=sum(alpha); %total
           tau=-(1/alpha_o)*log(rand()); %determing time interval tau using random number
           i=sum(rand >= cumsum([0 alpha]/alpha_o));
           state = state + stm(:,i);
           times(j)=t+tau; %store times that results correspond to
           t_runs=times(1,:); %stores times for each run
           t=t+tau;
           while j<=step && t>out(j)
               results(:,j)=state;
               j = j+1;
           end
           t_runs=times(1,:); %stores times for each run
           [row col]=max(results(8,:)); %stores peak luciferase dimensions
           peak_luc(k)=row; %stores peak luciferase values
           peak_luc_time(k)=t_runs(col); %stores time of peak luciferase 
           end_luc(k)=results(8,step-1); %stores amount of luciferase after each simulation run
           cum_luc(k)=sum(results(8,:)); %stores cumulative luciferace output for each simulation run
           t=t+tau;
    end
    k=k+1;
end

results(:,end)=state; %stores the results of the last run
t_runs(:,end)=end_time; %inputs the last time at the end_time

%mean, SD, CV of luc
format shortG
end_stats=[mean(end_luc) std(end_luc) std(end_luc)/mean(end_luc)]
peak_luc_stats=[mean(peak_luc) std(peak_luc) std(peak_luc)/mean(peak_luc)]
peak_time_stats= [mean(peak_luc_time) std(peak_luc_time) std(peak_luc_time)/mean(peak_luc_time)]./60
cum_luc_stats=[mean(cum_luc) std(cum_luc) std(cum_luc)/mean(cum_luc)]

figure
subplot(3,3,1)
plot(out./60,results(1,:))
xlabel('Time (hours)');
ylabel('Number arabinose molecules')

subplot(3,3,2)
plot(out./60,results(2,:))
xlabel('Time (hours)');
ylabel('Number luxA mRNA')

subplot(3,3,3)
plot(out./60,results(3,:))
xlabel('Time (hours)');
ylabel('Number luxB mRNA')

subplot(3,3,4)
plot(out./60,results(4,:))
xlabel('Time (hours)');
ylabel('Number luxCDE mRNA')

subplot(3,3,5)
plot(out./60,results(5,:))
xlabel('Time (hours)');
ylabel('Number luxA proteins')

subplot(3,3,6)
plot(out./60,results(6,:))
xlabel('Time (hours)');
ylabel('Number luxB proteins')

subplot(3,3,7)
plot(out./60,results(7,:))
xlabel('Time (hours)');
ylabel('Number luxCDE proteins')

subplot(3,3,8)
plot(out./60,results(8,:))
xlabel('Time (hours)');
ylabel('Number luciferase dimers')

%% AIM 1b: Constitutive Promoter
clear all
close all

%initialize values, initial conditions of molecule states
mol = [1 0 0 0 0 0 0 0]; %no arabinose molecules
active_gene=mol(1); %one active copy of gene instead of arabinose molecules
a_m=mol(2);
b_m=mol(3);
cde_m=mol(4);
a_p=mol(5);
b_p=mol(6);
cde_p=mol(7);
luc=mol(8);

%defining reaction rate
rate = [3.054 0.1386 0.1386 0.1386 2.0016 2.0016 2.0016 0.000578 0.000578 0.000578 0.1695 .00385]; %rates all in molecules/minute
r.lux_mtran=rate(1); %luxABCDE transcription rate
r.A_mde=rate(2); %luxA mRNA degradation rate
r.B_mde=rate(3); %luxB mRNA degradation rate
r.CDE_mde=rate(4); %luxCDE mRNA degradation rate
r.A_ptran=rate(5); %luxA mRNA translation rate
r.B_ptran=rate(6); %luxB mRNA translation rate
r.CDE_ptran=rate(7); %luxCDE mRNA translation rate
r.A_pde=rate(8); %luxA protein degradation rate
r.B_pde=rate(9); %luxB protein degradation rate
r.CDE_pde=rate(10); %luxCDE protein degradation rate
r.AB_dim=rate(11); %luxA and luxB dimerization rate
r.AB_deg=rate(12); %luxAB dimer degradation/dissocation rate

%propensity functions
prop=@(m,r)[r.lux_mtran*m(1) r.A_mde*m(2) r.B_mde*m(3) r.CDE_mde*m(4) r.A_ptran*m(2) r.B_ptran*m(3) r.CDE_ptran*m(4) r.A_pde*m(5) r.B_pde*m(6) r.CDE_pde*m(7) r.AB_dim*(m(5)*m(6)*m(7)) r.AB_deg*m(8)];

%stoich matrix
stm=[0	0	0	0	0	0	0	0	0	0	0	0; 1	-1	0	0	0	0	0	0	0	0	0	0; 1	0	-1	0	0	0	0	0	0	0	0	0; 1	0	0	-1	0	0	0	0	0	0	0	0; 0	0	0	0	1	0	0	-1	0	0	-1	0; 0	0	0	0	0	1	0	0	-1	0	-1	0; 0	0	0	0	0	0	1	0	0	-1	0	0; 0	0	0	0	0	0	0	0	0	0	1	-1 ];

runs=10; %put how many runs of simulation you want
step=100; %put how many timesteps you want each run to have
end_time=360; %how many seconds you want sim to run for
%720 minutes chosen because want to run simulation for 6 h

peak_luc=zeros(runs,1);
peak_luc_time=zeros(runs,1);
end_luc=zeros(runs,1);
cum_luc=zeros(runs,1);

t_runs=zeros(runs,step);

out=linspace(0,end_time,step); %output trajectory
t=0;

for k=1:length(end_luc)
    t=0; %initial time
    j=1; %index
    times=zeros(1,length(out)); %matrix to collect time snapshots
    state=mol(:); %system state
    results=zeros(length(mol),length(out)); %setup output matrix
    results(:,1)=mol(:); % initialize reactants
    while t<out(end)
           %prop=[mol(1)*lambda mol(2)*gamma mol(2)*mu mol(3)*delta mol(3)*mu_p mol(4)*delta_p];
           alpha=feval(prop,state,r); %calculates propensity
           alpha_o=sum(alpha); %total
           tau=-(1/alpha_o)*log(rand()); %determing time interval tau using random number
           i=sum(rand >= cumsum([0 alpha]/alpha_o));
           state = state + stm(:,i);
           times(j)=t+tau; %store times that results correspond to
           t_runs=times(1,:); %stores times for each run
           t=t+tau;
           while j<=step && t>out(j)
               results(:,j)=state;
               j = j+1;
           end
            t_runs=times(1,:); %stores times for each run
           [row col]=max(results(8,:)); %stores peak luciferase dimensions
           peak_luc(k)=row; %stores peak luciferase values
           peak_luc_time(k)=t_runs(col); %stores time of peak luciferase 
           end_luc(k)=results(8,step-1); %stores amount of luciferase after each simulation run
           cum_luc(k)=sum(results(8,:)); %stores cumulative luciferace output for each simulation run
           t=t+tau;
    end
    k=k+1;
end

results(:,end)=state; %stores the results of the last run
t_runs(:,end)=end_time; %inputs the last time at the end_time

%mean, SD, CV of luc
format shortG
end_stats=[mean(end_luc) std(end_luc) std(end_luc)/mean(end_luc)]
peak_luc_stats=[mean(peak_luc) std(peak_luc) std(peak_luc)/mean(peak_luc)]
peak_time_stats= [mean(peak_luc_time) std(peak_luc_time) std(peak_luc_time)/mean(peak_luc_time)]./60
cum_luc_stats=[mean(cum_luc) std(cum_luc) std(cum_luc)/mean(cum_luc)]

results(:,end)=state; %stores the results of the last run
t_runs(:,end)=end_time; %inputs the last time at the end_time

figure
subplot(3,3,1)
plot(out./60,results(1,:))
xlabel('Time (hours)');
ylabel('Number active genes')

subplot(3,3,2)
plot(out./60,results(2,:))
xlabel('Time (hours)');
ylabel('Number luxA mRNA')

subplot(3,3,3)
plot(out./60,results(3,:))
xlabel('Time (hours)');
ylabel('Number luxB mRNA')

subplot(3,3,4)
plot(out./60,results(4,:))
xlabel('Time (hours)');
ylabel('Number luxCDE mRNA')

subplot(3,3,5)
plot(out./60,results(5,:))
xlabel('Time (hours)');
ylabel('Number luxA proteins')

subplot(3,3,6)
plot(out./60,results(6,:))
xlabel('Time (hours)');
ylabel('Number luxB proteins')

subplot(3,3,7)
plot(out./60,results(7,:))
xlabel('Time (hours)');
ylabel('Number luxCDE proteins')

subplot(3,3,8)
plot(out./60,results(8,:))
xlabel('Time (hours)');
ylabel('Number luciferase dimers')

%% AIM 2a: Split lux System - A/B ind, C/D/E const

%%AIM 1: Arabinose Inducible Promoter
clear all
close all

%initialize values, initial conditions of molecule states
mol = [1000 0 0 0 0 0 0 0 1]; %active gene added as last value in mol matrix
ara=mol(1);
a_m=mol(2);
b_m=mol(3);
cde_m=mol(4);
a_p=mol(5);
b_p=mol(6);
cde_p=mol(7);
luc=mol(8);
active_gene=mol(9);

%defining reaction rate
rate = [3.054 0.1386 0.1386 0.1386 2.0016 2.0016 2.0016 0.000578 0.000578 0.000578 0.1695 .00385 .0693 3.054]; %rates all in molecules/minute
r.AB_mtran=rate(1); %luxAB transcription rate
r.A_mde=rate(2); %luxA mRNA degradation rate
r.B_mde=rate(3); %luxB mRNA degradation rate
r.CDE_mde=rate(4); %luxCDE mRNA degradation rate
r.A_ptran=rate(5); %luxA mRNA translation rate
r.B_ptran=rate(6); %luxB mRNA translation rate
r.CDE_ptran=rate(7); %luxCDE mRNA translation rate
r.A_pde=rate(8); %luxA protein degradation rate
r.B_pde=rate(9); %luxB protein degradation rate
r.CDE_pde=rate(10); %luxCDE protein degradation rate
r.AB_dim=rate(11); %luxA and luxB dimerization rate
r.AB_deg=rate(12); %luxAB dimer degradation/dissocation rate
r.ara_deg=rate(13); %arabinose degradation rate
r.CDE_mtran=rate(14); %luxCDE transcription rate

%propensity functions
prop=@(m,r)[r.AB_mtran*m(1) r.A_mde*m(2) r.B_mde*m(3) r.CDE_mde*m(4) r.A_ptran*m(2) r.B_ptran*m(3) r.CDE_ptran*m(4) r.A_pde*m(5) r.B_pde*m(6) r.CDE_pde*m(7) r.AB_dim*(m(5)*m(6)*m(7)) r.AB_deg*m(8) r.ara_deg*m(1) r.CDE_mtran*m(9)];

%stoich matrix
stm=[0	0	0	0	0	0	0	0	0	0	0	0	-1	0; 1	-1	0	0	0	0	0	0	0	0	0	0	0	0; 1	0	-1	0	0	0	0	0	0	0	0	0	0	0; 0	0	0	-1	0	0	0	0	0	0	0	0	0	1; 0	0	0	0	1	0	0	-1	0	0	-1	0	0	0; 0	0	0	0	0	1	0	0	-1	0	-1	0	0	0; 0	0	0	0	0	0	1	0	0	-1	0	0	0	0; 0	0	0	0	0	0	0	0	0	0	1	-1	0	0; 0	0	0	0	0	0	0	0	0	0	0	0	0	0];

runs=5; %put how many runs of simulation you want
step=50; %put how many timesteps you want each run to have
end_time=360; %how many seconds you want sim to run for
%720 minutes chosen because want to run simulation for 6 h

peak_luc=zeros(runs,1);
peak_luc_time=zeros(runs,1);
end_luc=zeros(runs,1);
cum_luc=zeros(runs,1);

t_runs=zeros(runs,step);

out=linspace(0,end_time,step); %output trajectory
t=0;

for k=1:length(end_luc)
    t=0; %initial time
    j=1; %index
    times=zeros(1,length(out)); %matrix to collect time snapshots
    state=mol(:); %system state
    results=zeros(length(mol),length(out)); %setup output matrix
    results(:,1)=mol(:); % initialize reactants
    while t<out(end)
           %prop=[mol(1)*lambda mol(2)*gamma mol(2)*mu mol(3)*delta mol(3)*mu_p mol(4)*delta_p];
           alpha=feval(prop,state,r); %calculates propensity
           alpha_o=sum(alpha); %total
           tau=-(1/alpha_o)*log(rand()); %determing time interval tau using random number
           i=sum(rand >= cumsum([0 alpha]/alpha_o));
           state = state + stm(:,i);
           times(j)=t+tau; %store times that results correspond to
           t_runs=times(1,:); %stores times for each run
           t=t+tau;
           while j<=step && t>out(j)
               results(:,j)=state;
               j = j+1;
           end
            t_runs=times(1,:); %stores times for each run
           [row col]=max(results(8,:)); %stores peak luciferase dimensions
           peak_luc(k)=row; %stores peak luciferase values
           peak_luc_time(k)=t_runs(col); %stores time of peak luciferase 
           end_luc(k)=results(8,step-1); %stores amount of luciferase after each simulation run
           cum_luc(k)=sum(results(8,:)); %stores cumulative luciferace output for each simulation run
           t=t+tau;
    end
    k=k+1;
end

results(:,end)=state; %stores the results of the last run
t_runs(:,end)=end_time; %inputs the last time at the end_time

%mean, SD, CV of luc
format shortG
end_stats=[mean(end_luc) std(end_luc) std(end_luc)/mean(end_luc)]
peak_luc_stats=[mean(peak_luc) std(peak_luc) std(peak_luc)/mean(peak_luc)]
peak_time_stats= [mean(peak_luc_time) std(peak_luc_time) std(peak_luc_time)/mean(peak_luc_time)]./60
cum_luc_stats=[mean(cum_luc) std(cum_luc) std(cum_luc)/mean(cum_luc)]

results(:,end)=state; %stores the results of the last run
t_runs(:,end)=end_time; %inputs the last time at the end_time

figure
subplot(3,3,1)
plot(out./60,results(1,:))
xlabel('Time (hours)');
ylabel('Number arabinose molecules')

subplot(3,3,2)
plot(out./60,results(2,:))
xlabel('Time (hours)');
ylabel('Number luxA mRNA')

subplot(3,3,3)
plot(out./60,results(3,:))
xlabel('Time (hours)');
ylabel('Number luxB mRNA')

subplot(3,3,4)
plot(out./60,results(4,:))
xlabel('Time (hours)');
ylabel('Number luxCDE mRNA')

subplot(3,3,5)
plot(out./60,results(5,:))
xlabel('Time (hours)');
ylabel('Number luxA proteins')

subplot(3,3,6)
plot(out./60,results(6,:))
xlabel('Time (hours)');
ylabel('Number luxB proteins')

subplot(3,3,7)
plot(out./60,results(7,:))
xlabel('Time (hours)');
ylabel('Number luxCDE proteins')

subplot(3,3,8)
plot(out./60,results(8,:))
xlabel('Time (hours)');
ylabel('Number luciferase dimers')

subplot(3,3,9)
plot(out./60,results(9,:))
xlabel('Time (hours)');
ylabel('Number active genes')

%% AIM 3c: Split lux System - shifted induction
%induction after 2 hours
clear all
close all

%initialize values, initial conditions of molecule states
mol = [0 0 0 0 0 0 0 0 1]; %active gene added as last value in mol matrix
ara=mol(1);
a_m=mol(2);
b_m=mol(3);
cde_m=mol(4);
a_p=mol(5);
b_p=mol(6);
cde_p=mol(7);
luc=mol(8);
active_gene=mol(9);

%defining reaction rate
rate = [3.054 0.1386 0.1386 0.1386 2.0016 2.0016 2.0016 0.000578 0.000578 0.000578 0.1695 .00385 .0693 3.054]; %rates all in molecules/minute
r.AB_mtran=rate(1); %luxAB transcription rate
r.A_mde=rate(2); %luxA mRNA degradation rate
r.B_mde=rate(3); %luxB mRNA degradation rate
r.CDE_mde=rate(4); %luxCDE mRNA degradation rate
r.A_ptran=rate(5); %luxA mRNA translation rate
r.B_ptran=rate(6); %luxB mRNA translation rate
r.CDE_ptran=rate(7); %luxCDE mRNA translation rate
r.A_pde=rate(8); %luxA protein degradation rate
r.B_pde=rate(9); %luxB protein degradation rate
r.CDE_pde=rate(10); %luxCDE protein degradation rate
r.AB_dim=rate(11); %luxA and luxB dimerization rate
r.AB_deg=rate(12); %luxAB dimer degradation/dissocation rate
r.ara_deg=rate(13); %arabinose degradation rate
r.CDE_mtran=rate(14); %luxCDE transcription rate

%propensity functions
prop=@(m,r)[r.AB_mtran*m(1) r.A_mde*m(2) r.B_mde*m(3) r.CDE_mde*m(4) r.A_ptran*m(2) r.B_ptran*m(3) r.CDE_ptran*m(4) r.A_pde*m(5) r.B_pde*m(6) r.CDE_pde*m(7) r.AB_dim*(m(5)*m(6)*m(7)) r.AB_deg*m(8) r.ara_deg*m(1) r.CDE_mtran*m(9)];

%stoich matrix
stm=[0	0	0	0	0	0	0	0	0	0	0	0	-1	0; 1	-1	0	0	0	0	0	0	0	0	0	0	0	0; 1	0	-1	0	0	0	0	0	0	0	0	0	0	0; 0	0	0	-1	0	0	0	0	0	0	0	0	0	1; 0	0	0	0	1	0	0	-1	0	0	-1	0	0	0; 0	0	0	0	0	1	0	0	-1	0	-1	0	0	0; 0	0	0	0	0	0	1	0	0	-1	0	0	0	0; 0	0	0	0	0	0	0	0	0	0	1	-1	0	0; 0	0	0	0	0	0	0	0	0	0	0	0	0	0];

runs=5; %put how many runs of simulation you want
step=50; %put how many timesteps you want each run to have
end_time=360; %how many seconds you want sim to run for
%720 minutes chosen because want to run simulation for 6 h

peak_luc=zeros(runs,1);
peak_luc_time=zeros(runs,1);
end_luc=zeros(runs,1);
cum_luc=zeros(runs,1);

t_runs=zeros(runs,step);

out=linspace(0,end_time,step); %output trajectory
t=0;

for k=1:length(end_luc)
    t=0; %initial time
    j=1; %index
    times=zeros(1,length(out)); %matrix to collect time snapshots
    state=mol(:); %system state
    results=zeros(length(mol),length(out)); %setup output matrix
    results(:,1)=mol(:); % initialize reactants
    i=0; % counter to make sure induction only happens once
    while t<out(end)
           alpha=feval(prop,state,r);
           alpha_o=sum(alpha); %total
           tau=-(1/alpha_o)*log(rand()); %determing time interval tau using random number
           i=sum(rand >= cumsum([0 alpha]/alpha_o));
           state = state + stm(:,i);
           times(j)=t+tau; %store times that results correspond to
           t_runs=times(1,:); %stores times for each run
           t=t+tau;
           while j<=step && t>out(j)
               results(:,j)=state;
               j = j+1;
           end
           t_runs=times(1,:); %stores times for each run
           [row col]=max(results(8,:)); %stores peak luciferase dimensions
           peak_luc(k)=row; %stores peak luciferase values
           peak_luc_time(k)=t_runs(col); %stores time of peak luciferase 
           end_luc(k)=results(8,step-1); %stores amount of luciferase after each simulation run
           cum_luc(k)=sum(results(8,:)); %stores cumulative luciferace output for each simulation run
           t=t+tau;
           if t>120 & t<130 %small window after two-hour mark during which arabinose is added
               state(1) = 1000;
           end
    end
    k=k+1;
end

results(:,end)=state; %stores the results of the last run
t_runs(:,end)=end_time; %inputs the last time at the end_time

%mean, SD, CV of luc
format shortG
end_stats=[mean(end_luc) std(end_luc) std(end_luc)/mean(end_luc)]
peak_luc_stats=[mean(peak_luc) std(peak_luc) std(peak_luc)/mean(peak_luc)]
peak_time_stats= [mean(peak_luc_time) std(peak_luc_time) std(peak_luc_time)/mean(peak_luc_time)]./60
cum_luc_stats=[mean(cum_luc) std(cum_luc) std(cum_luc)/mean(cum_luc)]

results(:,end)=state; %stores the results of the last run
t_runs(:,end)=end_time; %inputs the last time at the end_time

figure
subplot(3,3,1)
plot(out./60,results(1,:))
xlabel('Time (hours)');
ylabel('Number arabinose molecules')

subplot(3,3,2)
plot(out./60,results(2,:))
xlabel('Time (hours)');
ylabel('Number luxA mRNA')

subplot(3,3,3)
plot(out./60,results(3,:))
xlabel('Time (hours)');
ylabel('Number luxB mRNA')

subplot(3,3,4)
plot(out./60,results(4,:))
xlabel('Time (hours)');
ylabel('Number luxCDE mRNA')

subplot(3,3,5)
plot(out./60,results(5,:))
xlabel('Time (hours)');
ylabel('Number luxA proteins')

subplot(3,3,6)
plot(out./60,results(6,:))
xlabel('Time (hours)');
ylabel('Number luxB proteins')

subplot(3,3,7)
plot(out./60,results(7,:))
xlabel('Time (hours)');
ylabel('Number luxCDE proteins')

subplot(3,3,8)
plot(out./60,results(8,:))
xlabel('Time (hours)');
ylabel('Number luciferase dimers')

subplot(3,3,9)
plot(out./60,results(9,:))
xlabel('Time (hours)');
ylabel('Number active genes')

%% AIM 3a: Split lux System - shifted induction
%induction after 1 hour
clear all
close all

%initialize values, initial conditions of molecule states
mol = [0 0 0 0 0 0 0 0 1]; %active gene added as last value in mol matrix
ara=mol(1);
a_m=mol(2);
b_m=mol(3);
cde_m=mol(4);
a_p=mol(5);
b_p=mol(6);
cde_p=mol(7);
luc=mol(8);
active_gene=mol(9);

%defining reaction rate
rate = [3.054 0.1386 0.1386 0.1386 2.0016 2.0016 2.0016 0.000578 0.000578 0.000578 0.1695 .00385 .0693 3.054]; %rates all in molecules/minute
r.AB_mtran=rate(1); %luxAB transcription rate
r.A_mde=rate(2); %luxA mRNA degradation rate
r.B_mde=rate(3); %luxB mRNA degradation rate
r.CDE_mde=rate(4); %luxCDE mRNA degradation rate
r.A_ptran=rate(5); %luxA mRNA translation rate
r.B_ptran=rate(6); %luxB mRNA translation rate
r.CDE_ptran=rate(7); %luxCDE mRNA translation rate
r.A_pde=rate(8); %luxA protein degradation rate
r.B_pde=rate(9); %luxB protein degradation rate
r.CDE_pde=rate(10); %luxCDE protein degradation rate
r.AB_dim=rate(11); %luxA and luxB dimerization rate
r.AB_deg=rate(12); %luxAB dimer degradation/dissocation rate
r.ara_deg=rate(13); %arabinose degradation rate
r.CDE_mtran=rate(14); %luxCDE transcription rate

%propensity functions
prop=@(m,r)[r.AB_mtran*m(1) r.A_mde*m(2) r.B_mde*m(3) r.CDE_mde*m(4) r.A_ptran*m(2) r.B_ptran*m(3) r.CDE_ptran*m(4) r.A_pde*m(5) r.B_pde*m(6) r.CDE_pde*m(7) r.AB_dim*(m(5)*m(6)*m(7)) r.AB_deg*m(8) r.ara_deg*m(1) r.CDE_mtran*m(9)];

%stoich matrix
stm=[0	0	0	0	0	0	0	0	0	0	0	0	-1	0; 1	-1	0	0	0	0	0	0	0	0	0	0	0	0; 1	0	-1	0	0	0	0	0	0	0	0	0	0	0; 0	0	0	-1	0	0	0	0	0	0	0	0	0	1; 0	0	0	0	1	0	0	-1	0	0	-1	0	0	0; 0	0	0	0	0	1	0	0	-1	0	-1	0	0	0; 0	0	0	0	0	0	1	0	0	-1	0	0	0	0; 0	0	0	0	0	0	0	0	0	0	1	-1	0	0; 0	0	0	0	0	0	0	0	0	0	0	0	0	0];

runs=5; %put how many runs of simulation you want
step=50; %put how many timesteps you want each run to have
end_time=360; %how many seconds you want sim to run for
%720 minutes chosen because want to run simulation for 6 h

peak_luc=zeros(runs,1);
peak_luc_time=zeros(runs,1);
end_luc=zeros(runs,1);
cum_luc=zeros(runs,1);

t_runs=zeros(runs,step);

out=linspace(0,end_time,step); %output trajectory
t=0;

for k=1:length(end_luc)
    t=0; %initial time
    j=1; %index
    times=zeros(1,length(out)); %matrix to collect time snapshots
    state=mol(:); %system state
    results=zeros(length(mol),length(out)); %setup output matrix
    results(:,1)=mol(:); % initialize reactants
    i=0; % counter to make sure induction only happens once
    while t<out(end)
           alpha=feval(prop,state,r);
           alpha_o=sum(alpha); %total
           tau=-(1/alpha_o)*log(rand()); %determing time interval tau using random number
           i=sum(rand >= cumsum([0 alpha]/alpha_o));
           state = state + stm(:,i);
           times(j)=t+tau; %store times that results correspond to
           t_runs=times(1,:); %stores times for each run
           t=t+tau;
           while j<=step && t>out(j)
               results(:,j)=state;
               j = j+1;
           end
           t_runs=times(1,:); %stores times for each run
           [row col]=max(results(8,:)); %stores peak luciferase dimensions
           peak_luc(k)=row; %stores peak luciferase values
           peak_luc_time(k)=t_runs(col); %stores time of peak luciferase 
           end_luc(k)=results(8,step-1); %stores amount of luciferase after each simulation run
           cum_luc(k)=sum(results(8,:)); %stores cumulative luciferace output for each simulation run
           t=t+tau;
           if t>60 & t<70 %small window after one-hour mark during which arabinose is added
               state(1) = 1000;
           end
    end
    k=k+1;
end

results(:,end)=state; %stores the results of the last run
t_runs(:,end)=end_time; %inputs the last time at the end_time

%mean, SD, CV of luc
end_stats=[mean(end_luc) std(end_luc) std(end_luc)/mean(end_luc)]
peak_luc_stats=[mean(peak_luc) std(peak_luc) std(peak_luc)/mean(peak_luc)]
peak_time_stats= [mean(peak_luc_time) std(peak_luc_time) std(peak_luc_time)/mean(peak_luc_time)]./60
cum_luc_stats=[mean(cum_luc) std(cum_luc) std(cum_luc)/mean(cum_luc)]

results(:,end)=state; %stores the results of the last run
t_runs(:,end)=end_time; %inputs the last time at the end_time

figure
subplot(3,3,1)
plot(out./60,results(1,:))
xlabel('Time (hours)');
ylabel('Number arabinose molecules')

subplot(3,3,2)
plot(out./60,results(2,:))
xlabel('Time (hours)');
ylabel('Number luxA mRNA')

subplot(3,3,3)
plot(out./60,results(3,:))
xlabel('Time (hours)');
ylabel('Number luxB mRNA')

subplot(3,3,4)
plot(out./60,results(4,:))
xlabel('Time (hours)');
ylabel('Number luxCDE mRNA')

subplot(3,3,5)
plot(out./60,results(5,:))
xlabel('Time (hours)');
ylabel('Number luxA proteins')

subplot(3,3,6)
plot(out./60,results(6,:))
xlabel('Time (hours)');
ylabel('Number luxB proteins')

subplot(3,3,7)
plot(out./60,results(7,:))
xlabel('Time (hours)');
ylabel('Number luxCDE proteins')

subplot(3,3,8)
plot(out./60,results(8,:))
xlabel('Time (hours)');
ylabel('Number luciferase dimers')

subplot(3,3,9)
plot(out./60,results(9,:))
xlabel('Time (hours)');
ylabel('Number active genes')

%% AIM 3d: Split lux System - shifted induction
%induction after 4 hours
clear all
close all

%initialize values, initial conditions of molecule states
mol = [0 0 0 0 0 0 0 0 1]; %active gene added as last value in mol matrix
ara=mol(1);
a_m=mol(2);
b_m=mol(3);
cde_m=mol(4);
a_p=mol(5);
b_p=mol(6);
cde_p=mol(7);
luc=mol(8);
active_gene=mol(9);

%defining reaction rate
rate = [3.054 0.1386 0.1386 0.1386 2.0016 2.0016 2.0016 0.000578 0.000578 0.000578 0.1695 .00385 .0693 3.054]; %rates all in molecules/minute
r.AB_mtran=rate(1); %luxAB transcription rate
r.A_mde=rate(2); %luxA mRNA degradation rate
r.B_mde=rate(3); %luxB mRNA degradation rate
r.CDE_mde=rate(4); %luxCDE mRNA degradation rate
r.A_ptran=rate(5); %luxA mRNA translation rate
r.B_ptran=rate(6); %luxB mRNA translation rate
r.CDE_ptran=rate(7); %luxCDE mRNA translation rate
r.A_pde=rate(8); %luxA protein degradation rate
r.B_pde=rate(9); %luxB protein degradation rate
r.CDE_pde=rate(10); %luxCDE protein degradation rate
r.AB_dim=rate(11); %luxA and luxB dimerization rate
r.AB_deg=rate(12); %luxAB dimer degradation/dissocation rate
r.ara_deg=rate(13); %arabinose degradation rate
r.CDE_mtran=rate(14); %luxCDE transcription rate

%propensity functions
prop=@(m,r)[r.AB_mtran*m(1) r.A_mde*m(2) r.B_mde*m(3) r.CDE_mde*m(4) r.A_ptran*m(2) r.B_ptran*m(3) r.CDE_ptran*m(4) r.A_pde*m(5) r.B_pde*m(6) r.CDE_pde*m(7) r.AB_dim*(m(5)*m(6)*m(7)) r.AB_deg*m(8) r.ara_deg*m(1) r.CDE_mtran*m(9)];

%stoich matrix
stm=[0	0	0	0	0	0	0	0	0	0	0	0	-1	0; 1	-1	0	0	0	0	0	0	0	0	0	0	0	0; 1	0	-1	0	0	0	0	0	0	0	0	0	0	0; 0	0	0	-1	0	0	0	0	0	0	0	0	0	1; 0	0	0	0	1	0	0	-1	0	0	-1	0	0	0; 0	0	0	0	0	1	0	0	-1	0	-1	0	0	0; 0	0	0	0	0	0	1	0	0	-1	0	0	0	0; 0	0	0	0	0	0	0	0	0	0	1	-1	0	0; 0	0	0	0	0	0	0	0	0	0	0	0	0	0];

runs=5; %put how many runs of simulation you want
step=50; %put how many timesteps you want each run to have
end_time=360; %how many seconds you want sim to run for
%720 minutes chosen because want to run simulation for 6 h

peak_luc=zeros(runs,1);
peak_luc_time=zeros(runs,1);
end_luc=zeros(runs,1);
cum_luc=zeros(runs,1);

t_runs=zeros(runs,step);

out=linspace(0,end_time,step); %output trajectory
t=0;

for k=1:length(end_luc)
    t=0; %initial time
    j=1; %index
    times=zeros(1,length(out)); %matrix to collect time snapshots
    state=mol(:); %system state
    results=zeros(length(mol),length(out)); %setup output matrix
    results(:,1)=mol(:); % initialize reactants
    i=0; % counter to make sure induction only happens once
    while t<out(end)
           alpha=feval(prop,state,r);
           alpha_o=sum(alpha); %total
           tau=-(1/alpha_o)*log(rand()); %determing time interval tau using random number
           i=sum(rand >= cumsum([0 alpha]/alpha_o));
           state = state + stm(:,i);
           times(j)=t+tau; %store times that results correspond to
           t_runs=times(1,:); %stores times for each run
           t=t+tau;
           while j<=step && t>out(j)
               results(:,j)=state;
               j = j+1;
           end
           t_runs=times(1,:); %stores times for each run
           [row col]=max(results(8,:)); %stores peak luciferase dimensions
           peak_luc(k)=row; %stores peak luciferase values
           peak_luc_time(k)=t_runs(col); %stores time of peak luciferase 
           end_luc(k)=results(8,step-1); %stores amount of luciferase after each simulation run
           cum_luc(k)=sum(results(8,:)); %stores cumulative luciferace output for each simulation run
           t=t+tau;
           if t>240 & t<250 %small window after four-hour mark during which arabinose is added
               state(1) = 1000;
           end
    end
    k=k+1;
end

results(:,end)=state; %stores the results of the last run
t_runs(:,end)=end_time; %inputs the last time at the end_time

%mean, SD, CV of luc
end_stats=[mean(end_luc) std(end_luc) std(end_luc)/mean(end_luc)]
peak_luc_stats=[mean(peak_luc) std(peak_luc) std(peak_luc)/mean(peak_luc)]
peak_time_stats= [mean(peak_luc_time) std(peak_luc_time) std(peak_luc_time)/mean(peak_luc_time)]./60
cum_luc_stats=[mean(cum_luc) std(cum_luc) std(cum_luc)/mean(cum_luc)]

results(:,end)=state; %stores the results of the last run
t_runs(:,end)=end_time; %inputs the last time at the end_time

figure
subplot(3,3,1)
plot(out./60,results(1,:))
xlabel('Time (hours)');
ylabel('Number arabinose molecules')

subplot(3,3,2)
plot(out./60,results(2,:))
xlabel('Time (hours)');
ylabel('Number luxA mRNA')

subplot(3,3,3)
plot(out./60,results(3,:))
xlabel('Time (hours)');
ylabel('Number luxB mRNA')

subplot(3,3,4)
plot(out./60,results(4,:))
xlabel('Time (hours)');
ylabel('Number luxCDE mRNA')

subplot(3,3,5)
plot(out./60,results(5,:))
xlabel('Time (hours)');
ylabel('Number luxA proteins')

subplot(3,3,6)
plot(out./60,results(6,:))
xlabel('Time (hours)');
ylabel('Number luxB proteins')

subplot(3,3,7)
plot(out./60,results(7,:))
xlabel('Time (hours)');
ylabel('Number luxCDE proteins')

subplot(3,3,8)
plot(out./60,results(8,:))
xlabel('Time (hours)');
ylabel('Number luciferase dimers')

subplot(3,3,9)
plot(out./60,results(9,:))
xlabel('Time (hours)');
ylabel('Number active genes')
%% AIM 2b: Split lux System opposite - A/B const, C/D/E ind
%immediate induction
clear all
close all

%initialize values, initial conditions of molecule states
mol = [1000 0 0 0 0 0 0 0 1]; %active gene added as last value in mol matrix
ara=mol(1);
a_m=mol(2);
b_m=mol(3);
cde_m=mol(4);
a_p=mol(5);
b_p=mol(6);
cde_p=mol(7);
luc=mol(8);
active_gene=mol(9);

%defining reaction rate
rate = [3.054 0.1386 0.1386 0.1386 2.0016 2.0016 2.0016 0.000578 0.000578 0.000578 0.1695 .00385 .0693 3.054]; %rates all in molecules/minute
r.AB_mtran=rate(1); %luxAB transcription rate
r.A_mde=rate(2); %luxA mRNA degradation rate
r.B_mde=rate(3); %luxB mRNA degradation rate
r.CDE_mde=rate(4); %luxCDE mRNA degradation rate
r.A_ptran=rate(5); %luxA mRNA translation rate
r.B_ptran=rate(6); %luxB mRNA translation rate
r.CDE_ptran=rate(7); %luxCDE mRNA translation rate
r.A_pde=rate(8); %luxA protein degradation rate
r.B_pde=rate(9); %luxB protein degradation rate
r.CDE_pde=rate(10); %luxCDE protein degradation rate
r.AB_dim=rate(11); %luxA and luxB dimerization rate
r.AB_deg=rate(12); %luxAB dimer degradation/dissocation rate
r.ara_deg=rate(13); %arabinose degradation rate
r.CDE_mtran=rate(14); %luxCDE transcription rate

%propensity functions
prop=@(m,r)[r.AB_mtran*m(9) r.A_mde*m(2) r.B_mde*m(3) r.CDE_mde*m(4) r.A_ptran*m(2) r.B_ptran*m(3) r.CDE_ptran*m(4) r.A_pde*m(5) r.B_pde*m(6) r.CDE_pde*m(7) r.AB_dim*(m(5)*m(6)*m(7)) r.AB_deg*m(8) r.ara_deg*m(1) r.CDE_mtran*m(1)];

%stoich matrix
stm=[0	0	0	0	0	0	0	0	0	0	0	0	-1	0; 1	-1	0	0	0	0	0	0	0	0	0	0	0	0; 1	0	-1	0	0	0	0	0	0	0	0	0	0	0; 0	0	0	-1	0	0	0	0	0	0	0	0	0	1; 0	0	0	0	1	0	0	-1	0	0	-1	0	0	0; 0	0	0	0	0	1	0	0	-1	0	-1	0	0	0; 0	0	0	0	0	0	1	0	0	-1	0	0	0	0; 0	0	0	0	0	0	0	0	0	0	1	-1	0	0; 0	0	0	0	0	0	0	0	0	0	0	0	0	0];

runs=5; %put how many runs of simulation you want
step=50; %put how many timesteps you want each run to have
end_time=360; %how many seconds you want sim to run for
%720 minutes chosen because want to run simulation for 6 h

peak_luc=zeros(runs,1);
peak_luc_time=zeros(runs,1);
end_luc=zeros(runs,1);
cum_luc=zeros(runs,1);

t_runs=zeros(runs,step);

out=linspace(0,end_time,step); %output trajectory
t=0;

for k=1:length(end_luc)
    t=0; %initial time
    j=1; %index
    times=zeros(1,length(out)); %matrix to collect time snapshots
    state=mol(:); %system state
    results=zeros(length(mol),length(out)); %setup output matrix
    results(:,1)=mol(:); % initialize reactants
    i=0; % counter to make sure induction only happens once
    while t<out(end)
           alpha=feval(prop,state,r);
           alpha_o=sum(alpha); %total
           tau=-(1/alpha_o)*log(rand()); %determing time interval tau using random number
           i=sum(rand >= cumsum([0 alpha]/alpha_o));
           state = state + stm(:,i);
           times(j)=t+tau; %store times that results correspond to
           t_runs=times(1,:); %stores times for each run
           t=t+tau;
           while j<=step && t>out(j)
               results(:,j)=state;
               j = j+1;
           end
           t_runs=times(1,:); %stores times for each run
           [row col]=max(results(8,:)); %stores peak luciferase dimensions
           peak_luc(k)=row; %stores peak luciferase values
           peak_luc_time(k)=t_runs(col); %stores time of peak luciferase 
           end_luc(k)=results(8,step-1); %stores amount of luciferase after each simulation run
           cum_luc(k)=sum(results(8,:)); %stores cumulative luciferace output for each simulation run
           t=t+tau;
    end
    k=k+1;
end

results(:,end)=state; %stores the results of the last run
t_runs(:,end)=end_time; %inputs the last time at the end_time

%mean, SD, CV of luc
format shortG
end_stats=[mean(end_luc) std(end_luc) std(end_luc)/mean(end_luc)]
peak_luc_stats=[mean(peak_luc) std(peak_luc) std(peak_luc)/mean(peak_luc)]
peak_time_stats= [mean(peak_luc_time) std(peak_luc_time) std(peak_luc_time)/mean(peak_luc_time)]./60
cum_luc_stats=[mean(cum_luc) std(cum_luc) std(cum_luc)/mean(cum_luc)]

results(:,end)=state; %stores the results of the last run
t_runs(:,end)=end_time; %inputs the last time at the end_time

figure
subplot(3,3,1)
plot(out./60,results(1,:))
xlabel('Time (hours)');
ylabel('Number arabinose molecules')

subplot(3,3,2)
plot(out./60,results(2,:))
xlabel('Time (hours)');
ylabel('Number luxA mRNA')

subplot(3,3,3)
plot(out./60,results(3,:))
xlabel('Time (hours)');
ylabel('Number luxB mRNA')

subplot(3,3,4)
plot(out./60,results(4,:))
xlabel('Time (hours)');
ylabel('Number luxCDE mRNA')

subplot(3,3,5)
plot(out./60,results(5,:))
xlabel('Time (hours)');
ylabel('Number luxA proteins')

subplot(3,3,6)
plot(out./60,results(6,:))
xlabel('Time (hours)');
ylabel('Number luxB proteins')

subplot(3,3,7)
plot(out./60,results(7,:))
xlabel('Time (hours)');
ylabel('Number luxCDE proteins')

subplot(3,3,8)
plot(out./60,results(8,:))
xlabel('Time (hours)');
ylabel('Number luciferase dimers')

subplot(3,3,9)
plot(out./60,results(9,:))
xlabel('Time (hours)');
ylabel('Number active genes')

%% AIM 3b: Split lux System - shifted induction
%induction after 1.8776 hours
clear all
close all

%initialize values, initial conditions of molecule states
mol = [0 0 0 0 0 0 0 0 1]; %active gene added as last value in mol matrix
ara=mol(1);
a_m=mol(2);
b_m=mol(3);
cde_m=mol(4);
a_p=mol(5);
b_p=mol(6);
cde_p=mol(7);
luc=mol(8);
active_gene=mol(9);

%defining reaction rate
rate = [3.054 0.1386 0.1386 0.1386 2.0016 2.0016 2.0016 0.000578 0.000578 0.000578 0.1695 .00385 .0693 3.054]; %rates all in molecules/minute
r.AB_mtran=rate(1); %luxAB transcription rate
r.A_mde=rate(2); %luxA mRNA degradation rate
r.B_mde=rate(3); %luxB mRNA degradation rate
r.CDE_mde=rate(4); %luxCDE mRNA degradation rate
r.A_ptran=rate(5); %luxA mRNA translation rate
r.B_ptran=rate(6); %luxB mRNA translation rate
r.CDE_ptran=rate(7); %luxCDE mRNA translation rate
r.A_pde=rate(8); %luxA protein degradation rate
r.B_pde=rate(9); %luxB protein degradation rate
r.CDE_pde=rate(10); %luxCDE protein degradation rate
r.AB_dim=rate(11); %luxA and luxB dimerization rate
r.AB_deg=rate(12); %luxAB dimer degradation/dissocation rate
r.ara_deg=rate(13); %arabinose degradation rate
r.CDE_mtran=rate(14); %luxCDE transcription rate

%propensity functions
prop=@(m,r)[r.AB_mtran*m(1) r.A_mde*m(2) r.B_mde*m(3) r.CDE_mde*m(4) r.A_ptran*m(2) r.B_ptran*m(3) r.CDE_ptran*m(4) r.A_pde*m(5) r.B_pde*m(6) r.CDE_pde*m(7) r.AB_dim*(m(5)*m(6)*m(7)) r.AB_deg*m(8) r.ara_deg*m(1) r.CDE_mtran*m(9)];

%stoich matrix
stm=[0	0	0	0	0	0	0	0	0	0	0	0	-1	0; 1	-1	0	0	0	0	0	0	0	0	0	0	0	0; 1	0	-1	0	0	0	0	0	0	0	0	0	0	0; 0	0	0	-1	0	0	0	0	0	0	0	0	0	1; 0	0	0	0	1	0	0	-1	0	0	-1	0	0	0; 0	0	0	0	0	1	0	0	-1	0	-1	0	0	0; 0	0	0	0	0	0	1	0	0	-1	0	0	0	0; 0	0	0	0	0	0	0	0	0	0	1	-1	0	0; 0	0	0	0	0	0	0	0	0	0	0	0	0	0];

runs=5; %put how many runs of simulation you want
step=50; %put how many timesteps you want each run to have
end_time=360; %how many seconds you want sim to run for
%720 minutes chosen because want to run simulation for 6 h

peak_luc=zeros(runs,1);
peak_luc_time=zeros(runs,1);
end_luc=zeros(runs,1);
cum_luc=zeros(runs,1);

t_runs=zeros(runs,step);

out=linspace(0,end_time,step); %output trajectory
t=0;

for k=1:length(end_luc)
    t=0; %initial time
    j=1; %index
    times=zeros(1,length(out)); %matrix to collect time snapshots
    state=mol(:); %system state
    results=zeros(length(mol),length(out)); %setup output matrix
    results(:,1)=mol(:); % initialize reactants
    i=0; % counter to make sure induction only happens once
    while t<out(end)
           alpha=feval(prop,state,r);
           alpha_o=sum(alpha); %total
           tau=-(1/alpha_o)*log(rand()); %determing time interval tau using random number
           i=sum(rand >= cumsum([0 alpha]/alpha_o));
           state = state + stm(:,i);
           times(j)=t+tau; %store times that results correspond to
           t_runs=times(1,:); %stores times for each run
           t=t+tau;
           while j<=step && t>out(j)
               results(:,j)=state;
               j = j+1;
           end
           t_runs=times(1,:); %stores times for each run
           [row col]=max(results(8,:)); %stores peak luciferase dimensions
           peak_luc(k)=row; %stores peak luciferase values
           peak_luc_time(k)=t_runs(col); %stores time of peak luciferase 
           end_luc(k)=results(8,step-1); %stores amount of luciferase after each simulation run
           cum_luc(k)=sum(results(8,:)); %stores cumulative luciferace output for each simulation run
           t=t+tau;
           if t>(1.8776*60) & t<((1.8776*60)+10) %small window after 1.8776-hour mark during which arabinose is added
               state(1) = 1000;
           end
    end
    k=k+1;
end

results(:,end)=state; %stores the results of the last run
t_runs(:,end)=end_time; %inputs the last time at the end_time

%mean, SD, CV of luc
format shortG
end_stats=[mean(end_luc) std(end_luc) std(end_luc)/mean(end_luc)]
peak_luc_stats=[mean(peak_luc) std(peak_luc) std(peak_luc)/mean(peak_luc)]
peak_time_stats= [mean(peak_luc_time) std(peak_luc_time) std(peak_luc_time)/mean(peak_luc_time)]./60
cum_luc_stats=[mean(cum_luc) std(cum_luc) std(cum_luc)/mean(cum_luc)]

results(:,end)=state; %stores the results of the last run
t_runs(:,end)=end_time; %inputs the last time at the end_time

figure
subplot(3,3,1)
plot(out./60,results(1,:))
xlabel('Time (hours)');
ylabel('Number arabinose molecules')

subplot(3,3,2)
plot(out./60,results(2,:))
xlabel('Time (hours)');
ylabel('Number luxA mRNA')

subplot(3,3,3)
plot(out./60,results(3,:))
xlabel('Time (hours)');
ylabel('Number luxB mRNA')

subplot(3,3,4)
plot(out./60,results(4,:))
xlabel('Time (hours)');
ylabel('Number luxCDE mRNA')

subplot(3,3,5)
plot(out./60,results(5,:))
xlabel('Time (hours)');
ylabel('Number luxA proteins')

subplot(3,3,6)
plot(out./60,results(6,:))
xlabel('Time (hours)');
ylabel('Number luxB proteins')

subplot(3,3,7)
plot(out./60,results(7,:))
xlabel('Time (hours)');
ylabel('Number luxCDE proteins')

subplot(3,3,8)
plot(out./60,results(8,:))
xlabel('Time (hours)');
ylabel('Number luciferase dimers')

subplot(3,3,9)
plot(out./60,results(9,:))
xlabel('Time (hours)');
ylabel('Number active genes')

%% AIM 3e: Split lux System opposite - shifted induction - luxA/B const, luxC/D/E ind
%induction 1.8776h in
clear all
close all

%initialize values, initial conditions of molecule states
mol = [0 0 0 0 0 0 0 0 1]; %active gene added as last value in mol matrix
ara=mol(1);
a_m=mol(2);
b_m=mol(3);
cde_m=mol(4);
a_p=mol(5);
b_p=mol(6);
cde_p=mol(7);
luc=mol(8);
active_gene=mol(9);

%defining reaction rate
rate = [3.054 0.1386 0.1386 0.1386 2.0016 2.0016 2.0016 0.000578 0.000578 0.000578 0.1695 .00385 .0693 3.054]; %rates all in molecules/minute
r.AB_mtran=rate(1); %luxAB transcription rate
r.A_mde=rate(2); %luxA mRNA degradation rate
r.B_mde=rate(3); %luxB mRNA degradation rate
r.CDE_mde=rate(4); %luxCDE mRNA degradation rate
r.A_ptran=rate(5); %luxA mRNA translation rate
r.B_ptran=rate(6); %luxB mRNA translation rate
r.CDE_ptran=rate(7); %luxCDE mRNA translation rate
r.A_pde=rate(8); %luxA protein degradation rate
r.B_pde=rate(9); %luxB protein degradation rate
r.CDE_pde=rate(10); %luxCDE protein degradation rate
r.AB_dim=rate(11); %luxA and luxB dimerization rate
r.AB_deg=rate(12); %luxAB dimer degradation/dissocation rate
r.ara_deg=rate(13); %arabinose degradation rate
r.CDE_mtran=rate(14); %luxCDE transcription rate

%propensity functions
prop=@(m,r)[r.AB_mtran*m(9) r.A_mde*m(2) r.B_mde*m(3) r.CDE_mde*m(4) r.A_ptran*m(2) r.B_ptran*m(3) r.CDE_ptran*m(4) r.A_pde*m(5) r.B_pde*m(6) r.CDE_pde*m(7) r.AB_dim*(m(5)*m(6)*m(7)) r.AB_deg*m(8) r.ara_deg*m(1) r.CDE_mtran*m(1)];

%stoich matrix
stm=[0	0	0	0	0	0	0	0	0	0	0	0	-1	0; 1	-1	0	0	0	0	0	0	0	0	0	0	0	0; 1	0	-1	0	0	0	0	0	0	0	0	0	0	0; 0	0	0	-1	0	0	0	0	0	0	0	0	0	1; 0	0	0	0	1	0	0	-1	0	0	-1	0	0	0; 0	0	0	0	0	1	0	0	-1	0	-1	0	0	0; 0	0	0	0	0	0	1	0	0	-1	0	0	0	0; 0	0	0	0	0	0	0	0	0	0	1	-1	0	0; 0	0	0	0	0	0	0	0	0	0	0	0	0	0];

runs=5; %put how many runs of simulation you want
step=50; %put how many timesteps you want each run to have
end_time=360; %how many seconds you want sim to run for
%720 minutes chosen because want to run simulation for 6 h

peak_luc=zeros(runs,1);
peak_luc_time=zeros(runs,1);
end_luc=zeros(runs,1);
cum_luc=zeros(runs,1);

t_runs=zeros(runs,step);

out=linspace(0,end_time,step); %output trajectory
t=0;

for k=1:length(end_luc)
    t=0; %initial time
    j=1; %index
    times=zeros(1,length(out)); %matrix to collect time snapshots
    state=mol(:); %system state
    results=zeros(length(mol),length(out)); %setup output matrix
    results(:,1)=mol(:); % initialize reactants
    i=0; % counter to make sure induction only happens once
    while t<out(end)
           alpha=feval(prop,state,r);
           alpha_o=sum(alpha); %total
           tau=-(1/alpha_o)*log(rand()); %determing time interval tau using random number
           i=sum(rand >= cumsum([0 alpha]/alpha_o));
           state = state + stm(:,i);
           times(j)=t+tau; %store times that results correspond to
           t_runs=times(1,:); %stores times for each run
           t=t+tau;
           while j<=step && t>out(j)
               results(:,j)=state;
               j = j+1;
           end
           t_runs=times(1,:); %stores times for each run
           [row col]=max(results(8,:)); %stores peak luciferase dimensions
           peak_luc(k)=row; %stores peak luciferase values
           peak_luc_time(k)=t_runs(col); %stores time of peak luciferase 
           end_luc(k)=results(8,step-1); %stores amount of luciferase after each simulation run
           cum_luc(k)=sum(results(8,:)); %stores cumulative luciferace output for each simulation run
           t=t+tau;
           if t>(1.8776*60) & t<((1.8776*60)+10)%small window after 1.8776-hour mark during which arabinose is added
               state(1) = 1000;
           end
    end
    k=k+1;
end

results(:,end)=state; %stores the results of the last run
t_runs(:,end)=end_time; %inputs the last time at the end_time

%mean, SD, CV of luc
format shortG
end_stats=[mean(end_luc) std(end_luc) std(end_luc)/mean(end_luc)]
peak_luc_stats=[mean(peak_luc) std(peak_luc) std(peak_luc)/mean(peak_luc)]
peak_time_stats= [mean(peak_luc_time) std(peak_luc_time) std(peak_luc_time)/mean(peak_luc_time)]./60
cum_luc_stats=[mean(cum_luc) std(cum_luc) std(cum_luc)/mean(cum_luc)]

results(:,end)=state; %stores the results of the last run
t_runs(:,end)=end_time; %inputs the last time at the end_time

figure
subplot(3,3,1)
plot(out./60,results(1,:))
xlabel('Time (hours)');
ylabel('Number arabinose molecules')

subplot(3,3,2)
plot(out./60,results(2,:))
xlabel('Time (hours)');
ylabel('Number luxA mRNA')

subplot(3,3,3)
plot(out./60,results(3,:))
xlabel('Time (hours)');
ylabel('Number luxB mRNA')

subplot(3,3,4)
plot(out./60,results(4,:))
xlabel('Time (hours)');
ylabel('Number luxCDE mRNA')

subplot(3,3,5)
plot(out./60,results(5,:))
xlabel('Time (hours)');
ylabel('Number luxA proteins')

subplot(3,3,6)
plot(out./60,results(6,:))
xlabel('Time (hours)');
ylabel('Number luxB proteins')

subplot(3,3,7)
plot(out./60,results(7,:))
xlabel('Time (hours)');
ylabel('Number luxCDE proteins')

subplot(3,3,8)
plot(out./60,results(8,:))
xlabel('Time (hours)');
ylabel('Number luciferase dimers')

subplot(3,3,9)
plot(out./60,results(9,:))
xlabel('Time (hours)');
ylabel('Number active genes')

%% AIM 4: Split lux System - shifted induction - changing ratio arabinose to active gene
%induction after 1.8776 hours
clear all
close all

%initialize values, initial conditions of molecule states
mol = [0 0 0 0 0 0 0 0 1000]; %active gene added as last value in mol matrix
%ratio adjusted by increasing number of active genes 10X each run
ara=mol(1);
a_m=mol(2);
b_m=mol(3);
cde_m=mol(4);
a_p=mol(5);
b_p=mol(6);
cde_p=mol(7);
luc=mol(8);
active_gene=mol(9);

%defining reaction rate
rate = [3.054 0.1386 0.1386 0.1386 2.0016 2.0016 2.0016 0.000578 0.000578 0.000578 0.1695 .00385 .0693 3.054]; %rates all in molecules/minute
r.AB_mtran=rate(1); %luxAB transcription rate
r.A_mde=rate(2); %luxA mRNA degradation rate
r.B_mde=rate(3); %luxB mRNA degradation rate
r.CDE_mde=rate(4); %luxCDE mRNA degradation rate
r.A_ptran=rate(5); %luxA mRNA translation rate
r.B_ptran=rate(6); %luxB mRNA translation rate
r.CDE_ptran=rate(7); %luxCDE mRNA translation rate
r.A_pde=rate(8); %luxA protein degradation rate
r.B_pde=rate(9); %luxB protein degradation rate
r.CDE_pde=rate(10); %luxCDE protein degradation rate
r.AB_dim=rate(11); %luxA and luxB dimerization rate
r.AB_deg=rate(12); %luxAB dimer degradation/dissocation rate
r.ara_deg=rate(13); %arabinose degradation rate
r.CDE_mtran=rate(14); %luxCDE transcription rate

%propensity functions
prop=@(m,r)[r.AB_mtran*m(1) r.A_mde*m(2) r.B_mde*m(3) r.CDE_mde*m(4) r.A_ptran*m(2) r.B_ptran*m(3) r.CDE_ptran*m(4) r.A_pde*m(5) r.B_pde*m(6) r.CDE_pde*m(7) r.AB_dim*(m(5)*m(6)*m(7)) r.AB_deg*m(8) r.ara_deg*m(1) r.CDE_mtran*m(9)];

%stoich matrix
stm=[0	0	0	0	0	0	0	0	0	0	0	0	-1	0; 1	-1	0	0	0	0	0	0	0	0	0	0	0	0; 1	0	-1	0	0	0	0	0	0	0	0	0	0	0; 0	0	0	-1	0	0	0	0	0	0	0	0	0	1; 0	0	0	0	1	0	0	-1	0	0	-1	0	0	0; 0	0	0	0	0	1	0	0	-1	0	-1	0	0	0; 0	0	0	0	0	0	1	0	0	-1	0	0	0	0; 0	0	0	0	0	0	0	0	0	0	1	-1	0	0; 0	0	0	0	0	0	0	0	0	0	0	0	0	0];

runs=3; %put how many runs of simulation you want
step=50; %put how many timesteps you want each run to have
end_time=360; %how many seconds you want sim to run for
%720 minutes chosen because want to run simulation for 6 h

peak_luc=zeros(runs,1);
peak_luc_time=zeros(runs,1);
end_luc=zeros(runs,1);
cum_luc=zeros(runs,1);

t_runs=zeros(runs,step);

out=linspace(0,end_time,step); %output trajectory
t=0;

for k=1:length(end_luc)
    t=0; %initial time
    j=1; %index
    times=zeros(1,length(out)); %matrix to collect time snapshots
    state=mol(:); %system state
    results=zeros(length(mol),length(out)); %setup output matrix
    results(:,1)=mol(:); % initialize reactants
    i=0; % counter to make sure induction only happens once
    while t<out(end)
           alpha=feval(prop,state,r);
           alpha_o=sum(alpha); %total
           tau=-(1/alpha_o)*log(rand()); %determing time interval tau using random number
           i=sum(rand >= cumsum([0 alpha]/alpha_o));
           state = state + stm(:,i);
           times(j)=t+tau; %store times that results correspond to
           t_runs=times(1,:); %stores times for each run
           t=t+tau;
           while j<=step && t>out(j)
               results(:,j)=state;
               j = j+1;
           end
           t_runs=times(1,:); %stores times for each run
           [row col]=max(results(8,:)); %stores peak luciferase dimensions
           peak_luc(k)=row; %stores peak luciferase values
           peak_luc_time(k)=t_runs(col); %stores time of peak luciferase 
           end_luc(k)=results(8,step-1); %stores amount of luciferase after each simulation run
           cum_luc(k)=sum(results(8,:)); %stores cumulative luciferace output for each simulation run
           t=t+tau;
           if t>(1.8776*60) & t<((1.8776*60)+10)
               state(1) = 1000;
           end
    end
    k=k+1;
end

results(:,end)=state; %stores the results of the last run
t_runs(:,end)=end_time; %inputs the last time at the end_time

%mean, SD, CV of luc
format shortG
end_stats=[mean(end_luc) std(end_luc) std(end_luc)/mean(end_luc)]
peak_luc_stats=[mean(peak_luc) std(peak_luc) std(peak_luc)/mean(peak_luc)]
peak_time_stats= [mean(peak_luc_time) std(peak_luc_time) std(peak_luc_time)/mean(peak_luc_time)]./60
cum_luc_stats=[mean(cum_luc) std(cum_luc) std(cum_luc)/mean(cum_luc)]

results(:,end)=state; %stores the results of the last run
t_runs(:,end)=end_time; %inputs the last time at the end_time

figure
subplot(3,3,1)
plot(out./60,results(1,:))
xlabel('Time (hours)');
ylabel('Number arabinose molecules')

subplot(3,3,2)
plot(out./60,results(2,:))
xlabel('Time (hours)');
ylabel('Number luxA mRNA')

subplot(3,3,3)
plot(out./60,results(3,:))
xlabel('Time (hours)');
ylabel('Number luxB mRNA')

subplot(3,3,4)
plot(out./60,results(4,:))
xlabel('Time (hours)');
ylabel('Number luxCDE mRNA')

subplot(3,3,5)
plot(out./60,results(5,:))
xlabel('Time (hours)');
ylabel('Number luxA proteins')

subplot(3,3,6)
plot(out./60,results(6,:))
xlabel('Time (hours)');
ylabel('Number luxB proteins')

subplot(3,3,7)
plot(out./60,results(7,:))
xlabel('Time (hours)');
ylabel('Number luxCDE proteins')

subplot(3,3,8)
plot(out./60,results(8,:))
xlabel('Time (hours)');
ylabel('Number luciferase dimers')

subplot(3,3,9)
plot(out./60,results(9,:))
xlabel('Time (hours)');
ylabel('Number active genes')

