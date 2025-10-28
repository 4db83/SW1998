function [sim] = simulate_from_full_Stage3(Theta3, Lambda_g, Lambda_z, fitRR, Z)
% function [data_in_with_intercept, YY, XX] = simulate_from_full_Stage3(Theta3, Lambda_g, Lambda_z, fitRR, Z)
% load the hlw_data, if we want to compare it to real data and then also compute HP fitlered trend HP_trend)
% ---------------------------------------------------------------------------------------------------
% Simulate from the Stage 3 model, with or without z(t) process
% NOTE: having a unit-root process for inflation is crap. In the simulations, it just drifts of
% way too much either up or down. Empri
% ---------------------------------------------------------------------------------------------------
% if == 1 plot the simulated series
plot_ = 0; 
% Lambda_g = 0.053869037947199;
% Lambda_z = 0.0302172212554928;

% THETA_DB FROM	STAGE 2 FOR NOW
ay1 = Theta3(1);                        %  1.529572488611140
ay2 = Theta3(2);                        % -0.587564151816759
ar	= Theta3(3);                        % -0.071195686242373
bpi = Theta3(4);                        %  0.668207053291876
by	= Theta3(5);                        %  0.078957783248161
sigma_ycycl = Theta3(6);                %  0.353468454246346  STDEV NOT VAR! 
sigma_infl  = Theta3(7);                %  0.789194865945594  STDEV NOT VAR! 
sigma_ytrnd = Theta3(8);                %  0.572419245803987  STDEV NOT VAR! 
sigma_g = Lambda_g*sigma_ytrnd;	        %  0.030835674073922  STDEV NOT VAR! 
sigma_z = abs(Lambda_z*sigma_ycycl/ar); %  0.150020807842209  STDEV NOT VAR! 

TT = 225; % sample size in HLW
Ts = TT + 4;	
% REAL RATE: SIMULATE ARMA PROCESS FOR RR (correct)
rr = filter( fitRR.bL, fitRR.aL, fitRR.c / sum(fitRR.bL) + randn(Ts, 1), fitRR.mu*ones(fitRR.max_pq,1) );

% SPACE ALLOCATION
g = nan(Ts,1);
z = nan(Ts,1);

% r = demean(RR);
rstar	= nan(Ts,1);
ycycl	= nan(Ts,1);
ytrnd	= nan(Ts,1);
infl	= nan(Ts,1);

% SHOCKS
e = randn(Ts,5);

% INITIALISATION OF SOME VARIABLES: THIS IS IMPORTANT
% inflation is initialised at observed/empirical inflation data from 1960:Q1:1960:Q4
infl(1:4)		= [ 1.280768901620000;
                1.538793385097180;
                1.601692078959750;
                1.337329153375280;];
						
% output trend initialized at HP_trend estimate
ytrnd(1:4)	= [ 802.9645190625;
                804.1247787453;
                805.2850855620;
                806.4454909469;];

ycycl(1:4)	= 0;										% start the cycle at 0
g(1:4)			= .75;									% start g(t) at 3% trend growth when annualized. I do not use diff(HP_trend) here as it is 1.16, which is much too high
% g(1:4)			= rr(1:4)/4;					% start g(t) at first few Real rate observations, to make sure that g(t) starts somewhere close to r(t), rather than way above. 
z(1:4)			= 0;										% start z(t) at zero
rstar(1:4)	= 4*g(1:4) + z(1:4);
% rr(1:4) = hlw_data.real_rate(1:4);
% ir(1:4)	= hlw_data.int_rate(1:4);

% STAGE 3 MODEL NO Z(T) 
for t = 5:Ts
  g(t) = g(t-1) + sigma_g*e(t,1);		% these are kept at quarterly frequency exactly as in paper
  z(t) = z(t-1) + sigma_z*e(t,2);
	
	% r*(t) without z(t) process
	if Z == 0 
		rstar(t)	= 4*g(t);								% now annualized exactly as in paper
	elseif Z == 1
		% r*(t) with z(t) process in it
		rstar(t)	= 4*g(t) + z(t);			% now annualized exactly as in paper
	end

% output gap/cycle 
  ycycl(t)	= ay1*ycycl(t-1) + ay2*ycycl(t-2) ...
						+	ar*( rr(t-1)-rstar(t-1) + rr(t-2)-rstar(t-2) )/2	+ sigma_ycycl*e(t,3);
  % precompute inflation(t-2):(t-4) as in HLW
	infl_24		= ( infl(t-2) + infl(t-3) + infl(t-4) )/3;
	%	inflation
  infl(t)		= bpi*infl(t-1) + (1-bpi)*infl_24 + by*ycycl(t-1) + sigma_infl*e(t,4);
	% trend 
	ytrnd(t)	= ytrnd(t-1) + g(t-1) + sigma_ytrnd*e(t,5);
	%	real rate in model with only nominal interest rate exogenous
% 	rr(t)			= ir(t) - ( infl(t) + infl(t-1) +  infl(t-2) + infl(t-3) )/4;
end

% make 100*log(GDP) of the simulated data by combining Cycle and Trend.
ySim		= ycycl + ytrnd;
piSim		= infl;

% analogue to YY = [lnGDP INFL] on simulated data 
sim.YY = [ySim piSim];             % YY	  = [lnGDP INFL]; 

% analogue to XX = [ GPD_2_lags RR_2_lags INFL_4_lags(:,1) mean(INFL_4_lags(:,2:4),2) ]; on simulated data 
%								   [GDP(t-1) GDP(t-2) r(t-1) r(t-2) pi(t-1) mean(pi(t-2:t-4))]
INFL_4_lags = mlag(piSim,4);
sim.XX = [mlag(ySim,2) mlag(rr,2) lag(piSim,1) mean(INFL_4_lags(:,2:4),2) ];

% now return data with intercept from the simulated model
data_in_with_intercept = [sim.YY(5:end,:) sim.XX(5:end,:) ones(TT,1)];			% WITH INTERCEPT

sim.with_intercept	= data_in_with_intercept';
sim.no_intercept		= data_in_with_intercept(:,1:end-1)';

% now add the full set of data to it
sim.g = g;
sim.z = z;
sim.cycle	= ycycl;
sim.trend	= ytrnd;
sim.gdp		= ycycl + ytrnd;
sim.rstar	=	rstar;
sim.realrate	= rr;
sim.inflation = infl;


% ---------------------------------------------------------------------------------------------
%% if needed, Plotting 
% ---------------------------------------------------------------------------------------------
if plot_ == 1
% plotting parameters
clf; ps = -1.175;
% Load the real data file to also have the date vector and comparions to observed data.
data_dir	= './_data/';
load([data_dir 'hlw_data_2017Q1'])
% run HP Filter on lnGDP to get initial crude trend and cycle decomposition
[HP_cycle] = hp_filter(100*hlw_data.gdp, 36000); % Lambda = 36000 is in HLW paper

subplot(4,1,1);
	plot(100*hlw_data.gdp); hold on;
	plot(sim.gdp) % simulated series
setdateticks(hlw_data.Time)
setyticklabels(800:50:1050,0,11)
subtitle('Output (blue data, red simulated)',ps)

subplot(4,1,2); 
	plot(hlw_data.infl); hold on;
	plot(sim.inflation) % simulated series
setdateticks(hlw_data.Time)
setyticklabels(-10:5:20,0,11)
hline(0)
subtitle('Inflation (blue data, red simulated)',ps)

subplot(4,1,3); 
	plot(hlw_data.real_rate); hold on;
	plot(sim.realrate); % simulated series
	lg(1) = plot(rstar);
	lg(2) =	plot(z);
	lg(3) =	plot(4*g);
	setdateticks(hlw_data.Time)
	setyticklabels(-8:4:16,0,11)
hline(0)
legend(lg,{'r*','z','4*g'},'Location','NorthWest')
subtitle('Real rate (blue data, red simulated), (4*g(t) and r*',ps)

subplot(4,1,4); 
	plot(HP_cycle); hold on;
	plot(sim.cycle)	% simulated series
	setdateticks(hlw_data.Time)
	setyticklabels(-10:5:10,0,11)
hline(0)
subtitle('Cycle or Output gap (blue data, red simulated)',ps)
end







%EOF
