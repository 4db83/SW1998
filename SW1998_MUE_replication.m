% This file exactly replicates the empirical part of SW(1998) with their data.
% See also the readme.file in the directory StockWatson\tvpci Stock_Watson_JASA_1998\
% and the corresonding PROCS lnair1.prc lnair2.prc which use a diffuse prior for initial
% Value and where the initial value is estimated as nuissance parameter.
% 	MLE_GDP1.gss:  MMLE Results (lnair1.prc --> diffuse prior for initial value p00 = vague=1e+6)
% 	MLE_GDP2.gss:  MPLE Results (lnair2.prc --> initial value is estimated)
% ---------------------------------------------------------------------------------------------
% MODEL IS 
% ---------------------------------------------------------------------------------------------
% 			GY(t) = beta(t) + u(t);
% 	D(beta(t) = beta(t-1) + (Lambda/T)*eta(t);
%    a(L)u(t) = epsilon(t);
% where a(L) = [1 -a1L -a2L^2 -a3L^3 -a4L^4] -> AR(4).
% ---------------------------------------------------------------------------------------------

clear; clc;
set(groot,'defaultLineLineWidth',2); % sets the default linewidth to 1.5;
% ADD LOCAL FUNCTION PATH WITHOUT SUBFOLDERS (IE. '_OLD') (only add what is needed)
addpath(genpath('./localFunctions/'));
% CALL: restoredefaultpath (to restore to default path)

% figdir = './UC.figures/';
data_dir			= './_data/';
SW98_data_dir = './_data/SW98.inputs/';
% if ~(exist(figdir)==7);mkdir(figdir);end; 
% if ~(exist(matdir)==7);mkdir(matdir);end; 
data_name = 'Gy_data.txt';

% GET THEIR DATA GY = 400*LOG(DIFF(REAL GDP PER CAPITA)) SERIES. 
GY = load([SW98_data_dir data_name]);
% make quarterly date vector from 1947:Q2 - 1995:Q4;
dates	= make_quarterly_dates(1947,2,1995,4);
% % convert dates to end of quarter -> not really needed.
% dates_eoq = make_end_of_quarter_dates(dates);

% PRINT DATES OF TIME PERIOD CONSIDERED (5:end) is 1961:Q1.
% Dates	= datenum(hlw_data.Time(5:end));
% sep('=')
fprintf( '					Sample period is: %s', datestr(dates(1),'yyyy:qq'));
fprintf( ' to %s \n', datestr(dates(end),'yyyy:qq'));
% sep('=')

%% 1) FIT AN AR(4) TO DEMEANED X, x_xbar = x-mean(x). Note X = GY(t), the growth rate series 
% X = X(11:end);
x_meanx = demean(GY);	
lags_GY	= mlag(GY,4);
% *********************************************************************************************
% THIS IS WHAT THEY DO IN THE CODE! SOME TYPO (SEE FILE TST_GDP1.GSS, lines 40 to 47)
% *********************************************************************************************
XX = [mlag(x_meanx,1) lags_GY(:,2:end)];
% using mlag(GY,4) correct series
% XX = lags_GY;
% ---------------------------------------------------------------------------------------------

AR4_out = ols(x_meanx,XX,1,{'AR(1)','AR(2)','AR(3)','AR(4)'});
% ALTERNATIVE IMPLEMENATION WITH CALL TO ARMA FUNCTION, but note that they have a mistake in
% the putting together of the XX variable in the AR(4) regression, so results are not excactly the same.
% AR4			= estimate_armax(GY,1,[1 2 3 4],[0]); 
% print_arma_results(AR4);

% a(L) lag ploynomial weights
aL			= [1 -AR4_out.bhat'];
% what they call the Long-run standard deviation: Long-Run variance is VAR(uhat)/a(1)^2;
LRstd_u = std(AR4_out.uhat)/sum(aL);

%% 2) INVERSE-FILTER X SERIES TO REMOVE MILD AR(4) CORRLEATION, -> USE aL IN FITLER COMMAND.
% using the first 4 observations as initial values --> exactly the same SW xfilt series.
xfilt_0 = filter(aL,1,GY,GY(1:4));
% drop the first 4 observations
xfilt	= xfilt_0(5:end);
% NOTE: this is exactly the same as AR(4) resids + mean(xfilt_0), also we do not need to add the mean
% back in; it does not matter for the testing, as e below is demeaned for Nyblom's test
% xfilt	= AR4_out.uhat + mean(xfilt_0);
% xfilt	= AR4_out.uhat ; % all give same results

%% 3) Do rolling 'Chow' Break point test, by partitioning the data for every time period. 
% Note: SW98 use threshold of 15% for minium sample sizes in the two samples that are tested.
CC	= .15;
T		= size(xfilt,1);
t1	= floor(T*CC);
t2	= T-t1;
Tvec= (1:T)';

% space for Fstat.
Fstat_0 = nan(T,1);
% now loop trought t0 to tT
for tt = t1:t2
 	ols_out			= fullols(xfilt,Tvec > tt);
	Fstat_0(tt) = ols_out.Fstat;  % ols_out.tstat(2).^2 SAME 
end

% Wald versions of the test
Fstat			 = removenans(Fstat_0);
sup_Fstat  = max(Fstat);
exp_Wald	 = log(mean(exp(Fstat/2)));
mean_Fstat = mean(Fstat);

% Nyblom's test is simply: e_cumsum'e_cumsum/T divided by Var(e), where e = xfilt-mean(xfilt).
e			= demean(xfilt);                     % e=y-meanc(y);               
e_cs	= cumsum(e)/sqrt(T);                 % es=cumsumc(e)/sqrt(rows(e));
Lstat	= (e_cs'*e_cs/T)/var(xfilt);         % seesq=e'e/(rows(e)-1);      
                                           % l=(es'es/rows(e))/seesq;    
% Combine tests for printing later
stats.L  = Lstat; 
stats.MW = mean_Fstat; 
stats.EW = exp_Wald; 
stats.SW = sup_Fstat;	 % SW = QLR in their notation

% ---------------------------------------------------------------------------------------------
% TO GET THE MU ESTIMATOR LAMBDA VALUES, READ IN THE LOOKUP TABLE FROM THE GAUSS .OUT FILE
% READ IN THE DATA FROM 
% ---------------------------------------------------------------------------------------------
lookup_table_name = [SW98_data_dir 'lookup_table_original_grid.out'];
% check if exists, otherwise load the .out file from GAUSS
if exist([SW98_data_dir 'SW1998_MUE_lookup_table.mat'],'file')
	load([SW98_data_dir 'SW1998_MUE_lookup_table.mat']);
else
	fid		= fopen(lookup_table_name); 
	% this defines the read format '%[^\n\r]' this reads the last line and carrgie return
	frmt	= [repmat('%q',1,6) '%[^\n\r]']; 
	% include 'headerlines', 0, if you don't want to read the header
	scanned_data = textscan(fid, frmt, 'delimiter', ' ', 'MultipleDelimsAsOne', true, 'CollectOutput',1);
	fclose(fid);
	% now make a table object
	for ii = 1:size(scanned_data{1},2)
		lookup_tmp(:,ii) = str2double(scanned_data{1}(2:end,ii));
	end
	% make a lookup table 
	var_names = scanned_data{1}(1,:);
	var_names{1} = 'Lambda';
	var_names{5} = 'SW';
	
	lookup_table = array2table(lookup_tmp, 'VariableNames', var_names);
	save([SW98_data_dir 'SW1998_MUE_lookup_table.mat'], 'lookup_table');
end
% ---------------------------------------------------------------------------------------------
% MAKE TEST CDF TABLE/STRUCTURE (201X199) FOR EACH TEST. NEEDED FOR THE P-VALUES AND CIS.
% ---------------------------------------------------------------------------------------------
if exist([SW98_data_dir 'SW1998_MUE_lookup_TSTCDF.mat'],'file')
	load([SW98_data_dir 'SW1998_MUE_lookup_TSTCDF.mat']);
else 
	TSTCDF.L  = xlsread([SW98_data_dir 'TSTCDF_Ltab.csv' ]);
	TSTCDF.MW = xlsread([SW98_data_dir 'TSTCDF_mwtab.csv']);
	TSTCDF.EW = xlsread([SW98_data_dir 'TSTCDF_ewtab.csv']);
	TSTCDF.SW = xlsread([SW98_data_dir 'TSTCDF_swtab.csv']);
	TSTCDF.pdvec = (0.005:0.005:0.995)';
	save([SW98_data_dir 'SW1998_MUE_lookup_TSTCDF.mat'], 'TSTCDF');
end
% ---------------------------------------------------------------------------------------------

%% COMPUTE MUE LAMBDA HAT FROM
Lambdas.L	 	= lookup_table_mue(stats.L , lookup_table.Lambda, lookup_table.L );
Lambdas.MW  = lookup_table_mue(stats.MW, lookup_table.Lambda, lookup_table.MW);
Lambdas.EW  = lookup_table_mue(stats.EW, lookup_table.Lambda, lookup_table.EW);
Lambdas.SW  = lookup_table_mue(stats.SW, lookup_table.Lambda, lookup_table.SW);

% COMPUTE THE IMPLIED STDEV(DELTA_BETA) AS:
% inline function Sigma_{\Detla \beta} = Lambda/T*LRStdev(eps_t), called sv=stdc(ee)/(1-sumc(beta[1:nlag]));" in SW1998
get_sigma_D_beta = @(lambda) ( lambda./T *LRstd_u );
sigma_D = structfun(get_sigma_D_beta, Lambdas ,'UniformOutput', false);

% COMPUTE THE P-VALUES
% inline function to find p-value
get_p_value = @(x) ( 1-TSTCDF.pdvec(x == min(x)) );
p_values.L	= get_p_value( abs((TSTCDF.L (1,:) - stats.L  )) );  
p_values.MW	= get_p_value( abs((TSTCDF.MW(1,:) - stats.MW )) );  
p_values.EW	= get_p_value( abs((TSTCDF.EW(1,:) - stats.EW )) );  
p_values.SW	= get_p_value( abs((TSTCDF.SW(1,:) - stats.SW )) );  

% COMPUTE THE LOWER CI BOUND ON LAMBDA 90% CI
CI 	= .90;
TAR = (1-CI)/2;
% find columns corresponding to this target 
[~,Iup ] = min(abs( TSTCDF.pdvec - TAR ));
[~,Ilow] = min(abs( TSTCDF.pdvec - 1-TAR ));

Lambda_up_CI.L	= lookup_table_mue(stats.L	, lookup_table.Lambda, TSTCDF.L	(:,Iup));
Lambda_up_CI.MW	= lookup_table_mue(stats.MW	, lookup_table.Lambda, TSTCDF.MW(:,Iup));
Lambda_up_CI.EW	= lookup_table_mue(stats.EW	, lookup_table.Lambda, TSTCDF.EW(:,Iup));
Lambda_up_CI.SW	= lookup_table_mue(stats.SW	, lookup_table.Lambda, TSTCDF.SW(:,Iup));

Lambda_low_CI.L		= lookup_table_mue(stats.L	, lookup_table.Lambda, TSTCDF.L	(:,Ilow));
Lambda_low_CI.MW	= lookup_table_mue(stats.MW	, lookup_table.Lambda, TSTCDF.MW(:,Ilow));
Lambda_low_CI.EW	= lookup_table_mue(stats.EW	, lookup_table.Lambda, TSTCDF.EW(:,Ilow));
Lambda_low_CI.SW	= lookup_table_mue(stats.SW	, lookup_table.Lambda, TSTCDF.SW(:,Ilow));

% make the upper bound for Sigma_Delta_beta
up_sigma_D	= structfun(get_sigma_D_beta, Lambda_up_CI ,'UniformOutput', false);
low_sigma_D	= structfun(get_sigma_D_beta, Lambda_low_CI,'UniformOutput', false);

% PRINT THE RESULTS TO SCREEN (THESE ARE EXACTLY THE SAME AS IN SW1998)
MUE_stats			= [stats; p_values; Lambdas; Lambda_low_CI; Lambda_up_CI; sigma_D; low_sigma_D; up_sigma_D];
MUE_rownames	= {'Tests', 'p-vals', 'Lambdas', 'Low CI', 'Up CI', 'Sigma_D', 'Low CI', 'Up CI'};
print2screen(	MUE_stats,[], ...
							MUE_rownames,[4 10], [] ,...
							100);					
%  							MUE_rownames,[4 10], '\_output\Table4.xls',...						

%  FROM SW98 Median Unbiased Estimators
% L-Robust 		4.0558657193 0.1303031273 
% MW-Robust		3.4335430139 0.1103097152 
% EW-Robust		3.0712033611 0.0986687998 
% SW-Robust		0.7786146171 0.0250146150 

%% 4) NOW SET UP KALMAN FILTER.
% ---------------------------------------------------------------------------------------------
% STATE SPACE MODEL:
% 		Observed:	y_t			= D_t + M*alpha_t			+ e_t;		Var(e_t) = H.
% 		State:		alpha_t = C_t + Phi*alpha_t-1	+ R*n_t;	Var(n_t) = Q.
% CALL AS: 
% 		[LogLik, att, Ptt] = kalmanfilter(y, D, M, H, C, Phi, Q, R, a1, P1) 
% Pmean.att = kalmanfilter(Y, Pmean.Dt, Pmean.M, Pmean.H, Ct, Phi, Pmean.Q, R, a00, P00); 
% ---------------------------------------------------------------------------------------------
% SOME NUMERICAL OPTIMISATON SETTINGS (fair robust to different settings)
% ---------------------------------------------------------------------------------------------
numOptions = optimoptions(@fminunc, 'Display','off', 'Algorithm','quasi-newton', ...
	'FiniteDifferenceType','central', 'HessUpdate','bfgs', ...
	'TolFun',1e-08, 'TolX',1e-08, 'MaxFunEvals',5e5, 'MaxIter',5e5);

% ---------------------------------------------------------------------------------------------
% DATA TO BE PASSED INTO LOGLIKELIHOOD FUNCTION
data_in = GY'; % transpose to make (1xT) vector for KF routine

% *********************************************************************************************
% FROM RUNNING THE SW GAUSS FILES IN FOLDER D ... ecb\StockWats\tvpci
% Stock_Watson_JASA_1998\
% NOTE THAT I HAVE MODIFIED THESE TO PRINT WITH MORE ACCURACY AND HAVE INCREASE THE MAXLIK
% TOLERENACE FOR MORE OPTIMISTION PRECISION format /rd 14,14;
% @ _max_Diagnostic = 2; @
% _max_GradTol    = 1e-08;
% _max_MaxTime    = 1;
% _max_MaxIters   = 300;
% *********************************************************************************************
% THESE ARE THE G013 MLE ESTIMATES AS WELL AS THE INITIAL VALUES OF THE STATE VECTOR FROM SW1998
% MLE_GDPC.GSS CALLING LNAIRC.PRC
% ---------------------------------------------------------------------------------------------
SW_G13_bhat  = [ 2.44099940415309        % b00                
							 	 3.84661916899307        % sigma_eps          
							 	 0.33501453946479        % AR(1)              
							 	 0.12742309144336        % AR(2)              
							 	-0.01017052375416        % AR(3)              
							 	-0.08680297677321        % AR(4)              
							 	 0.13000000000000 ] ;    % sigma_Delta\beta   

% ---------------------------------------------------------------------------------------------
% MLE_GDP1.GSS CALLING LNAIR1.PRC USING A DIFFUSE_PRIOR. p00 = vague=1e+6
% ---------------------------------------------------------------------------------------------
SW_MMLE_bhat = [ 3.85859383842156        % sigma_eps        
								 0.34025236220875        % AR(1)            
								 0.13074613643110        % AR(2)            
								-0.00725117007105        % AR(3)            
								-0.08247858673446        % AR(4)            
								 0.04440086128563 ];     % sigma_Delta\beta 
				 
% ---------------------------------------------------------------------------------------------
% MLE_GDP2.GSS CALLING LNAIR2.PRC ESTIMATING THE INTIAL CONDITION.
% ---------------------------------------------------------------------------------------------
SW_MPLE_bhat = [ 1.79589935216043        % b00
							   3.85199477183636        % sigma_eps
							   0.33708320734919        % AR(1)
							   0.12890328120825        % AR(2)
							  -0.00917382000830        % AR(3)
							  -0.08564442026999        % AR(4)
							   0.00000000364121 ] ;    % sigma_Delta\beta
% *********************************************************************************************
		 
% re-arrange the parameters for later use and printing.	 
SWG13_fittedVals	= [SW_G13_bhat(3:6);SW_G13_bhat(2);SW_G13_bhat(1)];
SWG13_otherpars		= SW_G13_bhat(end);

% EVALUATE THE LIKELIHOOD (note SW compute the mean of the Log-Likelihood)
[~,SWG13]		= SW1998_LL_wrapper_G13_AR4(SWG13_fittedVals, data_in, SWG13_otherpars);
% LOAD THE SAVED SMOOTHED AND FILTERED OUTPUT FROM SW FILES
sf_g013			= load([SW98_data_dir  'smoothed_filtered_G013.txt']);
% COLLECT PARAMETERS: ORIGINAL SW ESTIMATES
SW_G13_pars	= [SWG13.LL; SW_G13_bhat(3:6); SW_G13_bhat(2); SW_G13_bhat(end); SW_G13_bhat(1)];
% also load the smothed and filtered MMLE estimats for comparison of the volatility of the filtered output.
% sf_mmle			= xlsread([SW98_data_dir  'mmle_smoothed_filtered.xls']);
sf_mmle				= load([SW98_data_dir  'mmle_smoothed_filtered.txt']);

% ---------------------------------------------------------------------------------------------
% db G13 REPLICATION
% ---------------------------------------------------------------------------------------------
% intial values for optimisation 
%									[a(L) AR(4) Sigma_e					 a00]
db_G13_initvals	= [-aL(2:end) sqrt(AR4_out.sig2) mean(GY)] ;
% sigma_D_beta		= 0.13;
sigma_D_beta		= sigma_D.L;
Lambda_T				= Lambdas.L/T;
% otherParams	= sigma_D_betas_(3);
% CALL THE OPTIMIZER ROUTINE
[db_G13_bhat,~,~,~,g_G13,H_G13]	= fminunc(@SW1998_LL_wrapper_G13_AR4, db_G13_initvals, numOptions, data_in, sigma_D_beta);
[~, db_G13]	= SW1998_LL_wrapper_G13_AR4(db_G13_bhat, data_in, sigma_D_beta, 1);
% COLLECT PARAMETERS: MY ESTIMATES G13 RESTRICTION
db_G13_pars	= [db_G13.LL; db_G13_bhat(1:5)'; sigma_D_beta; db_G13_bhat(6)];

% ***************************************************************************************************
% ***************** IF USING LAMBDA/T AS INPUT (results are NOT the same). **************************
% ***************************************************************************************************
% [db_G13_bhat,~,~,~,g_G13,H_G13]	= fminunc(@SW1998_LL_wrapper_G13_AR4, db_G13_initvals, numOptions, data_in, Lambda_T);
% [~, db_G13]	= SW1998_LL_wrapper_G13_AR4(db_G13_bhat, data_in, Lambda_T, 1);
% % COLLECT PARAMETERS: MY ESTIMATES G13 RESTRICTION
% db_G13_pars	= [db_G13.LL; db_G13_bhat(1:5)'; Lambda_T*db_G13_bhat(5)/sum(db_G13_bhat(1:4)); db_G13_bhat(6)];
% % % make sure to change the Q(1,1) definition in SW1998_LL_wrapper_G13_AR4
% ***************************************************************************************************

% NOW GET FILTERS SMOOTHED STATES FOR MYG13
db_G13_KFS	= kalman_filter_smoother(data_in, db_G13);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% comparions to GAUSS OUTPUT from \ecb\StockWatson\tvpci Stock_Watson_JASA_1998\MLE_GDPC.GSS;
% clc % no need to minus the Hessian because it is coming from minimization problem
% G13_ = [db_G13_bhat' g_G13, diag(sqrt(inv(H_G13)))];
% G13_([6;5;1;2;3;4],:)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---------------------------------------------------------------------------------------------
% db G62 REPLICATION (62 MEANS UPPER CI LIMIT USED FOR LAMBDA AND HENCE SIGMA_\Delta\beta)
% ---------------------------------------------------------------------------------------------
% intial values for optimisation 
%									[a(L) AR(4) Sigma_e					 a00]
db_G62_initvals	= [-aL(2:end) sqrt(AR4_out.sig2) mean(GY)] ;
sigma_D_beta		= 0.62;
% CALL THE OPTIMIZER ROUTINE
db_G62_bhat	= fminunc(@SW1998_LL_wrapper_G13_AR4, db_G62_initvals, numOptions, data_in, sigma_D_beta);
[~, db_G62]	= SW1998_LL_wrapper_G13_AR4(db_G62_bhat, data_in, sigma_D_beta, 1);
% COLLECT PARAMETERS: MY ESTIMATES G62 RESTRICTION
db_G62_pars	= [db_G62.LL; db_G62_bhat(1:5)'; sigma_D_beta; db_G62_bhat(6)];

% ---------------------------------------------------------------------------------------------
% db MPLE (MLE2) REPLICATION
% ---------------------------------------------------------------------------------------------
% intial values for optimisation 
db_MPLE_initvals		= [-aL(2:end) sqrt(AR4_out.sig2) 0.5 mean(GY)] ;
db_MPLE_otherpars		= [];
% CALL THE OPTIMIZER ROUTINE
db_MPLE_bhat	= fminunc(@SW1998_LL_wrapper_MPLE_AR4, db_MPLE_initvals, numOptions, data_in, db_MPLE_otherpars);
[~, db_MPLE]	= SW1998_LL_wrapper_MPLE_AR4(db_MPLE_bhat, data_in, db_MPLE_otherpars, 1);
% COLLECT PARAMETERS: MY ESTIMATES G62 RESTRICTION
db_MPLE_pars	= [db_MPLE.LL; db_MPLE_bhat(1:7)'];

% ---------------------------------------------------------------------------------------------
% db MMLE (MLE1) REPLICATION
% ---------------------------------------------------------------------------------------------
% intial values for optimisation 
db_MMLE_initvals	= [-aL(2:end) sqrt(AR4_out.sig2) 1] ;
% FROM MEDIAN UNBIASED ESTIMATION
diffuse_prior.a00 = 0;
diffuse_prior.P00 = 1e6;
diffuse_prior.sigma_z = 0.1;
% CALL THE OPTIMIZER ROUTINE (ALLS SAVE 
[db_MMLE_bhat,~,~,~,g_MMLE,H_MMLE] = fminunc(@SW1998_LL_wrapper_MMLE_AR4, db_MMLE_initvals, numOptions, data_in, diffuse_prior);
[~, db_MMLE]	= SW1998_LL_wrapper_MMLE_AR4(db_MMLE_bhat, data_in, diffuse_prior, 1);
% COLLECT PARAMETERS: MY ESTIMATES G62 RESTRICTION
db_MMLE_pars		= [db_MMLE.LL; db_MMLE_bhat(1:6)'; NaN];
db_MMLE_stderr	= [NaN; sqrt(diag(pinv(H_MMLE))); NaN];

% ---------------------------------------------------------------------------------------------
% db MMLE (MLE1) REPLICATION WITH SIGMA_Z = 0
% ---------------------------------------------------------------------------------------------
[db_MMLE_sigma_z_0_bhat]= fminunc(@SW1998_LL_wrapper_MMLE_sigma_z_0_AR4, db_MMLE_bhat(1:5), numOptions, data_in, diffuse_prior);
[~, db_MMLE_sigma_z_0]	= SW1998_LL_wrapper_MMLE_sigma_z_0_AR4(db_MMLE_sigma_z_0_bhat, data_in, diffuse_prior, 1);
% COLLECT PARAMETERS: MY ESTIMATES G62 RESTRICTION
db_MMLE_sigma_z_0_pars	= [db_MMLE_sigma_z_0.LL; db_MMLE_sigma_z_0_bhat'; 0; NaN];

% ---------------------------------------------------------------------------------------------
% db MUE MW REPLICATION WITH LAMBDA / SIGMA_\DELTA\BETA based on ( mean Wald Statistic)
% ---------------------------------------------------------------------------------------------
% intial values for optimisation 
db_MUE_MW_initvals	= [-aL(2:end) sqrt(AR4_out.sig2) mean(GY)] ;
sigma_D_beta		= sigma_D.MW;
% CALL THE OPTIMIZER ROUTINE
db_MUE_MW_bhat	= fminunc(@SW1998_LL_wrapper_G13_AR4, db_MUE_MW_initvals, numOptions, data_in, sigma_D_beta);
[~, db_MUE_MW]	= SW1998_LL_wrapper_G13_AR4(db_MUE_MW_bhat, data_in, sigma_D_beta, 1);
% COLLECT PARAMETERS: MY ESTIMATES G62 RESTRICTION
db_MUE_MW_pars	= [db_MUE_MW.LL; db_MUE_MW_bhat(1:5)'; sigma_D_beta; db_MUE_MW_bhat(6)];

% ---------------------------------------------------------------------------------------------
% db MUE EW REPLICATION WITH LAMBDA / SIGMA_\DELTA\BETA based on (Expnential Wald Statistic)
% ---------------------------------------------------------------------------------------------
% intial values for optimisation 
db_MUE_EW_initvals	= [-aL(2:end) sqrt(AR4_out.sig2) mean(GY)] ;
sigma_D_beta		= sigma_D.EW;
% CALL THE OPTIMIZER ROUTINE
db_MUE_EW_bhat	= fminunc(@SW1998_LL_wrapper_G13_AR4, db_MUE_EW_initvals, numOptions, data_in, sigma_D_beta);
[~, db_MUE_EW]	= SW1998_LL_wrapper_G13_AR4(db_MUE_EW_bhat, data_in, sigma_D_beta, 1);
% COLLECT PARAMETERS: MY ESTIMATES G62 RESTRICTION
db_MUE_EW_pars	= [db_MUE_EW.LL; db_MUE_EW_bhat(1:5)'; sigma_D_beta; db_MUE_EW_bhat(6)];

% ---------------------------------------------------------------------------------------------
% db MUE QLR/SW REPLICATION WITH LAMBDA / SIGMA_\DELTA\BETA based on (Quandt LR statistic)
% ---------------------------------------------------------------------------------------------
% intial values for optimisation 
db_MUE_SW_initvals	= [-aL(2:end) sqrt(AR4_out.sig2) mean(GY)] ;
sigma_D_beta		= sigma_D.SW;
% CALL THE OPTIMIZER ROUTINE
db_MUE_SW_bhat	= fminunc(@SW1998_LL_wrapper_G13_AR4, db_MUE_SW_initvals, numOptions, data_in, sigma_D_beta);
[~, db_MUE_SW]	= SW1998_LL_wrapper_G13_AR4(db_MUE_SW_bhat, data_in, sigma_D_beta, 1);
% COLLECT PARAMETERS: MY ESTIMATES G62 RESTRICTION
db_MUE_SW_pars	= [db_MUE_SW.LL; db_MUE_SW_bhat(1:5)'; sigma_D_beta; db_MUE_SW_bhat(6)];

% % UNCOMMENT TO GET FILTERED/SMOOHTED STATES AS A CHECK
% SSM_MMLE = kfs2ssm(db_MMLE, data_in); % now compar to db_MMLE.KFS output
% CIint = 1.96*sqrt( squeeze(SSM_MMLE.PtT(1,1,:)) );	% 95%

%% PRINT RESULTS TO SCREEN
% ---------------------------------------------------------------------------------------------
% clc;
combined_pars = [	db_MPLE_pars	, ...
									db_MMLE_pars	, ...
									db_MMLE_stderr, ...
									db_MMLE_sigma_z_0_pars	, ...
									db_G13_pars		, ...
									db_G62_pars		, ...
									SW_G13_pars		, ...
% 									db_MUE_MW_pars, ...
% 									db_MUE_EW_pars, ...
% 									db_MUE_SW_pars, ...
								];
								
rowNames	= { 'Log-Like'				, ...
							'AR(1)'           , ...
							'AR(2)'           , ...
							'AR(3)'           , ...
							'AR(4)'           , ...
							'sigma_eps'       , ...
							'sigma_Delta_beta', ...
							'beta00'};

order4printing = [7;6;2;3;4;5;8;1];

modelNames = {	'MPLE' 				, ...
								'MMLE' 				, ...
								'stderr'			, ...
								'MMLE(s_z=0)' , ...
								'G13' 				, ...
								'G62' 				, ...
								'SW.File'			, ...
% 								'MUE.MW'			, ...
% 								'MUE.EW'			, ...
% 								'MUE.QLR'			, ...
							};

sep(128)
	print2screen(	combined_pars(order4printing,:) ,... 
								rowNames(order4printing), ...
								modelNames,'%14.8f');
% modelNames,'%14.8f','\_output\Table5.xls');							
sep(128)

% HP Filter results (not used at the moment)
[cycl_x,trnd_x,hp_wghts_x] = hp_filter(GY,(LRstd_u/sigma_D.EW)^2);

% ---------------------------------------------------------------------------------------------
% COMPUTE THE CONFIDENCE INTERVALS (CI) FOR THE EXTRACTED STATE VECTOR
% ---------------------------------------------------------------------------------------------
% CIint = 1.96*sqrt( squeeze(db_MMLE.KFS.PtT(1,1,:)) );		% 95%
CIint = 1.645*sqrt( squeeze(db_MMLE.KFS.PtT(1,1,:)) );		% 90%

% ---------------------------------------------------------------------------------------------
% LOAD THE NBER RECESSION INDICATORS. 
% ---------------------------------------------------------------------------------------------
load([data_dir 'us_gdp_NBER_recessions_data_2019Q4.mat']);
% trim the recession dates to size
recI = US_DATA(daterange_q(dates(1), dates(end)), 3);
% recession bar color
rec_CLR = .9*ones(3,1); 

%% PLOT THE RESULTS
% ---------------------------------------------------------------------------------------------
clf; % sets the default linewidth to 1.5;
set(groot,'defaultLineLineWidth',1.75); 
% plotting controls for font size, date frequncy, subplot dimension, and subplot space
Fns = 17;
Dfq = 19;
figDims = [.86 .22];
offSet  = .26;
fanCLR	= [.88 .77 .99];

subplot(2,1,1)
	% RECESSION BARS IN THE BACKGROUND
	bar( recI.NBER_rec_I*100,1,'FaceColor',rec_CLR,'EdgeColor',rec_CLR,'ShowBaseLine','off'); hold on;
	bar(-recI.NBER_rec_I*100,1,'FaceColor',rec_CLR,'EdgeColor',rec_CLR,'ShowBaseLine','on' );

% uncertainty bands for the Smoothed state vector from P.tT(1,1)
	fanplot([db_MMLE.KFS.atT(:,1)+CIint  db_MMLE.KFS.atT(:,1)-CIint],[], fanCLR); 

	plot( db_G13.KFS.atT(:,1)		,'Color',clr('b'));
	plot( db_G62.KFS.atT(:,1)		,'Color',clr('g'));
	plot(db_MPLE.KFS.atT(:,1)		,'Color',clr('r'));
	plot(db_MMLE.KFS.atT(:,1)		,'Color',clr('p'));
 	plot(sf_g013(:,1),'--'			,'Color',clr('o'));
	% THESE ARE THE KF FILTERD SERIES IN TOP
%  	plot(db_MMLE.KFS.att(:,1)		,'Color',.3*ones(1,3));
%  	plot(sf_mmle(2:end,2),'--'	,'Color',.8*ones(1,3));
	% add any of the other structural break test series here for comparison
% 	plot(db_MUE_EW.KFS.atT(:,1)	,'Color',clr('c'));
hold off;
% get top subplot location
SL = gca;
box on; grid on;
setplot([ SL.Position(2) figDims])
setyticklabels(-0:1:5,0,Fns)
setdateticks(dates,Dfq ,'yyyy:qq', Fns);
set(gca,'GridLineStyle',':','GridAlpha',1/3)
% NOW MAKE OUTSIDE TICKS FOR TWO AXIS Y LABELS 
setoutsideTicks
add2yaxislabel

subplot(2,1,2)
	% RECESSION BARS IN THE BACKGROUND
	bar( recI.NBER_rec_I*100,1,'FaceColor',rec_CLR,'EdgeColor',rec_CLR,'ShowBaseLine','off'); hold on;
	bar(-recI.NBER_rec_I*100,1,'FaceColor',rec_CLR,'EdgeColor',rec_CLR,'ShowBaseLine','on' );

% uncertainty bands for the Smoothed state vector from P.tT(1,1)
	fanplot([db_MMLE.KFS.atT(:,1)+CIint  db_MMLE.KFS.atT(:,1)-CIint], [], fanCLR); 

	LG(1) = plot(GY,'Color',.0*ones(1,3),'Marker','.','MarkerSize',15,'LineStyle',':','LineWidth',3/2);	
	LG(2) =	plot( db_G13.KFS.atT(:,1)	,'Color',clr('b'));
	LG(3) =	plot( db_G62.KFS.atT(:,1)	,'Color',clr('g'));
	LG(4) =	plot(db_MPLE.KFS.atT(:,1)	,'Color',clr('r'));
	LG(5) =	plot(db_MMLE.KFS.atT(:,1)	,'Color',clr('p'));
% 					plot(db_MMLE.KFS.att(:,1)		,'Color',clr('p'));
% 					plot(sf_mmle(2:end,2),'--'	,'Color',clr('c'));
	LG(6) =	plot(sf_g013(:,1),'--'		,'Color',clr('o'));
	
hold off;

hline(0)
box on; grid on;
setplot([SL.Position(2)-offSet figDims])
setyticklabels(-12:4:16,0,Fns)
setdateticks(dates,Dfq ,'yyyy:qq', Fns);
set(gca,'GridLineStyle',':','GridAlpha',1/3)
% NOW MAKE OUTSIDE TICKS FOR TWO AXIS Y LABELS 
setoutsideTicks; add2yaxislabel;

% ADD LEGEND
legendflex(LG, {'GY';
								'$\hat\sigma_{\Delta\beta}=0.13$';
								'$\hat\sigma^{\mathrm{up}}_{\Delta\beta}=0.62$';
								'MPLE';
								'MMLE';
								'SW.GAUSS';}, 'fontsize',Fns-4, 'Interpreter','Latex', 'buffer', [0 0]);


% --------------------------------------------------------------------------------------------
% UNCOMMENT TO PRINT TO FILE.
% print2pdf('SW98trend_growth_with_X','../write.up/graphics');
% print2pdf('SW98trend_growth','../write.up/graphics');
% print2pdf('SW98trend_growth_all','../write.up/graphics');
print2pdf('SW98trend_growth_all');
						
						
% sqrt(diag(inv(H_MMLE)))
						
						
						
						
						
						
						
						
%% EOF


%% OUTUPT FROM \ecb\StockWatson\tvpci Stock_Watson_JASA_1998\MLE_GDPC.GSS;
% 
% % 
% % ===============================================================================
% %  MAXLIK Version 5.0.6                                      6/23/2020  11:07 am
% % ===============================================================================
% %                                Data Set:  dummy                                
% % -------------------------------------------------------------------------------
% % 
% % 
% % return code =    0
% % normal convergence
% % 
% % Mean log-likelihood        -1.85384
% % Number of cases     195
% % 
% % Covariance matrix of the parameters computed by the following method:
% % Inverse of computed Hessian
% % 
% % Parameters    Estimates     Std. err.  Est./s.e.  Prob.    Gradient
% % ------------------------------------------------------------------
% % P01              2.4410        0.8959    2.725   0.0064      0.0000
% % P02              3.8466        0.1961   19.615   0.0000      0.0000
% % P03              0.3350        0.0719    4.658   0.0000      0.0000
% % P04              0.1274        0.0755    1.688   0.0914      0.0000
% % P05             -0.0102        0.0753   -0.135   0.8926      0.0000
% % P06             -0.0868        0.0714   -1.215   0.2244      0.0000
% % 
% % Correlation matrix of the parameters
% %    1.000  -0.035  -0.049  -0.033  -0.030  -0.049
% %   -0.035   1.000   0.014   0.007   0.009   0.016
% %   -0.049   0.014   1.000  -0.306  -0.097   0.057
% %   -0.033   0.007  -0.306   1.000  -0.259  -0.100
% %   -0.030   0.009  -0.097  -0.259   1.000  -0.309
% %   -0.049   0.016   0.057  -0.100  -0.309   1.000
% % 
% % Number of iterations    38
% % Minutes to convergence     0.00987
% % Std.error omputed with sqrt(diag(h)), where h must already be the inverse of Hessian












