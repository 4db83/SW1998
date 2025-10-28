%% Median Unbiased Estimation as implemented in HLW(2017) STAGE 1. 
function [Lambdas_ols, Lambdas_chow, Stats_ols, Stats_chow, ExtraOutput] = MUE_HLW_Stage1(GY0)
 	% REMOVE ANY NANS
	GY = removenan(GY0);

	T			= size(GY,1);
	Tvec	= (1:T)';
	% lower and upper time thresholds over which to search for structural break
	t0 = 4;	
	t1 = T-4;


	
	% NOW DE-MEAN YY AND COMPUTE SSE AS IN THE MUE STAGE 1 MODEL
	e0	= demean(GY); 
	SS0 = e0'*e0;

	% space for SS1 and t-stat
	% Fstat_0 = nan(T,1); 
	tstat_t = nan(T,2);   beta_t = tstat_t; SS1 = nan(T,1);

	% NOW LOOP TROUGHT t0:tT
	for tt = t0:t1 % note it is 4:(T-4) as i in 4:(T-5) means that in R
		ols_out	= fullols(GY,Tvec > tt);					% with constant as well as 'dummy' included in OLS
	  % 	Fstat_0(tt) = ols_out.Fstat;				  % same as tstat(end)^2, not used            
		tstat_t(tt,:)	= ols_out.tstat';						% tstat on the last regressor dummy variable
		beta_t(tt,:)	= ols_out.bhat';						% store TVP beta coefficents, exclude break dummy
		
		% CHOW TEST LIKE IN SW98 'BY HAND'. TEST FOR BREAK IN MEAN,(1 PARAMETER ONLY)
		e11 = demean( GY(1:tt) );
		e12 = demean( GY(1+tt:T) );
		SS1(tt) = e11'*e11 + e12'*e12;
	end
	
	% SW98 version of the test (2 degrees of freedom) Classic Chow version (ONLY TESTING THE MEAN)
	Fstat_chow	= (T-2)*(SS0 - SS1(t0:t1,:) )./SS1(t0:t1,:);
	% Wald versions of the test from running regression on intercept plus break dummy
	Fstat_ols		= removenan(tstat_t(:,end).^2);
	% plot([LR Fstat])
	
	%% CALL THE STRUCTURAL BREAK TEST FUNCTION & MUE
	Stats_chow		= get_stage1_break_stats(Fstat_chow,GY,T);
	Lambdas_chow	= get_stage1_lambda_MUE(Stats_chow,T);
	
	% FROM THE OLS REGRESSIONS
	Stats_ols			= get_stage1_break_stats(Fstat_ols,GY,T);
	Lambdas_ols		= get_stage1_lambda_MUE(Stats_ols,T);

	% RETURN ALSO SOME OTHER TERMS
	ExtraOutput.tstat_t			= tstat_t;
	ExtraOutput.beta_t			= beta_t;
	ExtraOutput.Fstat_ols		= Fstat_ols;	
	ExtraOutput.Fstat_chow	= Fstat_chow;
	ExtraOutput.T						= T;
end

function Stats_out = get_stage1_break_stats(Fstat,GY,T)
	% COMPUTE THE REMAINING STRUCTURAL BREAK TEST
	mu_Fstat	= mean(Fstat);
	exp_Wald	= log(mean(exp(Fstat/2)));
	sup_Fstat = max(Fstat);

	% Nyblom's test is simply: e_cumsum'e_cumsum/T^2 divided by Var(e), where e = xfilt-mean(xfilt).
	e0		= demean(GY);
	e_cs	= cumsum(e0)/sqrt(T);
	Lstat	= (e_cs'*e_cs/T)/var(e0); 

	% Combine tests for printing later (IF ORDER CHANGED, MUST CHANGE ORDER OF CRITICALVALS AS WELL)
	Stats_out.L		= Lstat;			% Nybloms L-test
	Stats_out.MW	= mu_Fstat;		% mean Wald
	Stats_out.EW	= exp_Wald;		% exponetial Wald
	Stats_out.QLR = sup_Fstat;	% QLR sup-LR test
end

function Lambdas_out = get_stage1_lambda_MUE(Stats,T)
	% NOW GET/LOAD THE CRITICAL LAMBDA VALUES FROM THE COARSE GRID AS THE TABLE IN SW1998 AND AS USE IN HLW2017.
	Lambda_table = (0:30)';
	% L				MW				EW				QLR            % corresponding Lambda values from the SW GRID 
	CritVals = [                                                                                
	0.118	 	 0.689		 0.426	 	 3.198           %     0                                        
	0.127	 	 0.757		 0.476		 3.416           %     1                                        
	0.137	 	 0.806		 0.516		 3.594           %     2                                        
	0.169	 	 1.015		 0.661		 4.106           %     3                                        
	0.205	 	 1.234		 0.826		 4.848           %     4                                        
	0.266	 	 1.632		 1.111		 5.689           %     5                                        
	0.327	 	 2.018		 1.419		 6.682           %     6                                        
	0.387	 	 2.390		 1.762		 7.626           %     7                                        
	0.490	 	 3.081		 2.355		 9.160           %     8                                        
	0.593	 	 3.699		 2.910		10.660           %     9                                        
	0.670	 	 4.222		 3.413		11.841           %    10                                        
	0.768	 	 4.776		 3.868		13.098           %    11                                        
	0.908	 	 5.767		 4.925		15.451           %    12                                        
	1.036	 	 6.586		 5.684		17.094           %    13                                        
	1.214	 	 7.703		 6.670		19.423           %    14                                        
	1.360	 	 8.683		 7.690		21.682           %    15                                        
	1.471	 	 9.467		 8.477		23.342           %    16                                        
	1.576	 	10.101		 9.191		24.920           %    17                                        
	1.799	 	11.639		10.693		28.174           %    18                                        
	2.016	 	13.039		12.024		30.736           %    19                                        
	2.127	 	13.900		13.089		33.313           %    20                                        
	2.327	 	15.214		14.440		36.109           %    21                                        
	2.569	 	16.806		16.191		39.673           %    22                                        
	2.785	 	18.330		17.332		41.955           %    23                                        
	2.899	 	19.020		18.699		45.056           %    24                                        
	3.108	 	20.562		20.464		48.647           %    25                                        
	3.278	 	21.837		21.667		50.983           %    26                                        
	3.652	 	24.350		23.851		55.514           %    27                                        
	3.910	 	26.248		25.538		59.278           %    28                                        
	4.015	 	27.089		26.762		61.311           %    29                                        
	4.120	 	27.758		27.874		64.016           %    30                                        
	];

	% % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	% % ALTERNATIVELY: LOAD THE FINER GRID FROM ORIGINAL SW FILES AS USED IN SW98
	% % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	% SWdatadir = './_data/SW98.inputs/';
	% load([SWdatadir 'SW1998_MUE_lookup_table.mat']);
	% 
	% MUE_vals_tmp	= lookup_table.Variables;
	% Lambda_table	= MUE_vals_tmp(:,1);
	% CritVals			= MUE_vals_tmp(:,2:5);
	% % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	%% Final calculations for clean output
	testStats = struct2array(Stats);
	% testNames = {'L';'MW';'EW';'QLR'}; % SW = QLR
	testNames = fieldnames(Stats);
	% fprintf(' Test is: %s \n', char(testNames(ii)))

	IL_MAX = length(Lambda_table);

	% loop through the 4 main tests
	for ii = 1:4
		test_stat = testStats(ii);
		crit_val	= CritVals(:,ii);
		% get indicator
		isGreater = test_stat > crit_val;
		IL1 = find(isGreater,	1, 'Last');
		IL0 = find(isGreater);
		IL	= max(IL1, max(IL0));

		if(isempty(IL))
			Lambda_hat(ii,1) = 0;
		elseif(IL==(IL_MAX)) % if we are at the last table entry, just return that without interpolation
			Lambda_hat(ii,1) = Lambda_table(IL_MAX) / T; % if IL is greater than that largest lambda value of 30, return 30/T
		else
			% INTERPOLATE TWO ADJACENT OBSERVATIONS AND DIVIDE BY T. THIS IS A DOUBLE
			Lambda_hat(ii,1) = interp1(crit_val(IL:IL+1), Lambda_table(IL:IL+1), test_stat) / T;
		end
		% MAKE A STRUCTURE WITH THE FIELDS EQUAL TO THE TESTNAMES 
		Lambdas_out.(testNames{ii}) = Lambda_hat(ii);
	end
end

























% EOF 




















































