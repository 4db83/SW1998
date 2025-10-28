function YX_Stage2 = make_yx_4_stage2_M0g_unrestricted(YY, XX, Stage2_KFS_object)
% Function: makes the Y and X variables needed for M0.HLWs Stage 2 MUE inputs in unrestricted from.
% ---------------------------------------------------------------------------------------------
% INPUT: YY, XX data matrices, and Stage2_KFS_object obtained from running 
% [~,Stage2_KFS_object] = LogLike_Stage2_HLW_SSF(theta_db, data_in_with_intercept, other_paras, 1);
% OUTPUT: 
% YX_Stage2 structure with tt and tT elements referring to filtered and smoothed
% estimates.
% where - YY =	[y~(t)]
% and 
%				- XX =	[y~(t-1) y~(t-1) r(t-1) r(t-2) g(t-1) g(t-2) 1] (1 is an intercept term), 
%								g(t) NOT annualized.
% db
% --------------------------------------------------------------------------------------------- 

RR_2_lags = XX(:,3:4);
T = length(Stage2_KFS_object.D);

% MAKE THE STAGE 2 Y AND X VARIABLES FROM SMOOTHED KFS OUTPUT
% smoothed potential or y*(t) 
% note: xi(t) = [y*(t) y*(t-1) y*(t-2) g(t-1) g(t-2)], so take the first two row entries in y*(t-2) and stack on top of y*(t) 
ystar_tT	= [Stage2_KFS_object.KFS.atT(1:2,3); Stage2_KFS_object.KFS.atT(:,1)];
% smoothed output gap (or cycle) no make the cycle and then lag two times. This exactly
% reproduces HLW
cycle_tT  = YY(3:end,1)	-	ystar_tT;

% lag the cycle by two periods
cycle_2_lags = mlag(cycle_tT,2);

% Y = [ y~(t) ] output gap from 1960:Q3 (1961:Q1 is the start of the in-sample
YStage2 = cycle_tT(3:end);

% X = [y~(t-1) y~(t-1) r(t-1) r(t-2) g(t-1) g(t-2) 1]
XStage2 = [	cycle_2_lags(3:end,:)							, ... % y~_(t-1)	y~_(t-2)
	          RR_2_lags(5:end,:)								, ... % r_(t-1)		r_(t-2) 
	          Stage2_KFS_object.KFS.atT(:,4:5)	, ...	% g_(t-1)		g_(t-2)		% this is not annualized in any way.
	          ones(T,1)];															% Intercept term

YX_Stage2.tT = [YStage2 XStage2];					

% MAKE THE STAGE 2 Y AND X VARIABLES FROM FILTERED KFS OUTPUT
% filtered potential or y*(t)
ystar_tt	= [Stage2_KFS_object.KFS.att(1:2,3); Stage2_KFS_object.KFS.att(:,1)];
% filtered output gap (or cycle)
cycle_tt  = YY(3:end,1)	-	ystar_tt;

% lag the cycle by two periods
cycle_2_lags = mlag(cycle_tt,2);

% Y = [ y~(t) ] output gap from 1960:Q3 (1961:Q1 is the start of the in-sample
YStage2_tt = cycle_tt(3:end);

% X = [y~(t-1) y~(t-1) r(t-1) r(t-2) g(t-1) g(t-2) 1]
XStage2_tt = [	cycle_2_lags(3:end,:)							, ... % y~_(t-1)	y~_(t-2)
								RR_2_lags(5:end,:)								, ... % r_(t-1)		r_(t-2) 
								Stage2_KFS_object.KFS.att(:,4:5)	, ...	% g_(t-1)		g_(t-2)		% this is not annualized in any way.
								ones(T,1)];															% Intercept term

YX_Stage2.tt = [YStage2_tt XStage2_tt];										