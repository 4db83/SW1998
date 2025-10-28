function olsout = ols(y, X, no_const, Xnames, L, INCLUDE_PW, print_results, Recession_indicator)
% function olsout = ols(Y,X, no_const, Xnames, L, INCLUDE_PW, print_results, Recession_indicator)

% function olsout = ols(1  2         3       4  5  6           7              8)
% set L = NaN to use plain OLS Standard erros, L = 0 is White, L = [] HAC optmal truncation lag
% olsout = ols(y,X,L,INCLUDE_PW,Recession_indicator,no_const);
% wrapper function to fullols to print results to screen.

SetDefaultValue(3,'no_const',0)						% default is to include a constant
SetDefaultValue(4,'Xnames',[]);           % cellstring of regressor names
SetDefaultValue(5,'L',[]);								% set default value for truncation lag to L = 
SetDefaultValue(6,'INCLUDE_PW',[]);				% set to 1 to include pre-whitening in LRV computation
SetDefaultValue(7,'print_results',1)			% print results to screen is default
SetDefaultValue(8,'Recession_indicator',[]);

ytt = isa(y,"timetable");

if ytt 
  dates = y.Properties.RowTimes;
  y     = y.Variables; 
end

xtt = isa(X,"timetable"); 
if xtt
  if isempty(Xnames) 
    Xnames  = X.Properties.VariableNames;
  else 
    Xnames = Xnames;
  end
  dates   = X.Properties.RowTimes;
  X       = X.Variables; 
end

% chec if L = nan
isnanL = isnan(L);

II = anynans(y,X);
% olsout = fullols(y,x,L,INCLUDE_PW,Recession_indicator,no_const);
if isnanL
  olsout = fullols(y(~II),X(~II,:),no_const, [], INCLUDE_PW, Recession_indicator);
  % olsout = fullols(y(~II),x(~II,:),no_const, L, INCLUDE_PW, Recession_indicator);
else
  olsout = fullols(y(~II),X(~II,:),no_const, L,  INCLUDE_PW, Recession_indicator);
end

% add the dates vector
if ytt || xtt
  Dates = dates(~II);
  olsout.dates = Dates;
end

% sqrt(diag(olsout.HAC_VCV))
% sqrt(diag(olsout.HAC_VCV_pw_ARMA11))
% sqrt(diag(olsout.HAC_VCV_pw_AR1))

if ~isempty(Xnames)
  if ~iscell(Xnames); Xnames = cellstr(Xnames); end
end

% print sample period if it exists. 
if ytt || xtt
  fprintf( ['     Estimation period:  ' date_string(Dates(1),'m') ' to ' date_string(Dates(end),'m') '   Sample Size N = ' num2str(olsout.N) '\n' ])
end

if isnanL
  prnt_hac = 0;
else
  prnt_hac = [];
end

if ischar(print_results)
	print_fullols(olsout, Xnames, [], print_results);
elseif print_results==1
	print_fullols(olsout, Xnames, prnt_hac);
% 	print_fullols(olsout, Xnames, [], print_results);
%   print_fullols(fullols_output,variable_names,prnt_hac,prnt2xls)
end


