function q_dates = make_quarterly_dates(beg_Year,beg_Quarter,end_Year,end_Quarter)
% ------------------------------------------------------------------------------------------------------
% FUNCTION: makes a quarterly date vector.
% ------------------------------------------------------------------------------------------------------
% CALL AS: 
%			q_dates = make_quarterly_dates(beg_Year,beg_Quarter,end_Year,end_Quarter)
% ------------------------------------------------------------------------------------------------------
% INPUT: 
% beg_Year			= beginning year as number, ie. 1947.
% beg_Quarter		= beginning quarter as number, ie., 3.
% end_Year			= beginning year as number, ie. 2010.
% end_Quarter		= beginning quarter as number, ie., 4.
% 
% OUTPUT:
% q_dates				= a quarterly date vector with FRED quarterly date convention, that is:
%									Q1 = 1st January
%									Q2 = 1st April.
%									Q3 = 1st July.
%									Q4 = 1st October.
%									
% NOTES:
%			To convert to end of quarter day dates, ie., Q1 = 31st of March, etc., call the function, 
%			make_end_of_quarter_dates.m
% ------------------------------------------------------------------------------------------------------
% db (01.03.2018)
% ------------------------------------------------------------------------------------------------------

% make the quarterly date vector
dates_tmp = bsxfun(@(Month,Year) datenum(Year,Month,1),(1:3:12).',beg_Year:end_Year);
q_date0		= dates_tmp(:);

% return quarterly date vector
q_dates = q_date0(beg_Quarter:end-(4-end_Quarter));

