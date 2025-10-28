function date_range_out = daterange_q(beg_date, end_date)
% F: returns the data range need for a timetable date adjustment
% put in either datenum begin and end dates, or quarter strings, ie., QQ-yyyy.
%  
% CALL AS:	daterange_q(date_num_beg, date_num_date) or 
%						daterange_q('Q1-1997', 'Q3-2001').
% in a call to a timetable variable: ie.
%
% USE AS:
%		us_data(daterange_q('Q1-1997', 'Q3-2001'),:)	to truncate for all variables
% or 
% 
% db
% ---------------------------------------------------------------------------------------------


if isnumeric(beg_date)
	T0 = datestr( beg_date	,'QQ-yyyy');
else
	T0 = beg_date;
end

if isnumeric(end_date)
	T1 = datestr( end_date	,'QQ-yyyy');
else
	T1 = end_date;
end


date_range_out = timerange(T0, T1, 'closed');


% if isnumeric(beg_date)
% 	date_range_out = timerange( datestr( beg_date	,'QQ-yyyy'), ...
% 															datestr( end_date	,'QQ-yyyy'),'closed');
% else
% end	
	
	