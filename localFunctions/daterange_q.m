function date_range_out = daterange_q(beg_date, end_date)
% F: returns the date range needed for a timetable date adjustment.
% Put in either datenum begin and end dates, datetime inputs, or quarter
% strings, i.e. 'Q1-1997' or '1997:Q1'.
%
% CALL AS:	daterange_q(date_num_beg, date_num_date) or
%						daterange_q('Q1-1997', 'Q3-2001').
% In a call to a timetable variable, e.g.
%		us_data(daterange_q('Q1-1997', 'Q3-2001'),:)	to truncate for all variables.
% ---------------------------------------------------------------------------------------------

T0 = to_quarter_datetime(beg_date);
T1 = to_quarter_datetime(end_date);

date_range_out = timerange(T0, T1, 'closed');

end

% ---------------------------------------------------------------------------------------------
function dt_out = to_quarter_datetime(val_in)
% Convert numeric, datetime, or quarter-string inputs to a datetime value that
% points to the first day of the corresponding quarter.

if isnumeric(val_in)
	dt_out = datetime(val_in, 'ConvertFrom', 'datenum');
elseif isdatetime(val_in)
	dt_out = val_in;
elseif isduration(val_in) || iscalendarduration(val_in)
	error('daterange_q:InvalidType', 'Durations are not supported. Provide datenum, datetime, or quarter text.');
else
	val_str = strtrim(string(val_in));
	dt_out = quarter_text_to_datetime(val_str);
end
end

function dt_out = quarter_text_to_datetime(val_str)
% Parse quarter text such as 'Q1-1997', '1997:Q1', or '1997Q1'.

val_str = regexprep(val_str, '[\s_]', '');

% try Q#-YYYY, Q#YYYY, Q#/YYYY
tokens = regexp(val_str, '^[Qq]([1-4])[-/:]?(\d{4})$', 'tokens', 'once');

% try YYYY[-/:]?Q#
if isempty(tokens)
	tokens = regexp(val_str, '^(\d{4})[-/:]?[Qq]([1-4])$', 'tokens', 'once');
	if ~isempty(tokens)
		tokens = fliplr(tokens);
	end
end

if isempty(tokens)
	error('daterange_q:UnrecognizedQuarterFormat', ...
				'Cannot parse quarter string "%s". Expected formats like "Q1-1997" or "1997:Q1".', val_str);
end

quarter_num = str2double(tokens{1});
year_num    = str2double(tokens{2});

month_num = (quarter_num-1)*3 + 1;
dt_out    = datetime(year_num, month_num, 1);
end
