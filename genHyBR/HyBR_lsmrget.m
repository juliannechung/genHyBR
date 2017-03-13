function val = HyBR_lsmrget(options,name,default,flag)
%
%   VAL = HyBR_lsmrget(options,name,default,flag)
%
%   HyBR_lsmrget gets options parameters.
%   VAL = HyBR_lsmrget(OPTIONS,'NAME') extracts the value of the named parameter
%   from HyBR_lsmr options structure OPTIONS, returning an empty matrix if
%   the parameter value is not specified in OPTIONS.  It is sufficient to
%   type only the leading characters that uniquely identify the
%   parameter.  Case is ignored for parameter names.  [] is a valid OPTIONS
%   argument.
%
%   VAL = HyBR_lsmrget(OPTIONS,'NAME',DEFAULT) extracts the named parameter as
%   above, but returns DEFAULT if the named parameter is not specified
%   in OPTIONS.  For example
%
%     param = HyBR_lsmrget(opts,'RegPar','WGCV');
%
%   returns param = 'WGCV' if the RegPar property is not specified in opts.
%
%   See also HyBR_lsmrset and HyBR_lsmr.
%
%   J.Chung 3/2014
%

% undocumented usage for fast access with no error checking
if (nargin == 4) && isequal(flag,'fast')
    val = HyBR_lsmrgetfast(options,name,default);
    return
end

if nargin < 2
    error('Not enough input arguments.');
end
if nargin < 3
    default = [];
end

if ~isempty(options) && ~isa(options,'struct')
    error('First argument must be an options structure created with HyBR_lsmrset.');
end

if isempty(options)
    val = default;
    return;
end
allfields = {'InSolv'; 'RegPar';'nLevel';'Omega';'Iter';'Reorth'; ...
    'x_true';'BegReg'; 'FlatTol'; 'MinTol'; 'ResTol'};

Names = allfields;

name = deblank(name(:)'); % force this to be a row vector
j = find(strncmpi(name,Names,length(name)));
if isempty(j)               % if no matches
    error(['Unrecognized property name ''%s''.  ' ...
        'See HyBR_lsmrset for possibilities.'], name);
elseif length(j) > 1            % if more than one match
    % Check for any exact matches (in case any names are subsets of others)
    k = find(strcmpi(name,Names));
    if length(k) == 1
        j = k;
    else
        msg = sprintf('Ambiguous property name ''%s'' ', name);
        msg = [msg '(' Names{j(1),:}];
        for k = j(2:length(j))'
            msg = [msg ', ' Names{k,:}];
        end
        msg = sprintf('%s).', msg);
        error(msg);
    end
end

if any(strcmp(Names,Names{j,:}))
    val = options.(Names{j,:});
    if isempty(val)
        val = default;
    end
else
    val = default;
end

%------------------------------------------------------------------
function value = HyBR_lsmrgetfast(options,name,defaultopt)
%HYBRGETFAST- Get HyBR_lsmr OPTIONS parameter with no error checking.
%   VAL = HyBR_lsmrGETFAST(OPTIONS,FIELDNAME,DEFAULTOPTIONS) will get the
%   value of the FIELDNAME from OPTIONS with no error checking or
%   fieldname completion. If the value is [], it gets the value of the
%   FIELDNAME from DEFAULTOPTIONS, another OPTIONS structure which is
%   probably a subset of the options in OPTIONS.
%
if isempty(options)
     value = defaultopt.(name);
     return;
end
% We need to know if name is a valid field of options, but it is faster to use 
% a try-catch than to test if the field exists and if the field name is
% correct.
try
    value = options.(name);
catch
    value = [];
    lasterr('');  % clean up last error
end

if isempty(value)
    value = defaultopt.(name);
end


