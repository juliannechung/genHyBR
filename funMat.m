classdef funMat
  %funMat class
  %  A funMat object is used to represent a matrix A,
  %       where A is accessed via function evaluations for both
  %       matrix-vector and matrix-transpose-vector multiplications
  %
  %  The funMat class has inputs:
  %     Afun  - matvec function handle (required)
  %     Atfun - transpose function handle (required)
  %     Asize - size of the matrix (required)
  %     Ainvfun - inverse function handle (optional)
  %     Atinvfun - transpose inverse function handle (optional)
  %
  %   and is based on a structure with the following fields:
  %     Afun
  %     Atfun
  %     Asize
  %     Ainvfun
  %     Atinvfun
  %     transpose - indicates if the matrix has been transposed
  %
  %  Calling Syntax:
  %
  %    P = funMat    (returns object with empty fields)
  %    P = funMat(funMat) (returns funMat)
  %    P = funMat(Afun, Atfun)
  %    P = funMat(Afun, Atfun, Asize)
  %    P = funMat(Afun, Atfun, Asize, Ainvfun, Atinvfun)
  %
  % J. Chung, 3/4/2016
  
  properties
    Afun
    Atfun
    Asize
    Ainvfun
    Atinvfun
    transpose
  end
  
  methods
    function P = funMat(varargin) % Constructor
      switch nargin
        case 0
          P.transpose = false;
          P.Afun = [];
          P.Atfun = [];
          P.Ainvfun = [];
          P.Atinvfun = [];
        case 1
          if isa(varargin{1}, 'funMat')
            P = varargin{1};
          else
            error('Incorrect input arguments')
          end
        otherwise
          P.transpose = false;
          if nargin >= 3
            P.Afun = varargin{1};
            P.Atfun = varargin{2};
            P.Asize = varargin{3};
            if nargin >= 4
              P.Ainvfun = varargin{4};
            end
            if nargin == 5
              P.Atinvfun = varargin{5};
            end
          else
            error('Incorrect number of input arguments')
          end
      end
    end
    
    function P = ctranspose(P) % Overload transpose
      P.transpose = not(P.transpose); % switches booleen transpose flag
    end
    
    function y = mtimes(arg1, arg2) % Overload matrix vector multiplicaiton
      
      if isa(arg1,'funMat') && isa(arg2,'double')
        %   Implement A*s and A'*s for funMat object A.
        if isscalar(arg2)
          error('Matrix times a scalar not implemented yet')
        elseif isvector(arg2)          
          if arg1.transpose 
            % Matrix transpose times vector
            y = arg1.Atfun(arg2);
          else
            % Matrix times vector
            y = arg1.Afun(arg2);
          end
        elseif ismatrix(arg2)
          if arg1.transpose 
            % Matrix transpose times each column of the matrix: A'*M
            for i = 1:size(arg2,2)
              y(:,i) = arg1.Atfun(arg2(:,i));
            end
          else
            % Matrix times each column of the matrix A*M
            for i = 1:size(arg2,2)
              y(:,i) = arg1.Afun(arg2(:,i));
            end
          end
          
        end
        
      elseif isa(arg1,'double') && isa(arg2,'funMat')
        %   Implement s'*A and s'*A' for funMat object A.
        if isscalar(arg1)
          error('Scaler times a matrix not implemented yet')
        
        elseif isvector(arg1)          
          if arg2.transpose 
            % Vector transpose times funMat transpose x*A' = (A*x')'
            y = arg2.Afun(arg1')';
          else
            % Vector transpose times funMat x*A = (A'*x')'
            y = arg2.Atfun(arg1')';
          end
          
        elseif ismatrix(arg1)
          MT = arg1';
          if arg2.transpose 
            % Matrix transpose times funMat transpose M*A' = (A*M')'
            for i = 1:size(MT,2)
              y(:,i) = arg2.Afun(MT(:,i));
            end
          else
            % Vector transpose times funMat M*A = (A'*M')'
            for i = 1:size(MT,2)
              y(:,i) = arg2.Atfun(MT(:,i));
            end
          end
          y = y';
        end
        
      else
        error('Multiplication is not implemented.')
      end
    end
    
    function y = mldivide(arg1, arg2) % Overload matrix inverse       
      if arg1.transpose 
        if isempty(arg1.Atinvfun)
          error('Need to specify a transpose inverse function handle.')
        end
      else
        if isempty(arg1.Ainvfun)
          error('Need to specify an inverse function handle.')
        end
      end
      if isa(arg1,'funMat') && isa(arg2,'double')
        %   Implement A\s and A'\s for funMat object A.
        if isscalar(arg2)
          error('mldivide does not work for scalars yet')
        elseif isvector(arg2)          
          if arg1.transpose 
            % Matrix transpose inverse times vector
            y = arg1.Atinvfun(arg2);
          else
            % Matrix inverse times vector
            y = arg1.Ainvfun(arg2);
          end
        elseif ismatrix(arg2)
          if arg1.transpose 
            % Matrix inverse times each column of the matrix: A' \ M
            for i = 1:size(arg2,2)
              y(:,i) = arg1.Atinvfun(arg2(:,i));
            end
          else
            % Matrix times each column of the matrix A \M
            for i = 1:size(arg2,2)
              y(:,i) = arg1.Ainvfun(arg2(:,i));
            end
          end
        end
      else
        error('mldivide is not implemented for this case.')
      end
    end
    
    
    function varargout = size(A, dim) % Overload size
      d = A.Asize;
      if nargin == 2
        d = d(dim);
      end
      
      if nargout == 1 || nargout == 0
        varargout{1} = d;
      else
        for i = 1:length(d)
          varargout{i} = d(i);
        end
      end
      
    end
        
  end % methods
end % classdef

