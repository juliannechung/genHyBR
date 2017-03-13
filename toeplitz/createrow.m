function [Qr, pts] = createrow(ptsmin, ptsmax, N, kernel, theta)
% Example usage 
% Qr = createrow(0,1,100,@(r)exp(-r),1.0);
% dim = 3; Qr = CreateRow(zeros(3,1),ones(3,1), ...
%                   [20,20,10],@(r) exp(-r.^2./5^2), ones(dim,1));
%
% Distance r (b/w two points x and y) is defined as follows
% r(x,y)   - sqrt(sum_{i=1}^d (x_i - y_i).^2/theta_i.^2)
% 
% Input:
% N        - the number of points in each dimension [dim x 1]
% kernel   - radially symmetric covariance kernel 
% Output:
% Qr 	   - the row of the Covariance matrix
% pts      - the points involved in the grid
%
% Written by Arvind K. Saibaba arvind.s.krishna@gmail.com


    dim = numel(ptsmin);
    assert(dim == numel(ptsmax), 'posmin and posmax are incosistent');
    assert(dim == numel(theta),  'Check dimensions of theta');
    
    switch (dim)
        
        case 1  %1D
	    x = linspace(ptsmin(1), ptsmax(1), N(1))';
            
        case 2  %2D

	    x1 = linspace(ptsmin(1), ptsmax(1), N(1));
	    x2 = linspace(ptsmin(2), ptsmax(2), N(2));
	    [xx,yy] = meshgrid(x1,x2);
	    x = [xx(:) yy(:)];

        case 3  %3D 
	    
	    x1 = linspace(ptsmin(1), ptsmax(1), N(1));
	    x2 = linspace(ptsmin(2), ptsmax(2), N(2));
	    x3 = linspace(ptsmin(3), ptsmax(3), N(3));

	    [xx,yy,zz] = ndgrid(x1,x2, x3);
	    x = [xx(:) yy(:) zz(:)];

	otherwise
		error('Wrong dimension');
        
    end
    
    % Compute the distances 
    R = distances(x,x(1,:),theta);
    Qr = feval(kernel, R)';

    % Assign points
    pts = x;
    
end

function R = distances(x,y,theta)
    dim = size(x,2);

    assert(dim == size(x,2), 'Check dimensions of x and y');
    assert(dim == numel(theta), 'Check dimensions of theta');

    r = zeros(size(x,1),1);	 	
    for i = 1:dim
	r = r + (x(:,i) - y(1,i)).^2./theta(i).^2;
    end
    R = sqrt(r);	

end
