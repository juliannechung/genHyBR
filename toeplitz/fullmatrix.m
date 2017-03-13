function Q = fullmatrix(ptsmin, ptsmax, N, kernel, theta)

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
    R = distances(x,x,theta);
    Q = feval(kernel, R)';


    
end

function R = distances(x,y,theta)
    dim = size(x,2);

    assert(dim == size(x,2), 'Check dimensions of x and y');
    assert(dim == numel(theta), 'Check dimensions of theta');

    r = zeros(size(x,1),size(y,1));	 	
    for i = 1:dim
        [xx,yy] = meshgrid(x(:,i),y(:,i));
        r = r + (xx - yy).^2./theta(i).^2;
    end
    R = sqrt(r);	

end