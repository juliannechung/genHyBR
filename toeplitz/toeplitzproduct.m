function result = toeplitzproduct(x, Qr, N)
% Computes Q*x, when Q is Toeplitz, BTTB or BBTBTBTB
% 
% Input 
%  x    - the vector to be multiplied 
%  Qr   - the first row of the matrix, computed using CreateRow function
%  N    - the #elements in each dimension, dim x 1
% 
% Output
% result - the result of Q*x

% Written by Arvind K. Saibaba arvind.s.krishna@gmail.com

    dim = numel(N);

    assert(numel(x) == prod(N(1:dim)), 'Check dimensions of vector');
    assert(numel(Qr) == prod(N(1:dim)), 'Check dimensions of Toeplitz Row');
       
    
    switch dim
        case 1
            circ = [Qr, Qr((end-1):-1:2)];
            
            padded = [x; zeros(N(1)-2,1)];
        
            result = ifft(fft(circ').*fft(padded));
            result = result(1:N(1));	
        
        case 2
            circ = reshape(Qr, N(2), N(1));
            circ = [circ, circ(:,(end-1):-1:2)];
            circ = [circ; circ((end-1):-1:2,:)];
            
            padded = reshape(x, N(2), N(1));
            [n1,n2] = size(circ);
            
            
            result = ifft2(fft2(circ).*fft2(padded,n1,n2));
            result = result(1:N(2),1:N(1));
            result = reshape(result, prod(N(1:2)),1);
            
        case 3
            circ = reshape(Qr, N);
            circ = cat(1, circ, circ((end-1):-1:2,:,:));
            circ = cat(2, circ, circ(:,(end-1):-1:2,:));
            circ = cat(3, circ, circ(:,:,(end-1):-1:2));
            
            padded = reshape(x, N);
            [n1,n2,n3] = size(circ);
            
            result = ifftn(fftn(circ).*fftn(padded,[n1,n2,n3]));
            result = result(1:N(1),1:N(2),1:N(3));
            result = reshape(result, prod(N(1:3)),1);
            
            
        otherwise
            error('Wrong Dimension')
    end   	



end
