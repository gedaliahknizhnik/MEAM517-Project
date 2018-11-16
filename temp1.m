function [rows,cols] = genSparseStructure(nx,nu,N)

    dCeq = zeros(nx*(N-1),N*(nx+nu));

    part = ones(nx,2*(nx+nu));

    for ii=1:(N-1)
        x_inds = (1:2*(nx+nu)) + (nx+nu)*(ii-1);
        y_inds = (1:nx) + nx*(ii-1);

        dCeq(y_inds,x_inds) = part;
    end

    dCeq = sparse(dCeq);
    
    [rows,cols,~] = sparse(dCeq); 
    
end