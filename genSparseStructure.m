function [rows,cols] = genSparseStructure(params)
% GETSPARSESTRUCTURE(params) returns a sparse matrix containing ones at
%   every matrix location that can be non-zero. The constraints matrix C is
%   such that each constraint depends on w_i and w_{i+1}, so many of its
%   elements will be zero.
% -------------------------------------------------------------------------
% INPUTS:
%   params - a structure of problem parameters (masses, inertias, lengths,
%       etc.) 

    nw = params.nw; nu = params.nu; N = params.N;
    dCeq = zeros(nw*(N-1),N*(nw+nu));

    part = ones(nw,2*(nw+nu));

    for ii=1:(N-1)
        x_inds = (1:2*(nw+nu)) + (nw+nu)*(ii-1);
        y_inds = (1:nw) + nw*(ii-1);

        dCeq(y_inds,x_inds) = part;
    end

    dCeq = sparse(dCeq);
    
    [rows,cols,~] = find(dCeq); 
end