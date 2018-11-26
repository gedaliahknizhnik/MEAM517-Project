function [w_star,u_star] = getDes(t,t_stars,w_stars,u_stars,polys)
% GETDES(t,t_stars,w_stars,u_stars,polys) finds the the optimal 
%   trajectory point corresponding to the current time t.
% ---------------------------------------------------------------
% INPUTS:
%   t - the current time
%   t_stars - a vector of times corresponding to the SNOPT outputs w_stars
%       and u_stars
%   w_stars - a matrix of SNOPT state outputs
%   u_stars - a vector of SNOPT control outputs
%   polys - the polynomial trajectory output by SNOPT (currently unused)
% OUTPUTS:
%   w_star - the current optimal trajectory point.
%   u_star - the current optimal control point.

    % If using polynomial trajectories - figure out which polynomial to use and
    % evaluate it.
    %     for jj=size(td):-1:2
    %         if t <= td(jj) 
    %             ind = jj-1;
    %             break
    %         end
    %     end
    % 
    %     dt = t - td(ind);
    %     dts = [dt^3;dt^2;dt^1;1];
    %     xdNow = squeeze(polys(ind,:,:))'*dts;
 
    % Linearly interpolate in both state and control. (Linearly
    % interpolating in state is wrong but gives better results...)
    w_star = interp1(t_stars,w_stars,t)';
    u_star = interp1(t_stars,u_stars,t);
end