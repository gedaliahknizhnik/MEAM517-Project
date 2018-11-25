function [xdNow,udNow] = getDes(t,td,xd,ud,polys)
    
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
    
    xdNow = interp1(td,xd,t)';
    udNow = interp1(td,ud,t);
end