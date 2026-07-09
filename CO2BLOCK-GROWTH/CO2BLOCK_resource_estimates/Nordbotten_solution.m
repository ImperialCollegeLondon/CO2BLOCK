
function [PD] = Nordbotten_solution(r,R,psi,R_ext,gamma)

x = max(r,psi);
FD = (x < R_ext) * (log(min(R,R_ext)/x) + (R > R_ext)*(2/2.25*(R/R_ext)^2 -3/4));
% if (x >= R_ext)
%     assert(FD == 0);
% end
PD = ((r <= psi) || (r <= R)) * (gamma*log(max(psi/r, 1)) + FD);
% if ((r > psi) && (r > R))
%     assert(PD == 0);
% end
end
