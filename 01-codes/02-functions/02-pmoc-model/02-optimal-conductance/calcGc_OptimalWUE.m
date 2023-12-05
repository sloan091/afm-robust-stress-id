% =========================================================================
% Name   : calcMedlynGc.m
% Author : Brandon Sloan
% Date   : 3/25/21
%
% DESCRIPTION
% This function calculates the canopy conductance (G_c) using the
% formulations laid out in either Lin et al. (2018) or Medlyn et al.
% (2017).
%
% INPUTS
%   G_o   - Minimum stomatal conductance for H2O [moles air/m^2/s]
%   G_1   - Empirical Medlyn parameter [kPa^0.5]
%   m     - Exponent of VPD_l [-]
%   VPD_l - Vapor pressure difference between inside the leaf and at
%           the leaf surface [kPa]
%   GPP   - Estimate of net CO2 assimilation [micromoles CO2/m^2/s]
%   C_a   - CO2 concentration [micromoles CO2/moles air]
%
% OUTPUTS
%   G_c   - Canopy conductance of H20 [moles air/m^2/s]
%
% REFERENCES
%   (1) Lin, C., Gentine, P., Huang, Y., Guan, K., Kimm, H., & Zhou, S.
%   (2018). Diel ecosystem conductance response to vapor pressure deficit
%   is suboptimal and independent of soil moisture. Agricultural and
%   Forest Meteorology, 250–251, 24–34.
%   https://doi.org/10.1016/J.AGRFORMET.2017.12.078
%
%   (2) Medlyn, B. E., De Kauwe, M. G., Lin, Y.-S., Knauer, J., Duursma,
%   R. A., Williams, C. A., et al. (2017). How do leaf and ecosystem
%   measures of water-use efficiency compare? New Phytologist, 216(3),
%   758–770. https://doi.org/10.1111/nph.14626
% =========================================================================

function [G_c,G_1] = calcGc_OptimalWUE(G_o,G_1ww,m,a,VPD_l,GPP,C_a,SM,...
    GcType,StressType)

% Choose static or dynamic stress on G1
switch StressType
    case 0 % Well-watered
        G_1 = G_1ww;
    case 1 % Static
        G_1 = G_1ww.*exp(-a.*abs(max(SM) - SM));
    case 2 % Dynamic Beta stress formulation 1
        G_1 = G_1ww.*exp(-a(1).*abs(max(SM) - SM) - a(2)*VPD_l);
    case 3 % Dynamic Beta stress formulation 2
        G_1 = G_1ww.*exp(-(a(1) + a(2)*VPD_l).*abs(max(SM) - SM));
    case 4 % Static HESS version
        G_1 =  G_1ww.*2.^(-(abs(SM)./a(1)).^a(2));
    case 5 % Dynamic HESS version
        G_1 =  G_1ww.*2.^(-((abs(SM)  + (VPD_l - mean(VPD_l(:),'omitnan'))...
            *a(3))./a(1)).^(a(2)));
end

% Stomatal conductance of H20 [moles air/m^2/s]
switch GcType
    case 1 % Use Lin et al. (2018) G_c
        G_c = G_o + (G_1.*GPP)./(C_a.*VPD_l.^m);
        
    case 2 % Use Medlyn et al. (2017) G_c
        G_c = G_o + (1.6*GPP./C_a).*(1 + G_1./(VPD_l.^m));
        
end

end






