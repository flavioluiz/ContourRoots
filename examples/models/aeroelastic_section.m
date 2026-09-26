function model = aeroelastic_section(kind)
%AEROELASTIC_SECTION Published 2-DOF typical-section benchmarks, per span.
%   MODEL = AEROELASTIC_SECTION(KIND) returns a structure with the mass
%   matrix M, stiffness matrix K, structural damping Cs (zero), semichord b,
%   elastic-axis position a, air density rho, and the published reference
%   flutter values. Edit the fields to study other sections.
%   'nasa' (default): Perry, NASA/TP-2015-218765, Appendix C, standard case.
%     Consistent normalized units b=1, rho=1, omega_alpha=100. Velocities
%     numerically correspond to ft/s when b=1 ft, as in the NASA report.
%   'dlr': Kaiser & Quero (2022), Table 1, SI units; flutter 212.2 m/s.
%   h is positive downwards, alpha nose-up; a is the elastic-axis position
%   aft of midchord divided by b. Loads are [-lift; nose-up moment].
    if nargin==0, kind='nasa'; end
    switch lower(kind)
        case 'nasa'
            b=1; rho=1; m=10*pi*rho*b^2; a=-.4;
            S=.2*m*b; I=.25*m*b^2; kh=m*50^2; ka=I*100^2;
            ref=struct('U',173.26,'k',.4355,'velocityTolerance',.08, ...
                'frequencyTolerance',4e-4,'source','NASA/TP-2015-218765, Appendix C, Table CI');
            units='ft/s';
        case 'dlr'
            b=1; rho=1.225; m=292.4823; S=73.1206; I=113.482;
            kh=9.1396e5; ka=4.1965e5; a=-.15;
            ref=struct('U',212.2,'k',NaN,'velocityTolerance',.15, ...
                'frequencyTolerance',NaN,'source','Kaiser & Quero (2022), Table 1 / Figure 2');
            units='m/s';
        otherwise
            error('aeroelastic:Model','Use nasa or dlr.');
    end
    model=struct('name',lower(kind),'b',b,'rho',rho,'a',a, ...
        'M',[m S;S I],'K',diag([kh ka]),'Cs',zeros(2), ...
        'reference',ref,'velocityUnits',units);
end
