% function dedirs_ge_31 gives DE_scheme with b0s if amount of b0s is given
% SuS April 2008 (copied from original tensor.dat file)


function DE_scheme = dedirs_ge_31(nob0s)

if nargin == 0
    DE_scheme = [-0.404485801990 0.801145160285 0.441086916763; ...
        -0.748942567232 -0.662363704264 -0.018956641627; ...
        -0.236315790505 -0.912424523381 -0.334120242265; ...
        0.693310892542 0.482245749539 0.535498873327; ...
        -0.203714420443 -0.442765397309 0.873189119177; ...
        0.399692451068 -0.069366022027 -0.914020951372; ...
        0.245878456697 0.627777017425 -0.738538963715; ...
        0.847226441419 -0.286749277133 -0.447193704142; ...
        0.234753426819 0.909610476648 0.342782160227; ...
        0.985081948046 0.166407932485 0.043840114517; ...
        0.799395315780 -0.162097726540 0.578525242457; ...
        -0.420762355261 0.173028918374 0.890516722921; ...
        0.764559207832 0.336103454008 -0.549985168822; ...
        -0.760888489844 -0.344582260348 0.549828856897; ...
        0.084472616020 0.538920482568 0.838110428650; ...
        0.362971995096 -0.671052920431 0.646482257109; ...
        -0.077719061065 0.978066559969 -0.193249972363; ...
        -0.076072993014 -0.514704902253 -0.853985809795; ...
        0.273538375135 -0.959384592857 0.068977969797; ...
        0.301800326081 -0.036868539481 0.952658004729; ...
        -0.996891755320 -0.059015282365 -0.052192189278; ...
        -0.786920352434 0.187954609377 -0.587732442306; ...
        -0.659628679412 -0.424097649240 -0.620508814769; ...
        0.420802283612 -0.737952419329 -0.527590432925; ...
        0.806101843892 -0.582767614707 0.102867509555; ...
        -0.225406674727 0.112219296983 -0.967780274842; ...
        -0.737813919337 0.666892524247 -0.104331114898; ...
        0.615566050703 0.776673911335 -0.133626616635; ...
        -0.370512838694 0.665830974211 -0.647602771878; ...
        -0.287603733017 -0.874312423156 0.390975548430; ...
        -0.842393883403 0.314826265898 0.437329358156];
elseif nargin == 1
    DE_scheme = [zeros(nob0s,3);...
        -0.404485801990 0.801145160285 0.441086916763; ...
        -0.748942567232 -0.662363704264 -0.018956641627; ...
        -0.236315790505 -0.912424523381 -0.334120242265; ...
        0.693310892542 0.482245749539 0.535498873327; ...
        -0.203714420443 -0.442765397309 0.873189119177; ...
        0.399692451068 -0.069366022027 -0.914020951372; ...
        0.245878456697 0.627777017425 -0.738538963715; ...
        0.847226441419 -0.286749277133 -0.447193704142; ...
        0.234753426819 0.909610476648 0.342782160227; ...
        0.985081948046 0.166407932485 0.043840114517; ...
        0.799395315780 -0.162097726540 0.578525242457; ...
        -0.420762355261 0.173028918374 0.890516722921; ...
        0.764559207832 0.336103454008 -0.549985168822; ...
        -0.760888489844 -0.344582260348 0.549828856897; ...
        0.084472616020 0.538920482568 0.838110428650; ...
        0.362971995096 -0.671052920431 0.646482257109; ...
        -0.077719061065 0.978066559969 -0.193249972363; ...
        -0.076072993014 -0.514704902253 -0.853985809795; ...
        0.273538375135 -0.959384592857 0.068977969797; ...
        0.301800326081 -0.036868539481 0.952658004729; ...
        -0.996891755320 -0.059015282365 -0.052192189278; ...
        -0.786920352434 0.187954609377 -0.587732442306; ...
        -0.659628679412 -0.424097649240 -0.620508814769; ...
        0.420802283612 -0.737952419329 -0.527590432925; ...
        0.806101843892 -0.582767614707 0.102867509555; ...
        -0.225406674727 0.112219296983 -0.967780274842; ...
        -0.737813919337 0.666892524247 -0.104331114898; ...
        0.615566050703 0.776673911335 -0.133626616635; ...
        -0.370512838694 0.665830974211 -0.647602771878; ...
        -0.287603733017 -0.874312423156 0.390975548430; ...
        -0.842393883403 0.314826265898 0.437329358156];
else
    disp('Error: too many input arguments given');
end
DE_scheme = [DE_scheme(:,2) DE_scheme(:,1) -DE_scheme(:,3)];
