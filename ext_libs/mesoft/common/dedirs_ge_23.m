% function dedirs_ge_23 gives DE_scheme with b0s if amount of b0s is given
% SuS April 2008 (copied from original tensor.dat file)

function DE_scheme = dedirs_ge_23(nob0s)

if nargin == 0
    DE_scheme = [0.652755932306 -0.746414814390 -0.129516862607;...
        -0.571782213729 -0.244225346813 -0.783210750716;...
        -0.709688239042 -0.669143489929 -0.220430472601;...
        0.074260622701 0.149277887820 -0.986002774907;...
        -0.016875357639 -0.997959076091 -0.061586563074;...
        -0.009413501650 0.687898371154 0.725745972741;...
        0.772596564939 0.523362459372 0.359424935090;...
        0.263157839075 0.959418508514 0.101311782404;...
        0.400751973491 0.108379743006 0.909753640855;...
        0.839823682437 -0.276754944100 0.467014864147;...
        0.009413450197 -0.687898401663 -0.725745944490;...
        0.677548194000 0.538161884565 -0.501308518588;...
        0.609302219506 -0.266442386686 -0.746832819232;...
        -0.975664874445 0.152681067858 -0.157373899655;...
        0.030264912093 0.815154628228 -0.578452216846;...
        0.984938230299 -0.008800677927 -0.172682455866;...
        0.298116671784 -0.690842928425 0.658682395582;...
        -0.280144459258 -0.117059244474 0.952793899661;...
        -0.509780373863 0.859977207332 0.023731272475;...
        -0.684554448284 0.440828515136 0.580564748824;...
        -0.393141091044 -0.758368362464 0.519920676013;...
        -0.874060593742 -0.248142537028 0.417664171058;...
        -0.587825140081 0.476911936064 -0.653465079349];
elseif nargin == 1
    DE_scheme = [zeros(nob0s,3);...
        0.652755932306 -0.746414814390 -0.129516862607;...
        -0.571782213729 -0.244225346813 -0.783210750716;...
        -0.709688239042 -0.669143489929 -0.220430472601;...
        0.074260622701 0.149277887820 -0.986002774907;...
        -0.016875357639 -0.997959076091 -0.061586563074;...
        -0.009413501650 0.687898371154 0.725745972741;...
        0.772596564939 0.523362459372 0.359424935090;...
        0.263157839075 0.959418508514 0.101311782404;...
        0.400751973491 0.108379743006 0.909753640855;...
        0.839823682437 -0.276754944100 0.467014864147;...
        0.009413450197 -0.687898401663 -0.725745944490;...
        0.677548194000 0.538161884565 -0.501308518588;...
        0.609302219506 -0.266442386686 -0.746832819232;...
        -0.975664874445 0.152681067858 -0.157373899655;...
        0.030264912093 0.815154628228 -0.578452216846;...
        0.984938230299 -0.008800677927 -0.172682455866;...
        0.298116671784 -0.690842928425 0.658682395582;...
        -0.280144459258 -0.117059244474 0.952793899661;...
        -0.509780373863 0.859977207332 0.023731272475;...
        -0.684554448284 0.440828515136 0.580564748824;...
        -0.393141091044 -0.758368362464 0.519920676013;...
        -0.874060593742 -0.248142537028 0.417664171058;...
        -0.587825140081 0.476911936064 -0.653465079349];
elseif nargin > 1
    disp('Error: too many input arguments');
end
DE_scheme = [DE_scheme(:,2) DE_scheme(:,1) -DE_scheme(:,3)];

