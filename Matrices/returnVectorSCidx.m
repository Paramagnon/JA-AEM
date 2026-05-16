function [SCidx, SC_tape] = returnVectorSCidx(Geom)
    SCidx = [];
    SC_tape = [];

    Noffset = 0;
    for ii = 1:length(Geom)
        for nn = 1:length(Geom(ii).SC)
            Nsc = length(Geom(ii).SC(nn).idx);
            SCidx(end+1:end+Nsc) = Geom(ii).SC(nn).idx + Noffset;
            SC_tape(end+1).idx = Geom(ii).SC(nn).idx + Noffset;
        end
        Noffset = Noffset + Geom(ii).Nbnd;
    end
end