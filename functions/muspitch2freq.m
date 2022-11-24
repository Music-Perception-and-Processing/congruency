function f = muspitch2freq(notes)

noteTable = { 'C'   'C'  'C';
              'C#'  'Db' 'CIS';
              'D'   'D'  'D';
              'D#'  'Eb' 'DIS';
              'E'   'E'  'E';
              'F'   'F'  'F';
              'F#'  'Gb' 'FIS';
              'G'   'G'  'G';
              'G#'  'Ab' 'GIS';
              'A'   'A'  'A';
              'A#'  'Bb' 'AIS';
              'B'   'B'  'H';
              ''    ''  '';};

if ~ischar(notes)
    [nCol, nRow] = size(notes);
else
    notes = {notes};
    nCol = 1;
    nRow = 1;
end

f = zeros(nCol,nRow);

for ii = 1:nCol

    for jj = 1:nRow

        note = char(notes(ii,jj));
        if isnan(str2double(note(end-1)))
            L = 0;
        else
            L = 1;
        end
        oct = str2num(note(end-L:end));
        notename = note(1:end-(L+1));
          
        chroma = find(max(strcmp(notename, noteTable), [], 2)); 
        octave = oct - 4;
        T = chroma + 12*octave - 1; 
        fC4 = 261.625565300599;  % Middle C (C4) is 261.63 Hz
        f(ii,jj) = fC4 .* 2 .^ (T / 12);

    end

end

end