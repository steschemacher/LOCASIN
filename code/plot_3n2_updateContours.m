function plot_3n2_updateContours(hContour)
    % Update the text label colors
    drawnow  % very important!
    levels = hContour.LevelList;
    labels = hContour.TextPrims;  % undocumented/unsupported
    lines  = hContour.EdgePrims;  % undocumented/unsupported
    for idx = 1 : numel(labels)
        labelValue = str2double(labels(idx).String);
        lineIdx = find(abs(levels-labelValue)<10*eps, 1);  % avoid FP errors using eps
        labels(idx).ColorData = lines(lineIdx).ColorData;  % update the label color
%         labels(idx).Font.Size = 8;                        % update the label font size
    end
    drawnow  % optional
end