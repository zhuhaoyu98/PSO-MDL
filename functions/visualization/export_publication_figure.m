function export_publication_figure(figureHandle, outputBase, resolutionDpi)
%EXPORT_PUBLICATION_FIGURE Export a raster PNG and vector PDF of one figure.
%   OUTPUTBASE is an extension-free path. RESOLUTIONDPI applies to PNG only.
%   The caller controls figure geometry, fonts, and overwrite policy.
exportgraphics(figureHandle, [outputBase '.png'], ...
    'Resolution', resolutionDpi);
exportgraphics(figureHandle, [outputBase '.pdf'], ...
    'ContentType', 'vector');
end
