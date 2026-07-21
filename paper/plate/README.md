# Figure plates — shared house style

`style.py` is the shared house style for the manuscript figures: fonts, ink
colour, panel-letter size, white-border trimming (`load_trim`) and panel
placement (`place_panel`, `fit_rect`). It keeps every figure visually
consistent.

It is imported by each figure's composer under `paper/figureN_rebuild/`, which
is where the actual layouts live — one self-contained pipeline per figure (see
`paper/README.md`). The scientific plots stay 100 % the R tool's own output;
the Python side only does the layout, deriving each figure's height from its
panels' real aspect ratios so cells match their images with no stretching and
no white gaps. All geometry is in inches (isotropic), converted to figure
fractions only at draw time.
