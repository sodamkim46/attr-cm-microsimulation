Scripts that read the survival values in TRACE_validation_Brier_V.1.5.R from the published figures.
Python with pymupdf, numpy, pillow.
thaos_read_vector_figure.py <THAOS.pdf> <out.csv>             THAOS KM figure (page 6), vector paths.
masri_extract_image.py <Masri.pdf> <out_dir>                  Masri Figure 1 (page 5), saves the embedded image.
masri_read_pixels.py <masri_fig1_raw.png> <out.csv>           Reads the curve and CI band at 6-month ticks.
trials_read_figures.py <ATTR_ACT.pdf> <NEJM_Vutrisiran.pdf> <Long-term_ATTR-CM.pdf>   Trial curves (vector paths; raster for the ATTRibute-CM extension).
  ATTRibute-CM months 6 and 18: a step falls on the tick; the value before or after the step was chosen with the cumulative events in the figure table.
Fumagalli 2026 values and CIs are printed in its Figure 2; nothing is extracted from the curve.
