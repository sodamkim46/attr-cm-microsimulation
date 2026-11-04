Scripts that read the survival values in TRACE_validation_Brier_V.1.5.R from the published figures.
Python with pymupdf, numpy, pillow.
thaos_read_vector_figure.py <THAOS.pdf> <out.csv>             THAOS KM figure (page 6), vector paths.
masri_extract_image.py <Masri.pdf> <out_dir>                  Masri Figure 1 (page 5), saves the embedded image.
masri_read_pixels.py <masri_fig1_raw.png> <out.csv>           Reads the curve and CI band at 6-month ticks.
