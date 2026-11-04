# Reads the survival values for the trial curves in TRACE_validation_Brier_V.1.5.R from the published figures.
# ATTR-ACT Figure 2B and HELIOS-B Figure 2A are vector paths; ATTRibute-CM extension Figure 1 is a raster image.
# usage: trials_read_figures.py <ATTR_ACT.pdf> <NEJM_Vutrisiran.pdf> <Long-term_ATTR-CM.pdf>
import io, sys
import numpy as np
import pymupdf
from PIL import Image

path_attract, path_helios, path_ole = sys.argv[1:4]

def near(color, target, tol=0.02):
    return color is not None and len(color) == 3 and all(abs(a - b) < tol for a, b in zip(color, target))

# ATTR-ACT, page 6: tick marks x 184.7 = month 0, 487.2 = month 33; y 169.1 = 1.0, 319.5 = 0.0
page = pymupdf.open(path_attract)[5]
seg = [(i[1].x, i[1].y, i[2].x, i[2].y) for p in page.get_drawings()
       if near(p.get("color"), (.35, .53, .78)) and round(p["width"], 2) == .75 for i in p["items"] if i[0] == "l"]
def attract(month):
    x = 184.7 + month * (487.2 - 184.7) / 33
    y = max(max(a[1], a[3]) for a in seg if max(a[0], a[2]) <= x or a[0] <= x <= a[2])
    return 1 - (y - 169.1) / (319.5 - 169.1)
print("ATTR-ACT pooled tafamidis", {m: round(attract(m), 4) for m in range(6, 31, 6)})

# HELIOS-B, page 10 (rotated, so survival runs along x and time along y): x 54.86 = 100%, 113.23 = 70%; y 621.47 = month 0, 401.70 = month 42
page = pymupdf.open(path_helios)[9]
gold = [p for p in page.get_drawings() if near(p.get("color"), (.89, .62, .25)) and len(p["items"]) > 50 and p["rect"].y0 > 300][0]
def helios(month):
    y_month = 621.47 - month * (621.47 - 401.70) / 42
    xs = [x for i in gold["items"] if i[0] == "l" for x, y in ((i[1].x, i[1].y), (i[2].x, i[2].y)) if y >= y_month]
    return 1 - 0.3 * (max(xs) - 54.86) / (113.23 - 54.86)
print("HELIOS-B vutrisiran, overall", {m: round(helios(m), 4) for m in range(6, 43, 6)})

# ATTRibute-CM extension, page 5: largest embedded image; rows 29 = 1.0, 968 = 0.3; x ticks every 3 months
doc = pymupdf.open(path_ole)
xref = max(doc[4].get_images(full=True), key=lambda i: i[2] * i[3])[0]
img = np.asarray(Image.open(io.BytesIO(doc.extract_image(xref)["image"])).convert("RGB")).astype(int)
print("image", img.shape)
navy = (img[..., 0] < 90) & (img[..., 1] < 110) & (img[..., 2] > 110) & (img[..., 2] < 190)
navy[880:, :] = False  # legend
ticks = dict(zip(range(6, 43, 6), [386.5, 547, 707.5, 868, 1029, 1190, 1350.5]))
for month, x in ticks.items():
    for col in (int(round(x)) - 3, int(round(x)) - 1):
        rows = np.where(navy[:, col])[0]
        rows = rows[rows > 10]
        print(month, "col", col, "rows", rows.min(), rows.max(), "S", round(1 - .7 * ((rows.min() + rows.max()) / 2 - 29) / (968 - 29), 4))
