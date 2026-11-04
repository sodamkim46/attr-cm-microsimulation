import sys, io, pymupdf, numpy as np
from PIL import Image
d = pymupdf.open(sys.argv[1]); pg = d[4]; out = sys.argv[2]
xr = pg.get_images(full=True)[0][0]; b = d.extract_image(xr); im = Image.open(io.BytesIO(b["image"])).convert("RGB"); im.save(out + "/masri_fig1_raw.png")
a = np.asarray(im).astype(int); H, W, _ = a.shape; print("size", W, H, b["ext"])
R, G, B = a[..., 0], a[..., 1], a[..., 2]
line = (R > 190) & (G < 170) & (B < 170) & (R - G > 50)
print("line px", line.sum(), "x range", np.where(line.any(0))[0][[0, -1]], "y range", np.where(line.any(1))[0][[0, -1]])
dark = (R < 90) & (G < 90) & (B < 90)
cols = dark.sum(0); rows = dark.sum(1)
print("long dark cols", [(i, int(c)) for i, c in enumerate(cols) if c > H * 0.3][:10])
print("long dark rows", [(i, int(c)) for i, c in enumerate(rows) if c > W * 0.3][:10])
gray = (abs(R - G) < 8) & (abs(G - B) < 8) & (R > 215) & (R < 245)
print("band px", gray.sum())
ax_row = max(range(H), key=lambda i: rows[i]); ax_col = max(range(W), key=lambda j: cols[j]); print("x-axis row", ax_row, "y-axis col", ax_col)
# ticks: short dark runs just below the x axis and just left of the y axis
below = dark[ax_row + 3: ax_row + 8, :].all(0); xt = np.where(below)[0]; print("x tick px", xt.tolist()[:80])
left = dark[:, ax_col - 7: ax_col - 2].all(1); yt = np.where(left)[0]; print("y tick px", yt.tolist()[:40])
