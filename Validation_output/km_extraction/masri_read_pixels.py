import sys, csv, numpy as np
from PIL import Image
a = np.asarray(Image.open(sys.argv[1]).convert("RGB")).astype(int); out = sys.argv[2]
R, G, B = a[..., 0], a[..., 1], a[..., 2]; lum = a.mean(2)
X0, X65 = 156.0, 1665.5; Y100, Y0 = 27.0, 527.5
xm = lambda m: X0 + (X65 - X0) * m / 65
s_of = lambda y: 100 * (Y0 - y) / (Y0 - Y100)
line = (R - G > 45) & (R - B > 45)
gray = (abs(R - G) < 10) & (abs(G - B) < 10) & (lum > 205) & (lum < 246)
def curve(x):
    c = int(round(x)); ys = np.where(line[:540, c])[0]
    ys = ys[ys > 15]
    return (float(ys.mean()), int(ys.max() - ys.min() + 1)) if len(ys) else (None, 0)
def band(x):
    # walk outward from the curve, skipping JPEG gaps of up to 6 pixels; average 5 columns to the left
    tops, bots = [], []
    for c in range(int(round(x)) - 4, int(round(x)) + 1):
        yc = curve(c)[0]
        if yc is None: continue
        y0 = int(round(yc)); on = gray[:, c] | line[:, c]
        t = y0
        while t > 15 and on[t - 7:t].any(): t -= 1
        b = y0
        while b < 526 and on[b + 1:b + 8].any(): b += 1
        tops.append(t); bots.append(b)
    return (float(np.median(tops)), float(np.median(bots)))
print("px per month", round((X65 - X0) / 65, 2), "| px per pct", round((Y0 - Y100) / 100, 2))
risk_m = list(range(0, 70, 5)); risk = [624, 597, 567, 533, 500, 470, 446, 412, 365, 293, 229, 119, 78, 7]
rows = []
for m in [0] + list(range(6, 61, 6)) + [64, 65]:
    x = xm(m) - (2 if m else -3)      # just left of the tick; just right at month 0
    y, h = curve(x); t, b = band(x)
    s = s_of(y); up, lo = s_of(t), s_of(b); n = float(np.interp(m, risk_m, risk))
    peto = s * np.sqrt(max(1 - s / 100, 0) / n) ; 
    print(f"{m:>3} S={s:6.2f} (line px {h}) band {lo:6.2f}-{up:6.2f} n~{n:5.0f} peto95 +-{1.96*peto:4.2f} band half {0.5*(up-lo):4.2f}")
    rows.append(dict(month=m, s=round(s / 100, 4), lower=round(lo / 100, 4), upper=round(up / 100, 4), n_risk=round(n)))
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
print("Paper: 65 months 54.1 (47.4-60.4)")
