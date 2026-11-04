import sys, csv, pymupdf
doc = pymupdf.open(sys.argv[1]); pg = doc[5]; out = sys.argv[2]
dr = pg.get_drawings()
HW = 0.235                       # half of the 0.47 pt stroke width
X0 = (177.082 + 177.552) / 2     # month 0 tick centre
X42 = (403.342 + 403.812) / 2    # month 42 tick centre
Y1 = (79.299 + 79.769) / 2       # S = 1.0 tick centre
Y0 = (195.562 + 196.032) / 2     # S = 0.0 tick centre
xm = lambda m: X0 + (X42 - X0) * m / 42
s_of = lambda y: (Y0 - y) / (Y0 - Y1)
def col(d): return [round(x, 2) for x in (d.get("fill") or ())]
def find(c, n):
    return [d for d in dr if col(d) == c and len(d["items"]) >= n and d["rect"].width > 60][0]
def curve(d):
    # forward horizontal segments are the lower edge of the stroked step line
    seg = []
    for it in d["items"]:
        if it[0] == "l" and abs(it[1].y - it[2].y) < 1e-6 and it[2].x > it[1].x:
            seg.append((it[1].x, it[2].x, it[1].y - HW))
    return seg
def s_at(seg, x):
    hit = [s for s in seg if s[0] <= x <= s[1]]
    return s_of(max(h[2] for h in hit)) if hit else None
def band_at(d, x):
    ys = []
    for it in d["items"]:
        a, b = (it[1], it[2]) if it[0] == "l" else (it[1], it[4]) if it[0] == "c" else (None, None)
        if a is None: continue
        if min(a.x, b.x) <= x <= max(a.x, b.x):
            ys.append(a.y if a.x == b.x else a.y + (b.y - a.y) * (x - a.x) / (b.x - a.x))
    return (s_of(max(ys)), s_of(min(ys))) if ys else (None, None)
red = find([0.89, 0.0, 0.04], 400); blue = find([0.23, 0.28, 0.59], 2000)
pink = find([0.95, 0.6, 0.62], 300); lilac = find([0.63, 0.61, 0.8], 600)
sr, sb = curve(red), curve(blue)
print("red segments", len(sr), "blue segments", len(sb))
print("red path max x -> months", round((max(s[1] for s in sr) - X0) / (X42 - X0) * 42, 1), "S there", round(s_of(sorted(sr, key=lambda s: s[1])[-1][2]), 3))
atrisk = {0:587,6:549,12:475,18:401,24:317,30:232,36:145,42:65}
rows = []
print("month  treated  lo     hi     | untreated lo    hi")
# read just left of each tick so that a step exactly on the tick is not ambiguous
for m in range(0, 43, 6):
    x = xm(m) + (0.001 if m == 0 else -0.001)
    s = s_at(sr, x); u = s_at(sb, x)
    lo, hi = band_at(pink, x - 0.4) if m > 0 else (1, 1); ul, uh = band_at(lilac, x - 0.4) if m > 0 else (1, 1)
    print(f"{m:>4}  {s:.4f}  {lo:.4f} {hi:.4f} | {u:.4f} {ul:.4f} {uh:.4f}")
    rows.append(dict(month=m, s_treated=round(s,4), lo_treated=round(lo,4), hi_treated=round(hi,4), n_risk=atrisk[m], s_untreated=round(u,4), lo_untreated=round(ul,4), hi_untreated=round(uh,4)))
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
