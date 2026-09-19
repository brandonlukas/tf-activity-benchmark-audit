"""Static geometry check for the inline SVGs: text outside the viewBox, text-text overlap,
text overflowing the rect it sits in. Widths are conservative over-estimates so misses are unlikely."""
import re, sys
from xml.etree import ElementTree as ET

SANS, MONO = 0.545, 0.600          # em advance; Plex Mono is exactly 0.6
def adv(cls):
    if cls and ('mono' in cls or 'big' in cls): return MONO
    return SANS
def size(cls):
    cls = cls or ''
    if 'big' in cls: return 17.0
    if 'mono' in cls: return 12.5
    if 'lbl' in cls: return 11.0
    return 12.0
def width(t, cls):
    fs, a = size(cls), adv(cls)
    extra = 0.04*fs if (cls and 'lbl' in cls) else 0.0   # letter-spacing on .lbl
    return len(t)*fs*a + max(0, len(t)-1)*extra

html = open('writeup.html').read()
svgs = re.findall(r'(<svg class="dgm".*?</svg>)', html, re.S)
print(f"{len(svgs)} svg fragments\n")
bad = 0
for n, frag in enumerate(svgs, 1):
    root = ET.fromstring(frag)
    vb = [float(v) for v in root.get('viewBox').split()]
    W, H = vb[2], vb[3]
    rects, texts = [], []
    def walk(el, cls_inherit=''):
        cls = (el.get('class') or '') + ' ' + cls_inherit
        tag = el.tag.split('}')[-1]
        if tag == 'rect' and 'box' in cls:
            rects.append((float(el.get('x')), float(el.get('y')),
                          float(el.get('width')), float(el.get('height'))))
        if tag == 'text':
            if el.get('transform'):
                return                                    # rotated label, checked by hand
            t = ''.join(el.itertext())
            x, y = float(el.get('x')), float(el.get('y'))
            w = width(t, cls)
            anc = el.get('text-anchor') or 'start'
            x0 = x - w/2 if anc == 'middle' else (x - w if anc == 'end' else x)
            fs = size(cls)
            texts.append({'t': t, 'x0': x0, 'x1': x0+w, 'y0': y-fs*0.80, 'y1': y+fs*0.22,
                          'yb': y, 'fs': fs})
        for ch in el: walk(ch, cls)
    walk(root)

    issues = []
    for tx in texts:
        if tx['x1'] > W + 0.5: issues.append(f"CROPPED right by {tx['x1']-W:5.1f}px: {tx['t']!r}")
        if tx['x0'] < -0.5:    issues.append(f"CROPPED left  by {-tx['x0']:5.1f}px: {tx['t']!r}")
        if tx['y1'] > H + 0.5: issues.append(f"CROPPED bottom by {tx['y1']-H:4.1f}px: {tx['t']!r}")
        if tx['y0'] < -0.5:    issues.append(f"CROPPED top   by {-tx['y0']:5.1f}px: {tx['t']!r}")
    for i in range(len(texts)):
        for j in range(i+1, len(texts)):
            a, b = texts[i], texts[j]
            ox = min(a['x1'], b['x1']) - max(a['x0'], b['x0'])
            oy = min(a['y1'], b['y1']) - max(a['y0'], b['y0'])
            if ox > 1.5 and oy > 1.5:
                issues.append(f"OVERLAP {ox:4.1f}x{oy:4.1f}px: {a['t']!r} / {b['t']!r}")
    for tx in texts:
        for (rx, ry, rw, rh) in rects:
            cx = (tx['x0']+tx['x1'])/2
            inside_v = ry - 2 <= tx['yb'] <= ry + rh + 2
            if inside_v and rx - 2 <= cx <= rx + rw + 2:       # label belongs to this box
                if tx['x0'] < rx - 2 or tx['x1'] > rx + rw + 2:
                    over = max(rx - tx['x0'], tx['x1'] - (rx+rw))
                    issues.append(f"SPILLS box by {over:5.1f}px: {tx['t']!r}")
                break
    print(f"--- figure {n}  viewBox {W:.0f}x{H:.0f}  {len(texts)} labels, {len(rects)} boxes")
    if issues:
        bad += len(issues)
        for m in sorted(set(issues)): print("   ", m)
    else:
        print("    clean")
print(f"\n{bad} issues")
