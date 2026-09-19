"""JATS XML -> plain markdown: title, authors, abstract, body sections. Tables/figures reduced to captions."""
import sys, xml.etree.ElementTree as ET
root = ET.parse(sys.argv[1]).getroot()
def txt(e): return " ".join("".join(e.itertext()).split()) if e is not None else ""
front = root.find(".//front")
print("#", txt(root.find(".//article-title")))
authors = [txt(n.find("given-names")) + " " + txt(n.find("surname")) for n in root.findall(".//contrib[@contrib-type='author']/name")]
print("\n" + ", ".join(authors))
for aid in root.findall(".//article-id[@pub-id-type='doi']"): print("\nDOI:", txt(aid))
if (ab := root.find(".//abstract")) is not None: print("\n## Abstract\n\n" + txt(ab))
def walk(sec, depth):
    t = sec.find("title"); print(f"\n{'#' * min(depth, 6)} {txt(t)}\n") if t is not None else None
    for ch in sec:
        if ch.tag == "sec": walk(ch, depth + 1)
        elif ch.tag == "p": print(txt(ch) + "\n")
        elif ch.tag in ("fig", "table-wrap"): print(f"[{ch.tag}] {txt(ch.find('label'))} {txt(ch.find('caption'))}\n")
body = root.find(".//body")
if body is not None:
    for ch in body:
        if ch.tag == "sec": walk(ch, 2)
        elif ch.tag == "p": print(txt(ch) + "\n")
