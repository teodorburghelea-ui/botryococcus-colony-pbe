"""Rebuild the native 300-dpi scan of page 574 (PDF page 2) of Zhang & Kojima
(1998) from its embedded CCITT strips. Requires poppler-utils (pdfimages).
Usage: python3 extract_page.py /path/to/Zhang_Kojima_1998.pdf"""
import glob, subprocess, sys, tempfile, os
import numpy as np
from PIL import Image
pdf = sys.argv[1]
with tempfile.TemporaryDirectory() as t:
    subprocess.run(['pdfimages', '-f', '2', '-l', '2', '-png', pdf, os.path.join(t, 'p2')], check=True)
    strips = [np.array(Image.open(f).convert('L')) for f in sorted(glob.glob(os.path.join(t, 'p2-*.png')))]
Image.fromarray(np.vstack(strips)).save('page2_native.png')
print('page2_native.png', np.vstack(strips).shape)
