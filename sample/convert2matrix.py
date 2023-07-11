# about: converts an image to matrix data (cmatrix.h format)
# usage: python convert2matrix.py imagefile > matrix.dat

# requires PIL or Pillow library

import sys

from PIL import Image, ImageOps

image = Image.open(sys.argv[1])
image_gray = ImageOps.grayscale(image)

w, h = image.size
# prints "height width" to the 1st line
sys.stdout.write('{} {}\n'.format(h, w))
for y in range(h):
    for x in range(w):
        sys.stdout.write('{} '.format(image_gray.getpixel((x, y))))
    sys.stdout.write('\n')
