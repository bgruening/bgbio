#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os, sys
import argparse
from PIL import Image
import operator

__doc__ = """
Example:
./color_pixel_counter.py -i red2.png --mark
File: red2.png - colored pixel: 16.003125 (256000/40968)
./color_pixel_counter.py -i red2.png --mark --upper 256 40 40 --lower 190 0 0
File: red2.png - colored pixel: 16.003125 (256000/40968)

"""

def main(options):

    for infile in options.infiles:
        img = Image.open( infile ).convert('RGB')

        upper = options.upper
        lower = options.lower

        # Get the size of the image in pixels
        width, height = img.size
        count_pixel = 0
        counter = 0
        # Process every pixel
        for x in range(width):
            for y in range(height):
                pixel = img.getpixel( (x,y) )
                if False not in map(operator.le, lower, pixel) and False not in map(operator.ge, upper, pixel):
                    count_pixel += 1
                    if options.mark:
                        img.putpixel( (x,y), (0,0,0))
                counter += 1

        print( "File: %s - colored pixel: %s (%s/%s)" % (os.path.basename( infile ), (count_pixel/(counter/100.0)), counter, count_pixel) )

        if options.mark:
            file_ext = os.path.splitext( infile )[-1]
            img.save( infile.replace(file_ext, '_marked%s' % file_ext) )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Pixel counting of a specific color range.',
        epilog = __doc__,
        formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--infiles', nargs='+', required=True, help='Image file/files.')
    parser.add_argument('-m', '--mark', action='store_true', default=False, help='Replace all counted pixels with black.')
    parser.add_argument('-u', '--upper', metavar=' 256 40 40', default=[256,40,40], type=int, nargs='+', help='Upper color border - default: (256,40,40)')
    parser.add_argument('-l', '--lower', metavar='190 0 0', default=[190,0,0], type=int, nargs='+', help='Lower color border - default: (190,0,0)')

    options = parser.parse_args()
    main(options)

