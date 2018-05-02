import os, sys

from .external import run_bash

def rotate(image_in, image_out, angle=None, z=0.8, cube=False, clean=False, **kwargs):
    """Rotate a fits file using mRotate or mProject.

    If *angle* is given the function uses mRotate. If *angle=None* then the
    CROTA keywords are deleted from the header.

    By default the -X option of mProject is used, i.e. the rotated input image
    fits completely into the output image.

    Parameters:
        image_in (str): input FITS name.
        image_out (str): output FITS name.
        angle (float, default=None): rotation angle.
        z (float, default=0.8): drizzle algorithm scaling.
        cube (bool, default=False): if input image is a cube.
        clean (bool, default=False): delete the files created by this function.
        kwargs: other command line options for mProject. 
    """

    if angle:
        raise NotImplementedError
    else:
        # Write a header file
        header = mHdr(image_in)

        # Delete CROTA keys
        lines = ''
        with open(header, 'r') as inhead:
            for line in inhead:
                if line.startswith('CROTA'):
                    continue
                else:
                    lines += line
        
        # Save new header in the same file
        with open(header, 'w') as outhead:
            outhead.write(lines)

        if cube:
            cmd = 'mProjectCube -X -z %f %s %s %s'
        else:
            cmd = 'mProject -X -z %f %s %s %s'

        run_bash(cmd % (z, image_in, image_out, header))

def mHdr(image_in, header_file=None):
    """Create a text file with the header.

    If *header_file* is not given, the extension will be replaced by *.hdr*.

    Parameters:
        image_in (str): input image file name.
        header_file (str, default=None): header file name.
    """
    cmd = 'mGetHdr %s %s'
    if not header_file:
        ext = os.path.splitext(image_in)[1]
        header_file = image_in.replace(ext, '.hdr')
    else:
        pass

    run_bash(cmd % (image_in, header_file))
    return header_file

if __name__=='__main__':
    infile = 'Model_1000_PdBI1300D_raw.fits'
    outfile = 'Model_1000_PdBI1300D_projected.fits'
    rotate(infile, outfile, z=1.2, cube=True)
