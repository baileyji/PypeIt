""" Module for LRIS specific codes
"""
from __future__ import absolute_import, division, print_function

import glob

import numpy as np
from astropy.io import fits

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import parse
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph

from pypeit import debugger

class KeckKCWISpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Keck/LRIS specific code
    """
    def __init__(self):
        # Get it started
        super(KeckKCWISpectrograph, self).__init__()
        self.spectrograph = 'keck_kcwi_base'
        self.telescope = telescopes.KeckTelescopePar()

    @staticmethod
    def default_pypeit_par():
        """
        Set default parameters for Keck LRISr reductions.
        """
        par = pypeitpar.PypeItPar()
        # Set wave tilts order
        par['calibrations']['slits']['sigdetect'] = 30.
        par['calibrations']['slits']['pcapar'] = [3,2,1,0]
        # 1D wavelengths
        par['calibrations']['wavelengths']['rms_threshold'] = 0.20  # Might be grism dependent
        # Always sky subtract, starting with default parameters
        par['scienceimage'] = pypeitpar.ScienceImagePar()
        # Always flux calibrate, starting with default parameters
        par['fluxcalib'] = pypeitpar.FluxCalibrationPar()
        # Always correct for flexure, starting with default parameters
        par['flexure'] = pypeitpar.FlexurePar()
        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 1]
        par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['calibrations']['pixelflatframe']['exprng'] = [None, 60]
        par['calibrations']['traceframe']['exprng'] = [None, 60]
        par['scienceframe']['exprng'] = [59, None]
        return par

    def header_keys(self):
        """
        Return a dictionary with the header keywords to read from the
        fits file.

        Returns:
            dict: A nested dictionary with the header keywords to read.
            The first level gives the extension to read and the second
            level gives the common name for header values that is passed
            on to the PypeItMetaData object.
        """
        hdr_keys = {}
        hdr_keys[0] = {}
        hdr_keys[1] = {}
        hdr_keys[2] = {}
        hdr_keys[3] = {}
        hdr_keys[4] = {}

        # Copied over defaults
        hdr_keys[0]['idname'] = 'IMTYPE'
        hdr_keys[0]['time'] = 'MJD'
        #hdr_keys[0]['date'] = 'DATE'
        hdr_keys[0]['utc'] = 'UTC'
        hdr_keys[0]['ut'] = 'UT'
        hdr_keys[0]['ra'] = 'RA'
        hdr_keys[0]['dec'] = 'DEC'
        hdr_keys[0]['airmass'] = 'AIRMASS'
        hdr_keys[0]['binning'] = 'BINNING'
        hdr_keys[0]['decker'] = 'SLITNAME'
        hdr_keys[0]['dichroic'] = 'DICHNAME'

        hdr_keys[0]['target'] = 'TARGNAME'
        hdr_keys[0]['exptime'] = 'TTIME'
        hdr_keys[0]['hatch'] = 'TRAPDOOR'
        hdr_keys[0]['dispname'] = 'GRANAME'
        hdr_keys[0]['dispangle'] = 'GRANGLE'
        hdr_keys[0]['wavecen'] = 'WAVELEN'
        hdr_keys[0]['spectrograph'] = 'INSTRUME'
        #hdr_keys[1]['CCDGEOM'] = 'CCDGEOM'
        #hdr_keys[1]['CCDNAME01'] = 'CCDNAME'
        #hdr_keys[3]['CCDNAME02'] = 'CCDNAME'

        #lamp_names = ['MERCURY', 'NEON', 'ARGON', 'CADMIUM', 'ZINC', 'KRYPTON', 'XENON',
        #              'FEARGON', 'DEUTERI', 'FLAMP1', 'FLAMP2', 'HALOGEN']
        #for kk,lamp_name in enumerate(lamp_names):
        #    hdr_keys[0]['lampstat{:02d}'.format(kk+1)] = lamp_name

        return hdr_keys

    def metadata_keys(self):
        return super(KeckKCWISpectrograph, self).metadata_keys() \
                    + ['binning', 'dichroic', 'dispangle']

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        if ftype == 'science':
            return good_exp & (fitstbl['idname'] == 'OBJECT')
        if ftype == 'arc':
            return good_exp & (fitstbl['idname'] == 'ARCLAMP')
        if ftype == 'trace':
            return good_exp & (fitstbl['idname'] == 'CONTBARS')
        if ftype == 'pixelflat':
            return good_exp & (fitstbl['idname'] == 'FLATLAMP')
        if ftype == 'bias':
            return good_exp & (fitstbl['idname'] == 'BIAS')

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)
  
    def lamps(self, fitstbl, status):
        """
        Check the lamp status.

        Args:
            fitstbl (:obj:`astropy.table.Table`):
                The table with the fits header meta data.
            status (:obj:`str`):
                The status to check.  Can be `off`, `arcs`, or `dome`.
        
        Returns:
            numpy.ndarray: A boolean array selecting fits files that
            meet the selected lamp status.

        Raises:
            ValueError:
                Raised if the status is not one of the valid options.
        """
        pass
        '''
        if status == 'off':
            # Check if all are off
            return np.all(np.array([ (fitstbl[k] == 'off') | (fitstbl[k] == 'None')
                                        for k in fitstbl.keys() if 'lampstat' in k]), axis=0)
        if status == 'arcs':
            # Check if any arc lamps are on
            arc_lamp_stat = [ 'lampstat{0:02d}'.format(i) for i in range(1,9) ]
            return np.any(np.array([ fitstbl[k] == 'on' for k in fitstbl.keys()
                                            if k in arc_lamp_stat]), axis=0)
        if status == 'dome':
            # Check if any dome lamps are on
            dome_lamp_stat = [ 'lampstat{0:02d}'.format(i) for i in range(9,13) ]
            return np.any(np.array([ fitstbl[k] == 'on' for k in fitstbl.keys()
                                            if k in dome_lamp_stat]), axis=0)
        raise ValueError('No implementation for status = {0}'.format(status))
        '''

    def load_raw_img_head(self, raw_file, det=None, **null_kwargs):
        """
        Wrapper to the raw image reader for LRIS

        Args:
            raw_file:  str, filename
            det: int, REQUIRED
              Desired detector
            **null_kwargs:
              Captured and never used

        Returns:
            raw_img: ndarray
              Raw image;  likely unsigned int
            head0: Header

        """
        raw_img, head0, _, _ = read_kcwi(raw_file)

        return raw_img, head0

    def get_image_section(self, filename, det, section='datasec'):
        """
        Return a string representation of a slice defining a section of
        the detector image.

        Overwrites base class function to use :func:`read_kcwi` to get
        the image sections.

        .. todo::
            - It feels really ineffiecient to just get the image section
              using the full :func:`read_kcwi`.  Can we parse that
              function into something that can give you the image
              section directly?

        This is done separately for the data section and the overscan
        section in case one is defined as a header keyword and the other
        is defined directly.
        
        Args:
            filename (str):
                data filename
            det (int):
                Detector number
            section (:obj:`str`, optional):
                The section to return.  Should be either datasec or
                oscansec, according to the :class:`DetectorPar`
                keywords.

        Returns:
            list, bool: A list of string representations for the image
            sections, one string per amplifier, followed by three
            booleans: if the slices are one indexed, if the slices
            should include the last pixel, and if the slice should have
            their order transposed.
        """
        # Read the file
        temp, head0, secs, direc = read_kcwi(filename)
        if section == 'datasec':
            return secs[0], False, False, False
        elif section == 'oscansec':
            return secs[1], False, False, False
        else:
            raise ValueError('Unrecognized keyword: {0}'.format(section))

    def get_image_shape(self, filename=None, det=None, **null_kwargs):
        """
        Overrides :class:`Spectrograph.get_image_shape` for LRIS images.

        Must always provide a file.
        """
        # Cannot be determined without file
        if filename is None:
            raise ValueError('Must provide a file to determine the shape of an LRIS image.')

        # Use a file
        self._check_detector()
        self.naxis = (self.load_raw_frame(filename, det=det)[0]).shape
        return self.naxis

    def get_match_criteria(self):
        match_criteria = {}
        for key in framematch.FrameTypeBitMask().keys():
            match_criteria[key] = {}
        #
        match_criteria['standard']['match'] = {}
        match_criteria['standard']['match']['dispname'] = ''
        match_criteria['standard']['match']['dichroic'] = ''
        match_criteria['standard']['match']['binning'] = ''
        match_criteria['standard']['match']['decker'] = ''
        # Bias
        match_criteria['bias']['match'] = {}
        match_criteria['bias']['match']['binning'] = ''
        # Pixelflat
        match_criteria['pixelflat']['match'] = match_criteria['standard']['match'].copy()
        # Traceflat
        match_criteria['trace']['match'] = match_criteria['standard']['match'].copy()
        # Arc
        match_criteria['arc']['match'] = match_criteria['standard']['match'].copy()

        # Return
        return match_criteria


class KeckKCWIBSpectrograph(KeckKCWISpectrograph):
    """
    Child to handle Keck/KCWIb specific code
    """
    def __init__(self):
        # Get it started
        super(KeckKCWIBSpectrograph, self).__init__()
        self.spectrograph = 'keck_kcwi_blue'
        self.camera = 'LRISb'
        self.detector = [
                # Detector 1
                pypeitpar.DetectorPar(
                            dataext         = 1,
                            dispaxis        = 0,
                            dispflip        = False,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.135,
                            darkcurr        = 0.0,
                            saturation      = 65535.,
                            nonlinear       = 0.86,
                            numamplifiers   = 4,
                            gain            = [0.172, 0.173, 0.168, 0.174],
                            ronoise         = [4.0]*4,
                            suffix          = '_01blue'
                            ),
            ]
        self.numhead = 1
        # Uses default timeunit
        # Uses default primary_hdrext
        self.sky_file = 'sky_LRISb_600.fits'

    @staticmethod
    def default_pypeit_par():
        """
        Set default parameters for Keck LRISr reductions.
        """
        par = KeckKCWISpectrograph.default_pypeit_par()
        par['rdx']['spectrograph'] = 'keck_lris_blue'
        return par

    def check_headers(self, headers):
        """
        Check headers match expectations for an LRISb exposure.

        See also
        :func:`pypeit.spectrographs.spectrograph.Spectrograph.check_headers`.

        Args:
            headers (list):
                A list of headers read from a fits file
        """
        expected_values = { '0.INSTRUME': 'KCWI', '0.CAMERA': 'BLUE'}
        super(KeckKCWIBSpectrograph, self).check_headers(headers, expected_values=expected_values)

    def header_keys(self):
        hdr_keys = super(KeckKCWIBSpectrograph, self).header_keys()
        return hdr_keys

    def bpm(self, filename=None, det=None, **null_kwargs):
        """ Generate a BPM

        Parameters
        ----------
        det : int, REQUIRED
        **null_kwargs:
           Captured and never used

        Returns
        -------
        badpix : ndarray

        """
        # Get the empty bpm: force is always True
        self.empty_bpm(filename=filename, det=det)

        return self.bpm_img


'''
def bpm(slf, camera, fitsdict, det):
    """  Wrapper for core_bpm
    Will likely be deprecated

    Parameters
    ----------
    slf
    camera
    fitsdict
    det

    Returns
    -------
    badpix : ndarray

    """
    sidx = slf._idx_sci[0]
    # Binning
    xbin, ybin = [int(ii) for ii in fitsdict['binning'][sidx].split(',')]
    return core_bpm(xbin, ybin, camera, det)
'''

def read_kcwi(raw_file):
    # Read
    msgs.info("Reading LRIS file: {:s}".format(raw_file))
    hdu = fits.open(raw_file)
    head0 = hdu[0].header

    # Grab the sections
    ibsec, idsec, tsec, direc = map_ccd(head0)

    # Use em
    data = hdu[0].data.astype(float)

    # Turn dsec into IRAF format but with Python indexing
    bsec, dsec = [], []
    for kk, iidsec in enumerate(idsec):
        dsec.append('[{}:{},{}:{}]'.format(iidsec[0], iidsec[1], iidsec[2], iidsec[3]))
        bsec.append('[{}:{},{}:{}]'.format(ibsec[kk][0], ibsec[kk][1], ibsec[kk][2], ibsec[kk][3]))

    # Return
    return data, head0, (dsec, bsec), direc


def map_ccd(header):
    """Return CCD section variables useful for processing
    Uses FITS keyword NVIDINP to determine how many amplifiers were used
    to read out the CCD.  Then reads the corresponding BSECn, and
    DSECn keywords, where n is the amplifier number.  The indices are
    converted to Python (0-biased, y axis first) indices and an array
    is constructed for each of the two useful sections of the CCD as
    follows:
    Bsec[0][0] - First amp, y lower limit
    Bsec[0][1] - First amp, y upper limit
    Bsec[0][2] - First amp, x lower limit
    Bsec[0][3] - First amp, x upper limit
    Bsec[1][0] - Second amp, y lower limit
    etc.
    Bsec is the full overscan region for the given amplifier and is used
    to calculate and perform the overscan subtraction.
    Dsec is the full CCD region for the given amplifier and is used to
    trim the image after overscan subtraction has been performed.
    Tsec accounts for trimming the image according to Dsec.
    Amps are assumed to be organized as follows:
    (0,ny)	--------- (nx,ny)
            | 3 | 4 |
            ---------
            | 1 | 2 |
    (0,0)	--------- (nx, 0)
    Args:
    -----
        self: instance of CcdPrimitive class (automatic)
    Returns:
    --------
        list: (int) y0, y1, x0, x1 for bias section
        list: (int) y0, y1, x0, x1 for data section
        list: (int) y0, y1, x0, x1 for trimmed section
        list: (bool) y-direction, x-direction, True if forward, else False
    """

    namps = header['NVIDINP']
    # TODO: check namps
    # section lists
    bsec = []
    dsec = []
    tsec = []
    direc = []
    # loop over amps
    for i in range(namps):
        sec, rfor = parse_imsec(header, section_key='BSEC%d' % (i+1))
        bsec.append(sec)
        sec, rfor = parse_imsec(header, section_key='DSEC%d' % (i+1))
        dsec.append(sec)
        direc.append(rfor)
        if i == 0:
            y0 = 0
            y1 = sec[1] - sec[0]
            x0 = 0
            x1 = sec[3] - sec[2]
        elif i == 1:
            y0 = 0
            y1 = sec[1] - sec[0]
            x0 = tsec[0][3] + 1
            x1 = x0 + sec[3] - sec[2]
        elif i == 2:
            y0 = tsec[0][1] + 1
            y1 = y0 + sec[1] - sec[0]
            x0 = 0
            x1 = sec[3] - sec[2]
        elif i == 3:
            y0 = tsec[0][1] + 1
            y1 = y0 + sec[1] - sec[0]
            x0 = tsec[0][3] + 1
            x1 = x0 + sec[3] - sec[2]
        else:
            # should not get here
            y0 = -1
            y1 = -1
            x0 = -1
            x1 = -1
            #self.log.info("ERROR - bad amp number: %d" % i)
        tsec.append((y0, y1, x0, x1))

    return bsec, dsec, tsec, direc


def parse_imsec(header, section_key=None):
    if section_key is None:
        return None, None
    else:
        # forward read?
        xfor = True
        yfor = True
        section = header[section_key]
        p1 = int(section[1:-1].split(',')[0].split(':')[0])
        p2 = int(section[1:-1].split(',')[0].split(':')[1])
        p3 = int(section[1:-1].split(',')[1].split(':')[0])
        p4 = int(section[1:-1].split(',')[1].split(':')[1])
        # tests for individual axes
        if p1 > p2:
            x0 = p2 - 1
            x1 = p1 - 1
            xfor = False
        else:
            x0 = p1 - 1
            x1 = p2 - 1
        if p3 > p4:
            y0 = p4 - 1
            y1 = p3 - 1
            yfor = False
        else:
            y0 = p3 - 1
            y1 = p4 - 1
        # package output
        sec = (y0, y1, x0, x1)
        rfor = (yfor, xfor)
        # use python axis ordering
        return sec, rfor