""" Module for LRIS specific codes
"""
from __future__ import absolute_import, division, print_function

try:
    basestring
except NameError:  # For Python 3
    basestring = str

import glob

import numpy as np
from astropy.io import fits

from pypit import msgs
from pypit import arparse
from ..par.pypitpar import DetectorPar
from . import spectrograph

from pypit import ardebug as debugger

class KeckLRISSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Keck/LRIS specific code
    """

    # TODO: This should be an abstract class!

    def __init__(self):
        # Get it started
        super(KeckLRISSpectrograph, self).__init__()
        self.spectrograph = 'keck_lris'

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
        raw_img, head0, _ = read_lris(raw_file, det=det)

        return raw_img, head0

    def get_datasec(self, filename, det):
        """
        Load up the datasec and oscansec and also naxis0 and naxis1

        Args:
            filename (str):
                data filename
            det (int):
                Detector number

        Returns:
            datasec: list
            oscansec: list
            naxis0: int
            naxis1: int
        """
        # Check the detector
        if self.detector is None:
            raise ValueError('Must first define spectrograph detector parameters!')
        for d in self.detector:
            if not isinstance(d, DetectorPar):
                raise TypeError('Detectors must be specified using a DetectorPar instance.')
        
        # Read the file
        temp, head0, secs = read_lris(filename, det)

        # Get the data and overscan regions
        datasec, oscansec = [], []
        for kk in range(self.detector[det]['numamplifiers']):
            datasec.append(arparse.load_sections(secs[0][kk], fmt_iraf=False))
            oscansec.append(arparse.load_sections(secs[1][kk], fmt_iraf=False))

        # Return the sections and the shape of the image
        return (datasec, oscansec) + temp.shape


class KeckLRISBSpectrograph(KeckLRISSpectrograph):
    """
    Child to handle Keck/LRISb specific code
    """
    def __init__(self):
        # Get it started
        super(KeckLRISBSpectrograph, self).__init__()
        self.spectrograph = 'keck_lris_blue'
        self.detector = [
                # Detector 1
                DetectorPar(dataext=1, dispaxis=0, xgap=0., ygap=0., ysize=1., platescale=0.135,
                            darkcurr=0.0, saturation=65535., nonlinear=0.86, numamplifiers=2,
                            gain=[1.55, 1.56], ronoise=[3.9, 4.2], suffix='_01blue'),
                #Detector 2
                DetectorPar(dataext=2, dispaxis=0, xgap=0., ygap=0., ysize=1., platescale=0.135,
                            darkcurr=0., saturation=65535., nonlinear=0.86, numamplifiers=2,
                            gain=[1.63, 1.70], ronoise=[3.6, 3.6], suffix='_02blue')
            ]

    def setup_arcparam(self, arcparam, disperser=None, **null_kwargs):
        """
        Setup the arc parameters

        Args:
            arcparam: dict
            disperser: str, REQUIRED
            **null_kwargs:
              Captured and never used

        Returns:
            arcparam is modified in place

        """
        arcparam['lamps'] = ['NeI', 'ArI', 'CdI', 'KrI', 'XeI', 'ZnI', 'CdI', 'HgI']
        if disperser == '600/4000':
            arcparam['n_first']=2 # Too much curvature for 1st order
            arcparam['disp']=0.63 # Ang per pixel (unbinned)
            arcparam['b1']= 4.54698031e-04
            arcparam['b2']= -6.86414978e-09
            arcparam['wvmnx'][1] = 6000.
            arcparam['wv_cen'] = 4000.
        elif disperser == '400/3400':
            arcparam['n_first']=2 # Too much curvature for 1st order
            arcparam['disp']=1.02
            arcparam['b1']= 2.72694493e-04
            arcparam['b2']= -5.30717321e-09
            arcparam['wvmnx'][1] = 6000.
        elif disperser == '300/5000':
            arcparam['n_first'] = 2
            arcparam['wv_cen'] = 4500.
            arcparam['disp'] = 1.43
        else:
            msgs.error('Not ready for this disperser {:s}!'.format(disperser))


class KeckLRISRSpectrograph(KeckLRISSpectrograph):
    """
    Child to handle Keck/LRISr specific code
    """
    def __init__(self):
        # Get it started
        super(KeckLRISRSpectrograph, self).__init__()
        self.spectrograph = 'keck_lris_red'
        self.detector = [
                # Detector 1
                DetectorPar(dataext=1, dispaxis=0, xgap=0., ygap=0., ysize=1., platescale=0.135,
                            darkcurr=0.0, saturation=65535., nonlinear=0.76, numamplifiers=2,
                            gain=[1.255, 1.18], ronoise=[4.64, 4.76], suffix='_01red'),
                #Detector 2
                DetectorPar(dataext=2, dispaxis=0, xgap=0., ygap=0., ysize=1., platescale=0.135,
                            darkcurr=0., saturation=65535., nonlinear=0.76, numamplifiers=2,
                            gain=[1.191, 1.162], ronoise=[4.54, 4.62], suffix='_02red')
            ]

    # TODO: Anything that isn't general to the bpm methods for *all*
    # spectrograph should be held as part of the class (like detector)A
    # I think this means that bpm should be created when the data is
    # read using the binning from the fits headers.
    def bpm(self, binning=None, det=None, **null_kwargs):
        """ Generate a BPM

        Parameters
        ----------
        binning : str, REQUIRED
          Formatted like '1,1'
        det : int, REQUIRED
        **null_kwargs:
           Captured and never used

        Returns
        -------
        badpix : ndarray

        """
        xbin, ybin = [int(ii) for ii in binning.split(',')]
        xshp = 2048 // xbin
        yshp = 4096 // ybin
        badpix = np.zeros((yshp, xshp), dtype=np.int)
        # Do it
        if det == 2:
            msgs.info("Using hard-coded BPM for det=2 on LRISr")
            badc = 16 // xbin
            badpix[:, 0:badc] = 1.
        # Return
        return badpix

    def setup_arcparam(self, arcparam, disperser=None, fitstbl=None, arc_idx=None,
                       msarc_shape=None, binspectral=None, **null_kwargs):
        """
        Setup the arc parameters

        Args:
            arcparam: dict
            disperser: str, REQUIRED

        Returns:
            arcparam is modified in place

        """
        arcparam['wv_cen'] = fitstbl['wavecen'][arc_idx]
        arcparam['lamps'] = ['ArI','NeI','HgI','KrI','XeI']  # Should set according to the lamps that were on
        if disperser == '600/7500':
            arcparam['n_first']=3 # Too much curvature for 1st order
            arcparam['disp']=0.80 # Ang per pixel (unbinned)
            arcparam['b1']= 1./arcparam['disp']/msarc_shape[0] / binspectral
            arcparam['wvmnx'][1] = 11000.
        elif disperser == '600/10000':
            arcparam['n_first']=2 # Too much curvature for 1st order
            arcparam['disp']=0.80 # Ang per pixel (unbinned)
            arcparam['b1']= 1./arcparam['disp']/msarc_shape[0] / binspectral
            arcparam['wvmnx'][1] = 12000.
        elif disperser == '400/8500':
            arcparam['n_first']=2 # Too much curvature for 1st order
            arcparam['disp']=1.19 # Ang per pixel (unbinned)
            arcparam['b1']= 1./arcparam['disp']/msarc_shape[0] / binspectral
            arcparam['wvmnx'][1] = 11000.
            arcparam['min_ampl'] = 3000.  # Lines tend to be very strong
            arcparam['nsig_rej_final'] = 5.
        elif disperser == '900/5500':
            arcparam['n_first']=2 # Too much curvature for 1st order
            arcparam['disp']=0.53 # Ang per pixel (unbinned)
            arcparam['b1']= 1./arcparam['disp']/msarc_shape[0] / binspectral
            arcparam['wvmnx'][1] = 7000.
        else:
            msgs.error('Not ready for this disperser {:s}!'.format(disperser))


def read_lris(raw_file, det=None, TRIM=False):
    """
    Read a raw LRIS data frame (one or more detectors)
    Packed in a multi-extension HDU
    Based on readmhdufits.pro

    Parameters
    ----------
    raw_file : str
      Filename
    det : int, optional
      Detector number; Default = both
    TRIM : bool, optional
      Trim the image?

    Returns
    -------
    array : ndarray
      Combined image 
    header : FITS header
    sections : list
      List of datasec, oscansec, ampsec sections
    """

    # Check for file; allow for extra .gz, etc. suffix
    fil = glob.glob(raw_file+'*') 
    if len(fil) != 1:
        msgs.error("Found {:d} files matching {:s}".format(len(fil)))

    # Read
    msgs.info("Reading LRIS file: {:s}".format(fil[0]))
    hdu = fits.open(fil[0])
    head0 = hdu[0].header

    # Get post, pre-pix values
    precol = head0['PRECOL']
    postpix = head0['POSTPIX']
    preline = head0['PRELINE']
    postline = head0['POSTLINE']

    # Setup for datasec, oscansec
    dsec = []
    osec = []

    # get the x and y binning factors...
    binning = head0['BINNING']
    xbin, ybin = [int(ibin) for ibin in binning.split(',')]

    # First read over the header info to determine the size of the output array...
    n_ext = len(hdu)-1  # Number of extensions (usually 4)
    xcol = []
    xmax = 0
    ymax = 0
    xmin = 10000
    ymin = 10000
    for i in np.arange(1, n_ext+1):
        theader = hdu[i].header
        detsec = theader['DETSEC']
        if detsec != '0':
            # parse the DETSEC keyword to determine the size of the array.
            x1, x2, y1, y2 = np.array(arparse.load_sections(detsec, fmt_iraf=False)).flatten()

            # find the range of detector space occupied by the data
            # [xmin:xmax,ymin:ymax]
            xt = max(x2, x1)
            xmax = max(xt, xmax)
            yt =  max(y2, y1)
            ymax = max(yt, ymax)

            # find the min size of the array
            xt = min(x1, x2)
            xmin = min(xmin, xt)
            yt = min(y1, y2)
            ymin = min(ymin, yt)
            # Save
            xcol.append(xt)

    # determine the output array size...
    nx = xmax - xmin + 1
    ny = ymax - ymin + 1

    # change size for binning...
    nx = nx // xbin
    ny = ny // ybin

    # Update PRECOL and POSTPIX
    precol = precol // xbin
    postpix = postpix // xbin

    # Deal with detectors
    if det in [1,2]:
        nx = nx // 2
        n_ext = n_ext // 2
        det_idx = np.arange(n_ext, dtype=np.int) + (det-1)*n_ext
        ndet = 1
    elif det is None:
        ndet = 2
        det_idx = np.arange(n_ext).astype(int)
    else:
        raise ValueError('Bad value for det')

    # change size for pre/postscan...
    if not TRIM:
        nx += n_ext*(precol+postpix)
        ny += preline + postline

    # allocate output array...
    array = np.zeros( (nx, ny) )
    order = np.argsort(np.array(xcol))

    # insert extensions into master image...
    for kk, i in enumerate(order[det_idx]):

        # grab complete extension...
        data, predata, postdata, x1, y1 = lris_read_amp(hdu, i+1)
                            #, linebias=linebias, nobias=nobias, $
                            #x1=x1, x2=x2, y1=y1, y2=y2, gaindata=gaindata)
        # insert components into output array...
        if not TRIM:
            # insert predata...
            buf = predata.shape
            nxpre = buf[0]
            xs = kk*precol
            xe = xs + nxpre
            '''
            if keyword_set(VERBOSITY) then begin
                section = '['+stringify(xs)+':'+stringify(xe)+',*]'
                message, 'inserting extension '+stringify(i)+ $
                         ' predata  in '+section, /info
            endif 
            '''
            array[xs:xe, :] = predata

            # insert data...
            buf = data.shape
            nxdata = buf[0]
            nydata = buf[1]
            xs = n_ext*precol + kk*nxdata #(x1-xmin)/xbin
            xe = xs + nxdata
            # Data section
            section = '[{:d}:{:d},{:d}:{:d}]'.format(preline,nydata-postline, xs, xe)  # Eliminate lines
            dsec.append(section)
            #print('data',xs,xe)
            array[xs:xe, :] = data   # Include postlines

            #; insert postdata...
            buf = postdata.shape
            nxpost = buf[0]
            xs = nx - n_ext*postpix + kk*postpix
            xe = xs + nxpost 
            section = '[:,{:d}:{:d}]'.format(xs, xe)
            osec.append(section)
            '''
            if keyword_set(VERBOSITY) then begin
                section = '['+stringify(xs)+':'+stringify(xe)+',*]'
                message, 'inserting extension '+stringify(i)+ $
                         ' postdata in '+section, /info
            endif 
            '''
            array[xs:xe, :] = postdata
        else:
            buf = data.shape
            nxdata = buf[0]
            nydata = buf[1]

            xs = (x1-xmin)//xbin
            xe = xs + nxdata 
            ys = (y1-ymin)//ybin
            ye = ys + nydata - postline

            yin1 = preline
            yin2 = nydata - postline 

            '''
            if keyword_set(VERBOSITY) then begin
                section = '['+stringify(xs)+':'+stringify(xe)+ $
                          ','+stringify(ys)+':'+stringify(ye)+']'
                message, 'inserting extension '+stringify(i)+ $
                         ' data     in '+section, /info
            endif 
            '''
            array[xs:xe, ys:ye] = data[:, yin1:yin2]

    # make sure BZERO is a valid integer for IRAF
    obzero = head0['BZERO']
    head0['O_BZERO'] = obzero
    head0['BZERO'] = 32768-obzero

    # Return, transposing array back to goofy Python indexing
    return array.T, head0, (dsec, osec)


def lris_read_amp(inp, ext):
    """
    Read one amplifier of an LRIS multi-extension FITS image

    Parameters
    ----------
    inp: tuple 
      (str,int) filename, extension
      (hdu,int) FITS hdu, extension

    Returns
    -------
    data
    predata
    postdata
    x1
    y1

    ;------------------------------------------------------------------------
    function lris_read_amp, filename, ext, $
      linebias=linebias, nobias=nobias, $
      predata=predata, postdata=postdata, header=header, $
      x1=x1, x2=x2, y1=y1, y2=y2, GAINDATA=gaindata
    ;------------------------------------------------------------------------
    ; Read one amp from LRIS mHDU image
    ;------------------------------------------------------------------------
    """
    # Parse input
    if isinstance(inp, basestring):
        hdu = fits.open(inp)
    else:
        hdu = inp

    # Get the pre and post pix values
    # for LRIS red POSTLINE = 20, POSTPIX = 80, PRELINE = 0, PRECOL = 12
    head0 = hdu[0].header
    precol = head0['precol']
    postpix = head0['postpix']

    # Deal with binning
    binning = head0['BINNING']
    xbin, ybin = [int(ibin) for ibin in binning.split(',')]
    precol = precol//xbin
    postpix = postpix//xbin

    # get entire extension...
    temp = hdu[ext].data.transpose() # Silly Python nrow,ncol formatting
    tsize = temp.shape
    nxt = tsize[0]

    # parse the DETSEC keyword to determine the size of the array.
    header = hdu[ext].header
    detsec = header['DETSEC']
    x1, x2, y1, y2 = np.array(arparse.load_sections(detsec, fmt_iraf=False)).flatten()

    # parse the DATASEC keyword to determine the size of the science region (unbinned)
    datasec = header['DATASEC']
    xdata1, xdata2, ydata1, ydata2 = np.array(arparse.load_sections(datasec, fmt_iraf=False)).flatten()

    # grab the components...
    predata = temp[0:precol, :]
    # datasec appears to have the x value for the keywords that are zero
    # based. This is only true in the image header extensions
    # not true in the main header.  They also appear inconsistent between
    # LRISr and LRISb!
    #data     = temp[xdata1-1:xdata2-1,*]
    #data = temp[xdata1:xdata2+1, :]
    if (xdata1-1) != precol:
        msgs.error("Something wrong in LRIS datasec or precol")
    xshape = 1024 // xbin
    if (xshape+precol+postpix) != temp.shape[0]:
        msgs.error("Wrong size for in LRIS detector somewhere.  Funny binning?")
    data = temp[precol:precol+xshape,:]
    postdata = temp[nxt-postpix:nxt, :]

    # flip in X as needed...
    if x1 > x2:
        xt = x2
        x2 = x1
        x1 = xt
        data = np.flipud(data) #reverse(temporary(data),1)

    # flip in Y as needed...
    if y1 > y2:
        yt = y2
        y2 = y1
        y1 = yt
        data = np.fliplr(data)
        predata = np.fliplr(predata)
        postdata = np.fliplr(postdata)

    '''
    #; correct gain if requested...
    if keyword_set(GAINDATA) then begin
        gain = gainvalue( gaindata, header)
        data = FLOAT(temporary(data)) * gain
        predata = FLOAT(temporary(predata)) * gain
        postdata = FLOAT(temporary(postdata)) * gain
    endif
    '''

    '''
    ;; optional bias subtraction...
    if ~ keyword_set(NOBIAS) then begin
        if keyword_set( LINEBIAS) then begin
            ;; compute a bias for each line...
            bias = median( postdata, dim=1)

            ;; subtract for data...
            buf = size(data)
            nx = buf[1]
            ny = buf[2]
            data2 = fltarr(nx,ny)
            for i=0,nx-1 do begin
                data2[i,*] = float(data[i,*]) - bias
            endfor 
            data = data2
        endif else begin
            ;; compute a scalar bias....
            bias = median( postdata)
            data -= bias
        endelse
    endif
    '''

    return data, predata, postdata, x1, y1


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



def convert_lowredux_pixelflat(infil, outfil):
    """ Convert LowRedux pixelflat to PYPIT format
    Returns
    -------

    """
    # Read
    hdu = fits.open(infil)
    data = hdu[0].data

    #
    prihdu = fits.PrimaryHDU()
    hdus = [prihdu]
    prihdu.header['FRAMETYP'] = 'pixelflat'

    # Detector 1
    img1 = data[:,:data.shape[1]//2]
    hdu = fits.ImageHDU(img1)
    hdu.name = 'DET1'
    prihdu.header['EXT0001'] = 'DET1-pixelflat'
    hdus.append(hdu)

    # Detector 2
    img2 = data[:,data.shape[1]//2:]
    hdu = fits.ImageHDU(img2)
    hdu.name = 'DET2'
    prihdu.header['EXT0002'] = 'DET2-pixelflat'
    hdus.append(hdu)

    # Finish
    hdulist = fits.HDUList(hdus)
    hdulist.writeto(outfil, clobber=True)
    print('Wrote {:s}'.format(outfil))
