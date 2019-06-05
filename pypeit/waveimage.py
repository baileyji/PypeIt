"""
Module for guiding construction of the Wavelength Image

.. _numpy.ndarray: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html

"""
import os
import inspect

import numpy as np

from pypeit import msgs
from pypeit import utils
from pypeit.masterframe import MasterFrame
from pypeit import ginga
from pypeit.core import pixels
from pypeit.core import trace_slits

from IPython import embed

class WaveImage(MasterFrame):
    """
    Class to generate the Wavelength Image

    Args:
        tslits_dict (dict): dict from TraceSlits class (e.g. slitpix)
        tilts (np.ndarray): Tilt image
        wv_calib (dict): wavelength solution dictionary
            Parameters are read from wv_calib['par']
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            The `Spectrograph` instance that sets the
            instrument used to take the observations.  Used to set
            :attr:`spectrograph`.
        master_key (:obj:`str`, optional):
            The string identifier for the instrument configuration.  See
            :class:`pypeit.masterframe.MasterFrame`.
        master_dir (str, optional): Path to master frames
        maskslits (ndarray, optional): True = skip this slit
        reuse_masters (bool, optional):  Load from disk if possible

    Attributes:
        frametype : str
          Hard-coded to 'wave'
        mswave (ndarray): Wavelength image
        steps (list): List of the processing steps performed

    """
    # Frametype is a class attribute
#    frametype = 'wave'
    master_type = 'Wave'

    def __init__(self, tslits_dict, tilts, wv_calib, spectrograph, maskslits,
                 master_key=None, master_dir=None, reuse_masters=False):

        # MasterFrame
        MasterFrame.__init__(self, self.master_type, master_dir=master_dir,
                                         master_key=master_key, reuse_masters=reuse_masters)

        # Required parameters
        self.spectrograph = spectrograph

        self.tslits_dict = tslits_dict
        self.tilts = tilts
        self.wv_calib = wv_calib
        self.slitmask = pixels.tslits2mask(self.tslits_dict) if tslits_dict is not None else None
        # TODO: only echelle is ever used.  Do we need to keep the whole
        # thing?
        self.par = wv_calib['par'] if wv_calib is not None else None

        self.maskslits = maskslits

        # For echelle order, primarily
        self.slit_spat_pos = trace_slits.slit_spat_pos(self.tslits_dict)

        # List to hold ouptut from inspect about what module create the image?
        self.steps = []


        # Main output
        self.mswave = None

    def build_wave(self):
        """
        Main algorithm to build the wavelength image

        Returns:
            `numpy.ndarray`_: The wavelength image.

        """
        # Loop on slits
        ok_slits = np.where(~self.maskslits)[0]
        self.mswave = np.zeros_like(self.tilts)
        nspec = self.slitmask.shape[0]

        # Error checking on the wv_calib
        #if (nspec-1) != int(self.wv_calib[str(0)]['fmax']):
        #    msgs.error('Your wavelength fits used inconsistent normalization. Something is wrong!')

        # If this is echelle print out a status message and do some error checking
        if self.par['echelle']:
            msgs.info('Evaluating 2-d wavelength solution for echelle....')
            if len(self.wv_calib['fit2d']['orders']) != len(ok_slits):
                msgs.error('wv_calib and ok_slits do not line up. Something is very wrong!')

        # Unpack some 2-d fit parameters if this is echelle
        for slit in ok_slits:
            thismask = (self.slitmask == slit)
            if self.par['echelle']:
                order = self.spectrograph.slit2order(self.slit_spat_pos[slit])
                # evaluate solution
                self.mswave[thismask] = utils.func_val(self.wv_calib['fit2d']['coeffs'],
                                                       self.tilts[thismask],
                                                       self.wv_calib['fit2d']['func2d'],
                                                       x2=np.ones_like(self.tilts[thismask])*order,
                                                       minx=self.wv_calib['fit2d']['min_spec'],
                                                       maxx=self.wv_calib['fit2d']['max_spec'],
                                                       minx2=self.wv_calib['fit2d']['min_order'],
                                                       maxx2=self.wv_calib['fit2d']['max_order'])
                self.mswave[thismask] /= order
            else:
                iwv_calib = self.wv_calib[str(slit)]
                self.mswave[thismask] = utils.func_val(iwv_calib['fitc'], self.tilts[thismask],
                                                       iwv_calib['function'],
                                                       minx=iwv_calib['fmin'],
                                                       maxx=iwv_calib['fmax'])

        # Return
        self.steps.append(inspect.stack()[0][3])
        return self.mswave

    def show(self, item='wave'):
        """
        Show the image

        Args:
            item (str, optional):

        Returns:

        """
        if item == 'wave':
            if self.mswave is not None:
                ginga.show_image(self.mswave)
        else:
            msgs.warn("Not able to show this type of image")

    def __repr__(self):
        # Generate sets string
        txt = '<{:s}: >'.format(self.__class__.__name__)
        return txt

    def save(self, outfile=None, overwrite=True, mswave=None):
        """
        Save the master wavelength image.

        Args:
            outfile (:obj:`str`, optional):
                Name for the output file.  Defaults to
                :attr:`file_path`.
            overwrite (:obj:`bool`, optional):
                Overwrite any existing file.
        """
        _mswave = self.mswave if mswave is None else mswave
        super(WaveImage, self).save(_mswave, 'WAVE', outfile=outfile, overwrite=overwrite,
                                    steps=self.steps)

    # TODO: it would be better to have this instantiate the full class
    # as a classmethod.
    def load(self, ifile=None, return_header=False):
        """
        Load the wavelength image data from a saved master frame.

        Args:
            ifile (:obj:`str`, optional):
                Name of the master frame file.  Defaults to
                :attr:`file_path`.
            return_header (:obj:`bool`, optional):
                Return the header

        Returns:
            tuple: Returns an `numpy.ndarray`_ with the wavelength
            image.  Also returns the primary header, if requested.
        """
        return super(WaveImage, self).load('WAVE', ifile=ifile, return_header=return_header)

    @staticmethod
    def load_from_file(filename, return_header=False):
        """
        Load the wavelength image data from a saved master frame.

        Args:
            filename (:obj:`str`, optional):
                Name of the master frame file.
            return_header (:obj:`bool`, optional):
                Return the header

        Returns:
            tuple: Returns an `numpy.ndarray`_ with the wavelength
            image.  Also returns the primary header, if requested.
        """
        # Use of super() to call staticmethods of the base class seems
        # like a bit of a mess (bound vs. unbound methods).  There's a
        # syntax that works, but for now, I'm just going to call the
        # static method explicitly without using super().  See:
        # https://stackoverflow.com/questions/26788214/super-and-staticmethod-interaction
        return MasterFrame.load_from_file(filename, 'WAVE', return_header=return_header)

