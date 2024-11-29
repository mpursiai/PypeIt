"""
Module for NOT ALFOSC spectrograph

.. include:: ../include/links.rst
"""
from IPython import embed

import numpy as np

from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import units as u
from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.spectrographs import spectrograph
from pypeit.core import parse
from pypeit.images import detector_container


class LTSPRATSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle NOT ALFOSC spectrograph
    """
    ndet = 1
    name = 'lt_sprat'
    telescope = telescopes.LTTelescopePar()
    camera = 'SPRAT'
    url = 'https://telescope.livjm.ac.uk/TelInst/Inst/SPRAT/'
    header_name = 'SPRAT'
    supported = True
    comment = 'No comments'

    def get_detector_par(self, det, hdu=None):
        """
        Return metadata for the selected detector.

        Detector data from `here
        <http://www.not.iac.es/instruments/detectors/CCD14/>`__.

        .. warning::

            Many of the necessary detector parameters are read from the file
            header, meaning the ``hdu`` argument is effectively **required** for
            NOT/ALFOSC.  The optional use of ``hdu`` is only viable for
            automatically generated documentation.

        Args:
            det (:obj:`int`):
                1-indexed detector number.
            hdu (`astropy.io.fits.HDUList`_, optional):
                The open fits file with the raw image of interest.  If not
                provided, frame-dependent parameters are set to a default.

        Returns:
            :class:`~pypeit.images.detector_container.DetectorContainer`:
            Object with the detector metadata.
        """
        # http://www.not.iac.es/instruments/detectors/CCD14/

        if hdu is None:
            binning = '1,1'
            gain = None
            ronoise = None
        else:
            bin_x = np.atleast_1d(hdu[0].header['CCDXBIN']) # Checked!
            bin_y = np.atleast_1d(hdu[0].header['CCDYBIN']) # Checked!
            binning ="%i,%i" % (bin_x, bin_y) # Checked!
            gain =  np.atleast_1d(2.45) #np.atleast_1d(hdu[0].header['GAIN'])  # e-/ADU
            ronoise = np.atleast_1d(5.7)  # e-

        # Detector 1
        detector_dict = dict(
            binning         = binning,
            det             = 1,
            dataext         = 0,
            specaxis        = 1,
            specflip        = False,
            spatflip        = False,
            xgap            = 0.,
            ygap            = 0.,
            ysize           = 1.,
            platescale      = 0.44,
            mincounts       = -1e10,
            darkcurr        = 0,      # e-/pix/hr
            saturation      = 63038.,  # ADU
            nonlinear       = 1.0,
            datasec         = np.atleast_1d('[:,{}:{}]'.format(1, 1024)),  # Unbinned
            oscansec        = None,
            numamplifiers   = 1,
            gain            = gain,     # e-/ADU
            ronoise         = ronoise   # e-
        )

#        # Parse datasec, oscancsec from the header
#        head1 = hdu[1].header
#        detector_dict['gain'] = np.atleast_1d(head1['GAIN'])  # e-/ADU
#        detector_dict['ronoise'] = np.atleast_1d(head1['RDNOISE'])  # e-

        # Return
        return detector_container.DetectorContainer(**detector_dict)

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.

        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()


        par['reduce']['findobj']['find_min_max'] = [400,500]
        par['reduce']['findobj']['find_trim_edge'] = [50,50]
        par['reduce']['findobj']['snr_thresh'] = 4
        par['reduce']['skysub']['bspline_spacing'] = 0.6
        par['reduce']['skysub']['sky_sigrej'] = 20.
        par['reduce']['skysub']['no_local_sky'] = True

        #par['reduce']['skysub']['local_maskwidth'] = 40.
        #par['reduce']['skysub']['no_local_sky'] = True
        par['reduce']['extraction']['sn_gauss'] = 10 # One fiddly ass parameter. 4 is default for NOT, works mostly here except for very bright sources
        # Ignore PCA
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'
        par['calibrations']['slitedges']['bound_detector'] = True
        # Flats are sometimes quite ugly due to dust on the slit which leads to the erroneous detection of multiple slits. So set a higher edge_thresh and minimum_slit_gap.
        par['calibrations']['slitedges']['edge_thresh'] = 30
        par['calibrations']['slitedges']['minimum_slit_gap'] = 15

        # Set pixel flat combination method
        par['calibrations']['pixelflatframe']['process']['combine'] = 'median'
        # Wavelength calibration methods
        #par['calibrations']['wavelengths']['method'] = 'holy-grail'
        par['calibrations']['wavelengths']['method'] = 'full_template'
        par['calibrations']['wavelengths']['lamps'] = ['XeI']
        par['calibrations']['wavelengths']['sigdetect'] = 10.0
        # Set the default exposure time ranges for the frame typing
        #par['calibrations']['darkframe']['useframe'] = "False"
        #par['calibrations']['pinholeframe']['useframe'] = "False"

        par['calibrations']['biasframe']['exprng'] = [None, 1]
        par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['calibrations']['arcframe']['exprng'] = [None, None]  # Long arc exposures on this telescope
        par['calibrations']['standardframe']['exprng'] = [None, 120]
        par['scienceframe']['exprng'] = [10, None]

        # Multiple arcs with different lamps, so can't median combine nor clip, also need to remove continuum
        par['calibrations']['arcframe']['process']['clip'] = False
        par['calibrations']['arcframe']['process']['combine'] = 'mean'
        par['calibrations']['arcframe']['process']['subtract_continuum'] = True
        par['calibrations']['tiltframe']['process']['clip'] = False
        par['calibrations']['tiltframe']['process']['combine'] = 'mean'
        par['calibrations']['tiltframe']['process']['subtract_continuum'] = True


        par['scienceframe']['process']['use_biasimage'] = False
        par['scienceframe']['process']['use_darkimage'] = False
        par['scienceframe']['process']['use_illumflat'] = False
        par['scienceframe']['process']['use_overscan'] = False
        par['scienceframe']['process']['use_pattern'] = False
        par['scienceframe']['process']['use_pixelflat'] = False
        par['scienceframe']['process']['use_specillum'] = False

        par['calibrations']['standardframe']['process']['use_biasimage'] = False
        par['calibrations']['standardframe']['process']['use_darkimage'] = False
        par['calibrations']['standardframe']['process']['use_illumflat'] = False
        par['calibrations']['standardframe']['process']['use_overscan'] = False
        par['calibrations']['standardframe']['process']['use_pattern'] = False
        par['calibrations']['standardframe']['process']['use_pixelflat'] = False
        par['calibrations']['standardframe']['process']['use_specillum'] = False


        #par['baseprocess']['use_biasimage'] = False
        #par['calibrations']['biasframe'] = 
        #par['calibrations']['biasframe'] = 
        #par['calibrations']['biasframe'] = 
        #par['calibrations']['biasframe'] = 
        #par['calibrations']['biasframe'] = 
        #par['calibrations']['biasframe'] = 


        # No overscan region!
        turn_off = dict(use_overscan=False)
        par.reset_all_processimages_par(**turn_off)

        return par

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}
        # Required (core)
        self.meta['ra'] = dict(ext=0, card=None, compound=True)
        self.meta['dec'] = dict(ext=0, card=None, compound=True)
        self.meta['target'] = dict(ext=0, card='CAT-NAME')
        self.meta['decker'] = dict(ext=0, card='FOVTYPE')
        self.meta['binning'] = dict(card=None, compound=True)

        self.meta['mjd'] = dict(ext=0, card='MJD')
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        # Extras for config and frametyping
        self.meta['dispname'] = dict(ext=0, card='GRATSET')
        self.meta['idname'] = dict(ext=0, card='OBSTYPE')
        self.meta['instrument'] = dict(ext=0, card='INSTRUME')

    def compound_meta(self, headarr, meta_key):
        """
        Methods to generate metadata requiring interpretation of the header
        data, instead of simply reading the value of a header card.

        Args:
            headarr (:obj:`list`):
                List of `astropy.io.fits.Header`_ objects.
            meta_key (:obj:`str`):
                Metadata keyword to construct.

        Returns:
            object: Metadata value read from the header(s).
        """
        if meta_key == 'binning':
            # PypeIt frame
            binspatial = headarr[0]['CCDXBIN']
            binspec = headarr[0]['CCDYBIN']
            return parse.binning2string(binspec, binspatial)
    
        elif meta_key == 'ra':
            objra = headarr[0]['CAT-RA'] # Given in hours, not deg
            c = SkyCoord(objra, "+00:00:00", frame="icrs",unit=(u.hourangle, u.deg))
            return c.ra.degree

        elif meta_key == 'dec':
            objdec = headarr[0]['CAT-DEC'] # Given in hours, not deg
            c = SkyCoord("00:00:00", objdec, frame="icrs",unit=(u.hourangle, u.deg))
            return c.dec.degree
        msgs.error("Not ready for this compound meta")

    def configuration_keys(self):
        """
        Return the metadata keys that define a unique instrument
        configuration.

        This list is used by :class:`~pypeit.metadata.PypeItMetaData` to
        identify the unique configurations among the list of frames read
        for a given reduction.

        Returns:
            :obj:`list`: List of keywords of data pulled from file headers
            and used to constuct the :class:`~pypeit.metadata.PypeItMetaData`
            object.
        """
        return ['dispname', 'decker', 'binning']

    def raw_header_cards(self):
        """
        Return additional raw header cards to be propagated in
        downstream output files for configuration identification.

        The list of raw data FITS keywords should be those used to populate
        the :meth:`~pypeit.spectrographs.spectrograph.Spectrograph.configuration_keys`
        or are used in :meth:`~pypeit.spectrographs.spectrograph.Spectrograph.config_specific_par`
        for a particular spectrograph, if different from the name of the
        PypeIt metadata keyword.

        This list is used by :meth:`~pypeit.spectrographs.spectrograph.Spectrograph.subheader_for_spec`
        to include additional FITS keywords in downstream output files.

        Returns:
            :obj:`list`: List of keywords from the raw data files that should
            be propagated in output files.
        """
        return ['CONFNAME', 'FOVTYPE', 'CCDXBIN', 'CCDYBIN']

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.

        Args:
            ftype (:obj:`str`):
                Type of frame to check. Must be a valid frame type; see
                frame-type :ref:`frame_type_defs`.
            fitstbl (`astropy.table.Table`_):
                The table with the metadata for one or more frames to check.
            exprng (:obj:`list`, optional):
                Range in the allowed exposure time for a frame of type
                ``ftype``. See
                :func:`pypeit.core.framematch.check_frame_exptime`.

        Returns:
            `numpy.ndarray`_: Boolean array with the flags selecting the
            exposures in ``fitstbl`` that are ``ftype`` type frames.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        if ftype == 'science':
            return good_exp & (fitstbl['idname'] == 'EXPOSE')
        
        if ftype == 'standard':
            return good_exp & ((fitstbl['idname'] == 'STANDARD'))
        
        if ftype == 'bias':
            return good_exp & (fitstbl['idname'] == 'BIAS')
        
        if ftype in ['pixelflat', 'trace', 'illumflat']:
            return good_exp & (fitstbl['idname'] == 'LAMP-FLAT')
        if ftype in ['pinhole', 'dark']:
            # Don't type pinhole or dark frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc','tilt']:
            return good_exp & (fitstbl['idname'] == 'ARC')
        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

    def config_specific_par(self, scifile, inp_par=None):
        """
        Modify the PypeIt parameters to hard-wired values used for
        specific instrument configurations.

        Args:
            scifile (:obj:`str`):
                File to use when determining the configuration and how
                to adjust the input parameters.
            inp_par (:class:`~pypeit.par.parset.ParSet`, optional):
                Parameter set used for the full run of PypeIt.  If None,
                use :func:`default_pypeit_par`.

        Returns:
            :class:`~pypeit.par.parset.ParSet`: The PypeIt parameter set
            adjusted for configuration specific parameter values.
        """
        # Start with instrument wide
        par = super().config_specific_par(scifile, inp_par=inp_par)

        # Wavelength calibrations
        if self.get_meta_value(scifile, 'dispname') == 'blue':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'lt_sprat_blue_6000.fits'
 
        else:
            msgs.warn('lt_sprat.py: YOU NEED TO ADD IN THE WAVELENGTH SOLUTION FOR THIS GRISM')

        # Return
        return par
