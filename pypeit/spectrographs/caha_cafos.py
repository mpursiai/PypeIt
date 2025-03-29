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


class CAHACAFOSspectrograph(spectrograph.Spectrograph):
    """
    Child to handle NOT ALFOSC spectrograph
    """
    ndet = 1
    name = 'caha_cafos'
    telescope = telescopes.NOTTelescopePar()
    camera = 'CAFOS'
    url = 'https://w3.caha.es/CAHA/Instruments/CAFOS/cafos_overview.html'
    header_name = 'CAFOS'
    supported = True
    comment = 'Testing'

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
        # https://www.ing.iac.es/Engineering/detectors/g3_ultra_red%2B2.html
        # NOTE I AM ASSUMING RED+2

        #if hdu is None:
        #    binning = '1,1'
        #    gain = None
        #    ronoise = None
        #else:
        #    binning = self.get_meta_value(self.get_headarr(hdu), 'binning')
        gain = np.atleast_1d(hdu[0].header['CCDSENS'])  # e-/ADU # Checked!
        ronoise = np.atleast_1d(hdu[0].header['CCDRON'])  # e- # Checked!
        bin_x = np.atleast_1d(hdu[0].header['CCDBINX']) # Checked!
        bin_y = np.atleast_1d(hdu[0].header['CCDBINY']) # Checked!
        binning ="%i,%i" % (bin_x, bin_y) # Checked!

        # Detector 1
        detector_dict = dict(
            binning         = binning, # Checked!
            det             = 1, # Checked!
            dataext         = 0, # Checked!
            specaxis        = 0, # Checked!
            specflip        = True, # Checked!
            spatflip        = False, # Checked!
            xgap            = 0., # Checked! I assume
            ygap            = 0., # Checked! I assume
            ysize           = 1., # Checked! I assume
            platescale      = 0.53, # Checked! https://w3.caha.es/CAHA/Instruments/CAFOS/cafos_overview.html
            mincounts       = -1e10, # Typically set to that
            darkcurr        = 8,      # e-/pix/hr # Checked!
            saturation      = 180000.,  # ADU
            nonlinear       = 0.333, #  up to 60000 ADU so one third of saturation? 
            #datasec         = np.atleast_1d('[:,{}:{}]'.format(1, 2048)),  # Unbinned # Checked I think 
            datasec         = np.atleast_1d('[:,1:2048]'),  # Unbinned # Checked I think 
            #datasec         = np.atleast_1d('[:,{}:{}]'.format(550, 2650)),  # Unbinned # Checked I think 
            #oscansec        = np.atleast_1d('[2049:2167,1:2048]'),
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

        par['reduce']['findobj']['find_min_max'] = [690,1500]
        par['reduce']['findobj']['find_trim_edge'] = [250,250]
        par['reduce']['findobj']['snr_thresh'] = 20

        par['reduce']['skysub']['bspline_spacing'] = 0.6 # 0.6 by default, not sure if changing improves it
        par['reduce']['skysub']['sky_sigrej'] = 20.
        par['reduce']['skysub']['no_local_sky'] = False

        par['reduce']['extraction']['sn_gauss'] = 4 # One fiddly ass parameter. 4 is default for NOT, works mostly here except for very bright sources

        # Ignore PCA
        par['calibrations']['slitedges']['sync_predict'] = 'auto'
        par['calibrations']['slitedges']['bound_detector'] = True
        # Flats are sometimes quite ugly due to dust on the slit which leads to the erroneous detection of multiple slits. So set a higher edge_thresh and minimum_slit_gap.
        #par['calibrations']['slitedges']['edge_thresh'] = 100
        #par['calibrations']['slitedges']['minimum_slit_gap'] = 15
        par['calibrations']['slitedges']['exclude_regions'] = "1:0:100,"


        # Set pixel flat combination method
        par['calibrations']['pixelflatframe']['process']['combine'] = 'median'
        # Wavelength calibration methods
        #par['calibrations']['wavelengths']['method'] = 'holy-grail'
        par['calibrations']['wavelengths']['method'] = 'full_template'
        #par['calibrations']['wavelengths']['lamps'] = ['ArI','HeI', 'HgI_DeVeny', 'CdI', "CsI", "RbI"]
        par['calibrations']['wavelengths']['lamps'] = ["CAFOS"]
        par['calibrations']['wavelengths']['sigdetect'] = 10.0
        par['calibrations']['wavelengths']['n_first'] = 5
        par['calibrations']['wavelengths']['n_final'] = 5


        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 1]
        par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['calibrations']['arcframe']['exprng'] = [0, 5]  # Long arc exposures on this telescope
        par['calibrations']['standardframe']['exprng'] = [None, 300]
        par['scienceframe']['exprng'] = [300, None]

        # Multiple arcs with different lamps, so can't median combine nor clip, also need to remove continuum
        par['calibrations']['arcframe']['process']['clip'] = True
        par['calibrations']['arcframe']['process']['combine'] = 'mean'
        par['calibrations']['arcframe']['process']['subtract_continuum'] = True
        par['calibrations']['arcframe']['process']['trim'] = True

        par['calibrations']['tiltframe']['process']['clip'] = False
        par['calibrations']['tiltframe']['process']['combine'] = 'mean'
        par['calibrations']['tiltframe']['process']['subtract_continuum'] = True
       

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
        #self.meta['ra'] = dict(ext=0, card='CAT-RA') # apparently in wrong format
        self.meta['ra'] = dict(card=None, compound=True) 
        self.meta['dec'] = dict(card=None, compound=True) 
        #self.meta['RA'] = dict(ext=0, card='RA')
        #self.meta['dec'] = dict(ext=0, card='DEC')

        #self.meta['dec'] = dict(ext=0, card='CAT-DEC')
        self.meta['target'] = dict(ext=0, card='Object')
        #self.meta['decker'] = dict(ext=0, card='HIERARCH CAHA INS SLIT WID') # Maybe? should be WID
        self.meta['decker'] = dict(card=None, compound=True) #
        self.meta['binning'] = dict(card=None, compound=True)

       # self.meta['mjd'] = dict(ext=0, card=None, compound=True)
        self.meta['mjd'] = dict(ext=0, card='MJD-OBS')
       
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        # Extras for config and frametyping
        self.meta['dispname'] = dict(ext=0, card='INSGRNAM')
        self.meta['idname'] = dict(ext=0, card='IMAGETYP')
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
            binspatial = headarr[0]['CCDBINX']
            binspec = headarr[0]['CCDBINY']
            return parse.binning2string(binspec, binspatial)
        elif meta_key == 'mjd':
            time = headarr[0]['DATE-OBS'] + "T" + headarr[0]['UTSTART'] 
            ttime = Time(time, format='isot')
            return ttime.mjd

        elif meta_key == 'ra':
            objra = headarr[0]['RA'] # Given in hours, not deg
            #c = SkyCoord(objra, "+00:00:00", frame="icrs",unit=(u.hourangle, u.deg))
            return objra

        elif meta_key == 'dec':
            objdec = headarr[0]['dec'] # Given in hours, not deg
            #c = SkyCoord("00:00:00", objdec, frame="icrs",unit=(u.hourangle, u.deg))
            return objdec
       
        elif meta_key == 'decker':
            decker = headarr[0]['HIERARCH CAHA INS SLIT WID'] # Slit in microns
            px_size = headarr[0]['INSTRPIX'] # Pixel size in microns
            scale = headarr[0]['INSTRSCL'] # scale "/px
            slit_width = decker / px_size *scale
            #c = SkyCoord("00:00:00", objdec, frame="icrs",unit=(u.hourangle, u.deg))
            return np.around(slit_width, decimals=0)

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
        return ['GRATNAME', 'DEKKERID', 'CCDXBIN', 'CCDYBIN']

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
        if ftype in ['science', 'standard']:
            return good_exp & (fitstbl['idname'] == 'science')
        if ftype == 'bias':
            return good_exp & (fitstbl['idname'] == 'bias')
        if ftype in ['pixelflat', 'trace', 'illumflat']:
            return good_exp & (fitstbl['idname'] == 'flat')
        if ftype in ['pinhole', 'dark']:
            # Don't type pinhole or dark frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc','tilt']:
            return good_exp & (fitstbl['idname'] == 'arc')
        if ftype in ['sky']:
            return good_exp & (fitstbl['idname'] == 'sky')
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
        if self.get_meta_value(scifile, 'dispname') == 'blue-100' and self.get_meta_value(scifile, 'binning') == "1,1":
            par['calibrations']['wavelengths']['reid_arxiv'] = 'caha_cafos_b100.fits'
        elif self.get_meta_value(scifile, 'dispname') == 'green-100' and self.get_meta_value(scifile, 'binning') == "1,1":
            par['calibrations']['wavelengths']['reid_arxiv'] = 'caha_cafos_g100.fits'
        elif self.get_meta_value(scifile, 'dispname') == 'red-100' and self.get_meta_value(scifile, 'binning') == "1,1":
            par['calibrations']['wavelengths']['reid_arxiv'] = 'caha_cafos_r100.fits'
        elif self.get_meta_value(scifile, 'dispname') == 'blue-200' and self.get_meta_value(scifile, 'binning') == "1,1":
            par['calibrations']['wavelengths']['reid_arxiv'] = 'caha_cafos_b200.fits'        
        elif self.get_meta_value(scifile, 'dispname') == 'green-200' and self.get_meta_value(scifile, 'binning') == "1,1":
            par['calibrations']['wavelengths']['reid_arxiv'] = 'caha_cafos_g200.fits' 
        elif self.get_meta_value(scifile, 'dispname') == 'blue-400' and self.get_meta_value(scifile, 'binning') == "1,1":
            par['calibrations']['wavelengths']['reid_arxiv'] = 'caha_cafos_b400.fits'
        elif self.get_meta_value(scifile, 'dispname') == 'red-400' and self.get_meta_value(scifile, 'binning') == "1,1":
            par['calibrations']['wavelengths']['reid_arxiv'] = 'caha_cafos_r400.fits'
        else:
            msgs.warn('caha_cafos.py: YOU NEED TO ADD IN THE WAVELENGTH SOLUTION FOR THIS GRISM')

        # Return
        return par
