#re-writing PLs writehlsp function to do what I want it to

def writehlsp(star_or_spectbl, components=True, overwrite=False):
    """
    Writes spectbl to a standardized MUSCLES FITS file format that also
    includes all the keywords required for the archive.

    Parameters
    ----------
    spectbl : astropy table in MUSCLES format
        The spectrum to be written to a MSUCLES FITS file.
    name : str
        filename for the FITS output
    overwrite : {True|False}
        whether to overwrite if output file already exists

    Returns
    -------
    None
    """

    if type(star_or_spectbl) is str:
        star = star_or_spectbl
        pfs = db.allpans(star)
        pan = read(filter(lambda s: 'native_resolution' in s, pfs))[0]
        writehlsp(pan, components=components, overwrite=overwrite)
        dpan = read(filter(lambda s: 'dR=' in s, pfs))[0]
        writehlsp(dpan, components=False)
        return
    else:
        spectbl = star_or_spectbl

    star = spectbl.meta['STAR']
    srcspecs = spectbl.meta['SOURCESPECS']
    name = spectbl.meta['NAME']
    pan = 'panspec' in name
    mod = 'mod' in name

    hlspname = db.hlsppath(name)

    # add a wavelength midpoint column
    w = (spectbl['w0'] + spectbl['w1']) / 2.0
    spectbl['w'] = w
    spectbl['w'].description = 'midpoint of the wavelength bin'
    spectbl['w'].unit = 'Angstrom'

    # CREATE PRIMARY EXTENSION
    prihdr = fits.Header()
    if pan:
        prihdr['TELESCOP'] = 'MULTI'
        prihdr['INSTRUME'] = 'MULTI'
        prihdr['GRATING'] = 'MULTI'

        insts = []
        for specname in srcspecs:
            tel = db.parse_observatory(specname)
            spec = db.parse_spectrograph(specname)
            grating = db.parse_info(specname, 3, 4)
            insts.append([rc.HLSPtelescopes[tel], rc.HLSPinstruments[spec], rc.HLSPgratings[grating]])
        insts = set(tuple(inst) for inst in insts)

        for i,inst in enumerate(insts):
            tkey, ikey, gkey = 'TELESC{:02d}'.format(i), 'INSTRU{:02d}'.format(i), 'GRATIN{:02d}'.format(i)
            prihdr[tkey], prihdr[ikey], prihdr[gkey] = inst
    else:
        prihdr['TELESCOP'] = rc.HLSPtelescopes[db.parse_observatory(name)]
        prihdr['INSTRUME'] = rc.HLSPinstruments[db.parse_spectrograph(name)]
        prihdr['GRATING'] = rc.HLSPgratings[db.parse_grating(name)]


        if 'hst' in name:
            # clooge. meh.
            band = name[0]
            f = db.findfiles(band, name, fullpaths=True)[0]
            aper_key = 'APER_ID' if 'fos' in name else 'APERTURE'
            if 'custom' in name or 'coadd' in name:
                srcspecs = fits.getdata(f, 'sourcespecs')
                srcids = [db.parse_id(s) for s in srcspecs['sourcespecs']]
                srcpaths = [db.sourcespecfiles(star, id)[0] for id in srcids]
                apertures = [fits.getval(sf, aper_key) for sf in srcpaths]
                assert len(set(apertures)) == 1
                prihdr['APERTURE'] = apertures[0]
            else:
                prihdr['APERTURE'] = fits.getval(f, aper_key)
        if 'xmm' in name or 'cxo' in name:
            hdr = fits.getheader(db.name2path(name))
            prihdr['GRATING'] = 'NA'
            if 'multi' in name:
                prihdr['DETECTOR'] = 'MULTI'
                prihdr['DETECT00'] = 'PN'
                prihdr['DETECT01'] = 'MOS1'
                prihdr['DETECT02'] = 'MOS2'
                prihdr['FILTER'] = 'MULTI'
                prihdr['FILTER00'] = hdr['pn_filter']
                prihdr['FILTER01'] = hdr['mos1_filter']
                prihdr['FILTER02'] = hdr['mos2_filter']
            if 'pn' in name:
                prihdr['DETECTOR'] = 'PN'
                prihdr['FILTER'] = hdr['pn_filter']
            if 'acs' in name:
                prihdr['DETECTOR'] = hdr['DETNAM']
                prihdr['FILTER'] = 'OBF'

    prihdr['TARGNAME'] = star.upper()
    prihdr['RA_TARG'] = rc.starprops['RA'][star]
    prihdr['DEC_TARG'] = rc.starprops['dec'][star]
    prihdr['PROPOSID'] = 13650
    prihdr['HLSPNAME'] = 'Measurements of the Ultraviolet Spectral Characteristics of Low-mass Exoplanet Host Stars'
    prihdr['HLSPACRN'] = 'MUSCLES'
    prihdr['HLSPLEAD'] = 'David J. Wilson'
    prihdr['PR_INV_L'] = 'Froning'
    prihdr['PR_INV_F'] = 'Cynthia'

    if not (pan or mod):
        mjd0 = np.min(spectbl['minobsdate'])
        mjd1 = np.max(spectbl['maxobsdate'])
        if np.isfinite(mjd0):
            date0 = Time(mjd0, format='mjd')
            prihdr['DATE-OBS'] = date0.isot
            prihdr['EXPSTART'] = date0.mjd
        if np.isfinite(mjd1):
            date1 = Time(mjd1, format='mjd')
            prihdr['EXPEND'] =  date1.mjd
        expt = spectbl['exptime']
        if 'xmm' in name:
            prihdr['EXPTIME'] = expt[0]
            prihdr['EXPDEFN'] = 'MEAN'
        if 'cxo' in name:
            prihdr['EXPTIME'] = expt[0]
        if not np.allclose(expt, expt[0]):
            expmed = np.median(expt)
            prihdr['EXPTIME'] = expmed
            prihdr['EXPDEFN'] = 'MEDIAN'
            prihdr['EXPMAX'] = np.max(expt)
            prihdr['EXPMIN'] = np.min(expt[expt > 0])
            prihdr['EXPMED'] = expmed
        else:
            prihdr['EXPTIME'] = expt[0]

    if not (pan or mod) or 'phx' in name:
        try:
            inst = db.parse_instrument(name)
            normfac = rc.normfacs[star][inst][0]
            panspec = readpan(star)
            insti = rc.getinsti(inst)
            assert insti == spectbl['instrument'][0]
            normfac_vec = panspec['normfac'][panspec['normfac'] == insti]
            if len(normfac_vec) > 0:
                assert np.isclose(normfac_vec, normfac)
        except KeyError:
            normfac = 1.0
        prihdr['normfac'] = (normfac, 'normalization factor used by MUSCLES')


    prihdr['WAVEMIN'] = w[0]
    prihdr['WAVEMAX'] = w[-1]
    prihdr['WAVEUNIT'] = 'ang'
    prihdr['AIRORVAC'] = 'vac'

    if not pan or 'constant' in name:
        mid = len(w) / 2
        prihdr['SPECRES'] = w[mid]
        prihdr['WAVERES'] = w[mid+1] - w[mid]

    prihdr['FLUXMIN'] = np.min(spectbl['flux'])
    prihdr['FLUXMAX'] = np.max(spectbl['flux'])
    prihdr['FLUXUNIT'] = 'erg/s/cm2/ang' if 'phx' not in name else 'arbitrary'

    # CREATE SPECTRUM EXTENSION
    spechdr = fits.Header()
    spechdr['EXTNAME'] = 'SPECTRUM'
    spechdr['EXTNO'] = 2

    cols = ['w', 'w0', 'w1', 'flux']
    descriptions = ['midpoint of the wavelength bin',
                    'left/blue edge of the wavelength bin',
                    'right/red edge of the wavelength bin',
                    'average flux over the bin']
    fitsnames = ['WAVELENGTH', 'WAVELENGTH0', 'WAVELENGTH1', 'FLUX']
    fmts = ['D']*4

    if 'mod' not in name:
        cols.extend(['error', 'exptime', 'flags', 'minobsdate', 'maxobsdate'])
        descriptions.extend(['error on the flux',
                             'cumulative exposure time for the bin',
                             'data quality flags (HST data only)',
                             'modified julian date of start of first exposure',
                             'modified julian date of end of last exposure'])
        fitsnames.extend(['ERROR', 'EXPTIME', 'DQ', 'EXPSTART', 'EXPEND'])
        fmts.extend(['D']*2 + ['I'] + ['D']*2)

    if pan:
        # add a normalized flux column
        spectbl = utils.add_normflux(spectbl)
        spectbl['normflux'].unit = 'Angstrom-1'
        spectbl['normerr'].unit = 'Angstrom-1'
        prihdr['BOLOFLUX'] = utils.bolo_integral(spectbl.meta['STAR'])

        # add header keywords for lorentzian fit
        prihdr['LNZ_NORM'] = spectbl.meta['LNZ_NORM']
        prihdr['LNZ_GAM'] = spectbl.meta['LNZ_GAM']

        cols.extend(['instrument', 'normfac', 'normflux', 'normerr'])
        descriptions.extend(['bitmask identifying the source instrument(s). See "instlgnd" extension for a legend.',
                             'normalization factor applied to the source spectrum',
                             'flux density normalized by the bolometric flux',
                             'error on bolometrically-normalized flux density'])
        fitsnames.extend(['INSTRUMENT', 'NORMFAC', 'BOLOFLUX', 'BOLOERR'])
        fmts.extend(['J', 'D', 'D', 'D'])

    for i, desc in enumerate(descriptions):
        spechdr['TDESC' + str(i+1)] = desc

    if 'COMMENT' in spectbl.meta and len(spectbl.meta['COMMENT']) > 1 and not pan:
        spechdr['COMMENT'] = spectbl.meta['COMMENT']

    datas = [spectbl[col].data for col in cols]
    units = [spectbl[col].unit.to_string() for col in cols]
    fitscols = [fits.Column(array=a, name=n, format=fmt, unit=u) for a, n, fmt, u in zip(datas, fitsnames, fmts, units)]
    spechdu = fits.BinTableHDU.from_columns(fitscols, header=spechdr)
    spechdu.name = 'SPECTRUM'
    for fname, data in zip(fitsnames, datas):
        spechdu.data[fname] = data

    prihdu = fits.PrimaryHDU(header=prihdr)
    hdus = [prihdu, spechdu]

    # INSTRUMENT LEGEND
    if pan:
        lgndhdr = fits.Header()
        lgndhdr['comment'] = legendcomment
        lgndhdr['extno'] = 3
        lgndhdr['comment'] = ('Not all of these instruments were used to acquire data for this particular spectrum. '
                              'Therefore, not all the listed HLSP files will exist in the database. Also note that '
                              'polynomial fits for filling spectral gaps were not saved as separate spectra.')

        vals = rc.instvals
        instnames = rc.instruments
        pieces = [s.split('_') for s in instnames]
        tels, insts, gratings = zip(*pieces)
        tels = [rc.HLSPtelescopes[t] for t in tels]
        insts = [rc.HLSPinstruments[inst] for inst in insts]
        gratings = [rc.HLSPgratings[g] for g in gratings]
        dummynames = ['-_' + s + '_' + star for s in instnames]
        hlspnames = [path.basename(db.hlsppath(n)) for n in dummynames]

        names = ['BITVALUE', 'TELESCOPE', 'INSTRUMENT', 'GRATING', 'HLSP_FILE']
        datas = [vals, tels, insts, gratings, hlspnames]
        lens = [max(map(len, d)) for d in datas[1:]]
        fmts = ['J'] + [str(n) + 'A' for n in lens]
        fitscols = [fits.Column(n, fmt, array=a) for n, fmt, a in zip(names, fmts, datas)]

        lgndhdu = fits.BinTableHDU.from_columns(fitscols, header=lgndhdr)
        lgndhdu.name = 'INSTLGND'

        hdus.append(lgndhdu)

        if components:
            specs, lyaspec = read_panspec_sources(star)
            if lyaspec is not None: specs.append(lyaspec)
            for inst in instnames:
                spec = filter(lambda s: inst in s.meta['NAME'], specs)
                if len(spec) == 0:
                    continue
                assert len(spec) == 1
                spec = spec[0]
                writehlsp(spec, overwrite=overwrite)

    # SOURCE SPECTRA LIST
    if 'hst' in name:
        srchdr = fits.Header()
        srchdr['COMMENT'] = ('This extension contains a list of HST rootnames (9 character string in HST files '
                             'downloaded from MAST) and dataset IDs of the exposures used to create this spectrum '
                             'file. The dataset IDs can be used to directly locate the observations through the MAST '
                             'HST data archive search interface. Multiple identifiers indicate the spectra were '
                             'coadded.')
        srchdr['EXTNO'] = 3
        specnames = spectbl.meta['SOURCESPECS']
        if len(specnames) == 0: specnames = [name]
        rootnames = [s.split('_')[5] for s in specnames]
        files = [db.choosesourcespecs(db.findfiles(band, star, rn))[0] for rn in rootnames]
        id_key = 'ROOTNAME' if 'fos' in name else 'ASN_ID'
        dataids = [fits.getval(f, id_key) for f in files]
        custom = [('custom' in s) or ('x2d' in s) for s in specnames]
        assert all(custom) or (not any(custom))
        srchdr['CUSTOM'] = custom[0], 'spectrum extracted from x2d (bad x1d)'

        if 'gj551' in name:
            if 'g230lb' in name:
                rootnames = dataids = ['OCR7QQANQ', 'OCR7QQANQ']
            if 'g430l' in name:
                rootnames = dataids = ['OCR7QQAOQ', 'OCR7QQAPQ']
            if 'g750l' in name:
                rootnames = dataids = ['OCR7QQARQ', 'OCR7QQASQ', 'OCR7QQAQQ']

        fitscols = [fits.Column(name='ROOTNAME', format='9A', array=rootnames),
                    fits.Column(name='DATASET_ID', format='9A', array=dataids)]
        srchdu = fits.BinTableHDU.from_columns(fitscols, header=srchdr)
        srchdu.name = 'SRCSPECS'

        hdus.append(srchdu)
    if 'xmm' in name or 'cxo' in name:
        srchdr = fits.Header()
        if 'xmm' in name:
            srchdr['COMMENT'] = ('This extension contains a list of observation IDs (DATASET_ID used for consistency '
                                 'with HST data) that can be used to locate the data in the XMM archives. XMM data '
                                 'all come from only a single observation (unlike the HST observations), '
                                 'but this extension is retained in place of a keyword for consistency with the HST '
                                 'files.')
        if 'cxo' in name:
            srchdr['COMMENT'] = ('This extension contains a list of observation IDs (DATASET_ID used for consistency '
                                 'with HST data) that can be used to locate the data in the CXO archives.')
        srchdr['EXTNO'] = 3
        obsids = _parse_keys_sequential(hdr, 'OBS_ID')
        col = fits.Column(name='DATASET_ID', format='10A', array=obsids)
        srchdu = fits.BinTableHDU.from_columns([col], header=srchdr)
        srchdu.name = 'SRCSPECS'
        hdus.append(srchdu)

    hdus = fits.HDUList(hdus)
    hdus.writeto(hlspname, clobber=overwrite)
