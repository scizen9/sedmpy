"""Generate Makefile for reducing ifu data"""

import os
import sys

import astropy.io.fits as pf


def extract_info(infiles):

    headers = []

    print("-- Plan.py: Ingesting headers --")

    for ix, ifile in enumerate(infiles):
        FF = pf.open(ifile)
        FF[0].header['filename'] = ifile
        if 'JD' not in FF[0].header:
            # print("Skipping %s" % ifile)
            continue
        headers.append(FF[0].header)
        FF.close()

    return sorted(headers, key=lambda x: x['JD'])


def identify_observations(headers):
    """Return a list of object name, observation number, and list of files.

    e.g. will return:

    {'STD-BD+25d4655': {1: ['...']}, {2: ['...']}, 
           'ZTF14dvo': {1: ['...', '...']}}
    
    where STD-BD+25d4655 was observed at the beginning and end of night. SN
    14dov was observed once with A-B.
    """
    JD = 0.

    objcnt = {}
    objs = {}
    calibs = {}

    for header in headers:
        if header['JD'] < JD:
            raise Exception("Headers not sorted by JD")
        JD = header['JD']

        fname = header['filename']
        obj = header['OBJECT'].lstrip()
        name = header['NAME'].lstrip()
        exptime = header['exptime']
        adcspeed = header['ADCSPEED']
        if "test" in obj or "Test" in obj or "TEST" in obj:
            continue
        if "Calib" in obj or "bias" in obj:

            def appendToCalibs(Str):

                if Str in obj:
                    if "bias" in Str and exptime == 0.:
                        Str = "%s%1.1f" % (Str, adcspeed)
                        prefix = ""
                        suffix = ""
                    elif "Xe" in Str or "Hg" in Str or "Cd" in Str or \
                                    "Ne" in Str or "dome" in Str:
                        prefix = "b_"
                        suffix = ""
                    else:
                        prefix = "crr_b_"
                        suffix = ""

                    if "bias" in Str and exptime != 0.:
                        print("Mis-labeled bias with exptime > 0: %9.1f" %
                              exptime)
                    else:
                        calibs[Str] = calibs.get(Str, [])
                        calibs[Str].append(prefix + fname + suffix)

            appendToCalibs("bias")
            appendToCalibs("dome")
            appendToCalibs("Xe")
            appendToCalibs("Hg")
            appendToCalibs("Cd")
            appendToCalibs("Ne")
            appendToCalibs("twilight")

        if "Focus:" in obj:
            continue
        if "dark" in obj:
            continue
        if "Calib" in obj:
            continue
        if "STOW" in name:
            continue
        if obj.rstrip() == "":
            continue
        name = name.replace(" ", "_")
        name = name.replace(")", "_")
        name = name.replace("(", "_")
        name = name.replace("[", "_")
        name = name.replace("]", "_")
        name = name.replace("/", "_")
        name = name.replace(":", "_")
        name = name.replace('"', "")

        # The 'A' position defines the start of an object set
        if '[A]' in obj or name not in objcnt:
            cnt = objcnt.get(name, 0) + 1
            vals = objs.get(name, {})
            objcnt[name] = cnt
            vals[cnt] = [fname]
            objs[name] = vals
        else:
            cnt = objcnt[name]
            objs[name][cnt].append(fname)

    print("\n-- Calibrations --")
    for k, v in calibs.items():
        print("%15s : %2.0i" % (k, len(v)))

    print("\n-- Standard Star Sets --")
    for k, v in objs.items():
        if "STD-" in k:
            print("%20s : %2.0i" % (k, len(v[1])))

    print("\n-- Science Object Sets --")
    for k, v in objs.items():
        if "STD-" not in k:
            print("%20s : %2.0i" % (k, len(v)))

    return objs, calibs


make_preamble = """
PY = ~/spy
PYC = ~/sedmpy/drpifu
IMCOMBINE = $(PY) $(PYC)/Imcombine.py
REPORT = $(PY) $(PYC)/DrpReport.py
CLASS = $(PY) $(PYC)/Classify.py
ZTFUPLOAD = $(PY) $(PYC)/growth.py

BSUB = $(PY) $(PYC)/Debias.py
CRRSUB =  $(PY) $(PYC)/CosmicX.py

SRCS = $(wildcard ifu*fits)
BIAS = $(addprefix b_,$(SRCS))
CRRS = $(addprefix crr_,$(BIAS))

mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))
current_dir := $(notdir $(patsubst %/,%,$(dir $(mkfile_path))))

crr_b_% : b_%
	$(CRRSUB) --niter 4 --sepmed --gain 1.0 --readnoise 5.0 --objlim 4.8 \\
		--sigclip 8.0 --fsmode convolve --psfmodel gaussy --psffwhm=2 \\
		$< $@ mask$@

bs_crr_b_%.gz : crr_b_%
	$(BGDSUB) fine.npy $< --gausswidth=100

.PHONY: report finalreport

bias: bias0.1.fits bias2.0.fits $(BIAS)
crrs: $(CRRS)

$(BIAS): bias0.1.fits bias2.0.fits
	$(BSUB) $(subst b_,,$@)

$(CRRS): 
	$(CRRSUB) --niter 4 --sepmed --gain 1.0 --readnoise 5.0 --objlim 4.8 \\
		--sigclip 8.0 --fsmode convolve --psfmodel gaussy --psffwhm=2 \\
		$(subst crr_,,$@) $@ mask$@

calimgs: dome.fits Hg.fits Cd.fits Xe.fits

report:
	$(REPORT) | tee report.txt

ztfupload:
	$(ZTFUPLOAD) $(current_dir)

classify:
	$(CLASS) --specdir $(dir $(mkfile_path))
	$(REPORT) | tee report.txt

finalreport:
	cat report*.txt | mail -s "SEDM DRP Report for $(current_dir)" neill@srl.caltech.edu,rsw@astro.caltech.edu,nblago@caltech.edu,fremling@caltech.edu,ah@astro.caltech.edu,yashuvatsas@gmail.com,tda@lists.astro.caltech.edu

"""


def MF_imcombine(objname, files, dependencies=""):

    filelist = " ".join(["%s " % ifile for ifile in files])
    first = "%s.fits: %s %s\n" % (objname, filelist, dependencies)

    if len(files) > 7:
        reject = "--Nlo 1 --Nhi 1"
    else:
        reject = ""
    if "bias" in objname:
        second = "\t$(IMCOMBINE) --outname %s.fits --listfile %s.lst %s --files %s\n" % (
            objname, objname, reject, filelist)
    else:
        second = "\t$(IMCOMBINE) --outname %s.fits --listfile %s.lst %s --files %s\n" % (
            objname, objname, reject, filelist)

    if "bias" not in objname and "dome" not in objname:
        second += "\n%s.npy : cube.npy %s.fits\n\t$(EXTSINGLE) cube.npy --A %s.fits --outname %s.npy --flat_correction flat-dome-700to900.npy --nosky\n" % (
            objname, objname, objname, objname)

    return first + second + "\n"


def to_makefile(objs, calibs):

    MF = ""

    all = ""

    for calibname, files in calibs.items():

        if "bias" not in calibname:
            pass
        MF += MF_imcombine(calibname, files)
        all += "%s.fits " % calibname

    preamble = make_preamble

    f = open("Makefile", "w")
    clean = "\n\nclean:\n\trm %s" % all

    f.write(preamble + "\nall: %s%s" % (all, clean) + "\n" + MF)
    f.close()


def make_plan(headers):
    """Convert headers to a makefile, assuming headers sorted by JD."""

    objs, calibs = identify_observations(headers)
    to_makefile(objs, calibs)


if __name__ == '__main__':

    files = sys.argv[1:]
    to_process = extract_info(files)

    make_plan(to_process)
