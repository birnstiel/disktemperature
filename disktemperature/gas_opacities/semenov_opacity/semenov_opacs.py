#!/usr/bin/env python
"""
Python module to use the semenov opacities.
"""

import os
_data_path = os.path.dirname(os.path.abspath(__file__))


def _init_():

    import subprocess
    import glob

    # make sure the files are there

    needed_files = ['opacity.f', 'kP_h2001.dat', 'kR_h2001.dat']
    present_files = [os.path.basename(f) for f in glob.glob(_data_path + os.sep + '*.*')]

    if not all([fi in present_files for fi in needed_files]):
        print('Downloading data files')
        p = subprocess.Popen(['make', 'get'], stderr=subprocess.PIPE,
                             stdout=subprocess.PIPE, cwd=_data_path)
        ret_out, ret_err = p.communicate()
        ret_val = p.poll()
        if ret_val != 0:
            print('return value: {}'.format(ret_val) + '\n' + ret_out + '\n' + ret_err)

    # if the wrapper doesn't exist, strip the program part from the opacities

    if 'opacity_routine.f' not in present_files:
        print('Patching fortran code')
        with open(_data_path + os.sep + 'opacity_routine.f', 'w') as fo, open(_data_path + os.sep + 'opacity.f') as fi:
            write = False
            for line in fi.readlines():
                if write or (write is False and line[0] == 'C'):
                    fo.write(line)
                if line.strip() == 'END':
                    write = True

    # store the header
    header = []
    with open(_data_path + os.sep + 'opacity_routine.f') as fi:
        for line in fi:
            if line[0] != 'C':
                break
            else:
                header += [' ' + line[1:]]
    header = '\n'.join(header)

    # try importing, if it doesn't work, try compiling & then importing

    try:
        from . import opac # noqa
    except:
        print('Compiling fortran module')
        try:
            p = subprocess.Popen('make', stderr=subprocess.PIPE,
                                 stdout=subprocess.PIPE, cwd=_data_path, shell=False)
            ret_out, ret_err = p.communicate()
            ret_val = p.poll()
        except OSError as e:
            ret_out = "OSError({0}): {1}".format(e.errno, e.strerror)
            ret_err = 1
        if ret_val != 0:
            print(ret_err)

    return header


_header = _init_()

from . import opac # noqa
__all__ = ['opac','get_opac'] # noqa


def get_opac(model, top, shape, rosseland, rho, nt, t0, t1):
    """
    Get Opacities from Semenov et al. 2003.

    Arguments:
    ----------

    model : string
        silicate types: 'nrm', 'ips', 'irs'

    top : string
        dust topology: 'c', 'p', 'h'

    shape : string
        dust shape: 's', 'a', '5'

    rosseland : bool
        True:  return Rosseland mean opacities
        False: return Planck mean opacities

    nt : int
        number of temperature points

    t0 : float
        lowest temperature [K]

    t1 : float
         highest temperature [K]

    Output:
    -------
    T,kappa: temperature and opacity arrays
    """
    #
    # this is to circumvent what is likely a C-compiler error
    #
    import locale
    locale.setlocale(locale.LC_ALL, 'C')

    cwd = os.getcwd()
    kappa = None
    try:
        os.chdir(_data_path)
        kappa = opac.opacity(model, top, shape, rosseland, rho, nt, t0, t1)
    except:
        pass
    finally:
        os.chdir(cwd)
    return kappa


get_opac.__doc__ += '\nOriginal Header with more description:\n' + 38 * '=' + '\n' + _header


def _test(plot=True):
    import numpy as np
    with open(_data_path + os.sep + 'opacity.inp') as f:
        model = f.readline().split()[0]
        top = f.readline().split()[0]
        shape = f.readline().split()[0]
        ross = bool(f.readline().split()[0])
        rho = float(f.readline().split()[0])
        nt = int(f.readline().split()[0])
        t0 = float(f.readline().split()[0])
        t1 = float(f.readline().split()[0])
    k_orig = np.loadtxt(_data_path + os.sep + 'k' + ('R' * ross) + ('P' * (not ross)) + '.out').T
    k_test = get_opac(model, top, shape, ross, rho, nt, t0, t1)
    k_test = np.array([k_test[0], k_test[1]])
    if plot:
        import matplotlib.pyplot as plt
        from matplotlib.pyplot import style
        style.use(['seaborn-darkgrid', {'font.size': 15}])
        f, ax = plt.subplots()
        ax.loglog(*k_orig, label='original data')
        ax.loglog(*k_test, label='test data')
        ax.set_xlabel('T [K]')
        ax.legend(fontsize='small', loc='best')
        ax.set_ylabel('$\kappa_\mathrm{' + (ross * 'R' + (not ross) * 'P') + '}$ [cm$^2$ / g]')
    print('Test was ' + 'not ' * (not np.all(np.isclose(k_orig, k_test))) + 'passed')
