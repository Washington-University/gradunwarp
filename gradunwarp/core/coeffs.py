### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the gradunwarp package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
from __future__ import absolute_import
from collections import namedtuple
import numpy as np
import logging
import re
from . import globals
from .globals import siemens_cas, ge_cas


log = logging.getLogger('gradunwarp')


Coeffs = namedtuple('Coeffs', 'alpha_x, alpha_y, alpha_z, \
                        beta_x, beta_y, beta_z, R0_m')

try:
    advance_iterator = next
except NameError:
    def advance_iterator(it):
        return it.next()
next = advance_iterator


def get_coefficients(vendor, cfile):
    ''' depending on the vendor and the coefficient file,
    return the spherical harmonics coefficients as a named tuple.
    '''
    log.info('Parsing ' + cfile + ' for harmonics coeffs')
    if vendor == 'siemens' and cfile.endswith('.coef'):
        return get_siemens_coef(cfile)
    if vendor == 'siemens' and cfile.endswith('.grad'):
        return get_siemens_grad(cfile)


def coef_file_parse(cfile, txt_var_map):
    ''' a separate function because GE and Siemens .coef files
    have similar structure

    modifies txt_var_map in place
    '''
    # parse .coef file. Strip unneeded characters. a valid line in that file is
    # broken into validline_list
    coef_re = re.compile('^[^\#]')  # regex for first character not a '#'
    coef_file = open(cfile, 'r')
    for line in coef_file.readlines():
        if coef_re.match(line):
            validline_list = line.lstrip(' \t').rstrip(';\n').split()
            if validline_list:
                log.info('Parsed : %s' % validline_list)
                l = validline_list
                x = int(l[1])
                y = int(l[2])
                txt_var_map[l[0]][x, y] = float(l[3])


def get_siemens_coef(cfile):
    ''' Parse the Siemens .coef file.
    Note that R0_m is not explicitly contained in the file
    '''
    R0m_map = {'sonata': 0.25,
               'avanto': 0.25,
               'quantum': 0.25,
               'allegra': 0.14,
               'as39s': 0.25,
               'as39st': 0.25,
               'as39t': 0.25}
    for rad in R0m_map.keys():
        if cfile.startswith(rad):
            R0_m = R0m_map[rad]

    coef_array_sz = siemens_cas
    # allegra is slightly different
    if cfile.startswith('allegra'):
        coef_array_sz = 15

    txt_var_map = create_txt_var_map(coef_array_sz)

    coef_file_parse(cfile, txt_var_map)

    return Coeffs(
        txt_var_map['alpha_y'],
        txt_var_map['alpha_z'],
        txt_var_map['alpha_x'],
        txt_var_map['beta_y'],
        txt_var_map['beta_x'],
        txt_var_map['beta_z'],
        R0_m)


def get_ge_coef(cfile):
    ''' Parse the GE .coef file.
    '''
    txt_var_map = create_txt_var_map(coef_array_sz)

    coef_file_parse(cfile, txt_var_map)

    return Coeffs(
        txt_var_map['alpha_y'],
        txt_var_map['alpha_z'],
        txt_var_map['alpha_x'],
        txt_var_map['beta_y'],
        txt_var_map['beta_x'],
        txt_var_map['beta_z'],
        R0_m)

def grad_file_parse(gfile, txt_var_map):
    xmax = 0
    ymax = 0
    with open(gfile, 'r') as gf:
        for line in gf:
            re_search = re.search("(?P<no>\d+)[\s]+(?P<aorb>[AB])[\s]*\(\s*(?P<x>\d+),\s*(?P<y>\d+)\)\s+(?P<spectrum>[\-]?\d+\.\d+)\s+(?P<axis>[xyz])", line)
            if re_search:
                re_res = re_search.groupdict()
                alphabeta = 'alpha' if re_res['aorb'] == 'A' else 'beta'
                x, y = int(re_res['x']),int(re_res['y'])
                txt_var_map[f"{alphabeta}_{re_res['axis']}"][x, y] = float(re_res['spectrum'])
                xmax, ymax = max(x, xmax), max(y, ymax)
            else:
                re_search = re.search("(?P<R0>\d+\.\d+) m = R0", line)
                if re_search:
                    R0_m = float(re_search.group('R0'))
    return R0_m, (xmax, ymax)

def create_txt_var_map(coef_array_sz):
    return {f'{ab}_{axis}': np.zeros((coef_array_sz,)*2 ) for ab in ['alpha', 'beta'] for axis in 'xyz'}

def get_siemens_grad(gfile):
    ''' Parse the siemens .grad file
    '''
    # allegra is slightly different
    coef_array_sz = 15 if gfile.startswith('coef_AC44') else siemens_cas

    txt_var_map = create_txt_var_map(coef_array_sz)

    R0_m, max_ind = grad_file_parse(gfile, txt_var_map)
    ind = max(max_ind) + 1

    return Coeffs(
        txt_var_map['alpha_x'][:ind, :ind],
        txt_var_map['alpha_y'][:ind, :ind],
        txt_var_map['alpha_z'][:ind, :ind],
        txt_var_map['beta_x'][:ind, :ind],
        txt_var_map['beta_y'][:ind, :ind],
        txt_var_map['beta_z'],
        R0_m)
