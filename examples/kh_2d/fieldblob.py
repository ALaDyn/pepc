#!/usr/bin/env python

import numpy as np

#
# Functions for extracting values from a binary field file.
# For format description see write_quantity_on_grid() in
# src/frontends/pepc-kh/field_helper.f90
#

def FIELDBLOB_HEADER_DTYPE(endianness): return np.dtype([
  ('magic',  '{e}i4'.format(e = endianness), (1,)),
  ('n',      '{e}i4'.format(e = endianness), (2,)),
  ('offset', '{e}f8'.format(e = endianness), (2,)),
  ('extent', '{e}f8'.format(e = endianness), (2,)),
  ('nt'    , '{e}i4'.format(e = endianness), (1,)),
  ('t'     , '{e}f8'.format(e = endianness), (1,)),
  ('B0'    , '{e}f8'.format(e = endianness), (1,)),
  ('vte'   , '{e}f8'.format(e = endianness), (1,)),
  ('vti'   , '{e}f8'.format(e = endianness), (1,)),
  ('qe'    , '{e}f8'.format(e = endianness), (1,)),
  ('qi'    , '{e}f8'.format(e = endianness), (1,)),
  ('me'    , '{e}f8'.format(e = endianness), (1,)),
  ('mi'    , '{e}f8'.format(e = endianness), (1,)),
  ('min'   , '{e}f8'.format(e = endianness), (1,)),
  ('max'   , '{e}f8'.format(e = endianness), (1,))
])

FIELDBLOB_HEADER_SKIP = 512

def endianness_of_fieldblob(fname):
  with open(fname, 'r') as f:
    magic = np.fromfile(f, dtype = np.int8, count = 1)[0]
  
  return '<' if (magic == 1) else '>'

def header_of_fieldblob(fname):
  header_dtype = FIELDBLOB_HEADER_DTYPE(endianness_of_fieldblob(fname))
  with open(fname, 'r') as f:
    header = np.fromfile(f, dtype = header_dtype, count = 1)[0]
  
  return header

def n_of_fieldblob(fname):
  return header_of_fieldblob(fname)['n']

def nx_of_fieldblob(fname):
  return header_of_fieldblob(fname)['n'][0]

def ny_of_fieldblob(fname):
  return header_of_fieldblob(fname)['n'][1]

def offset_of_fieldblob(fname):
  return header_of_fieldblob(fname)['offset']

def offsetx_of_fieldblob(fname):
  return header_of_fieldblob(fname)['offset'][0]

def offsety_of_fieldblob(fname):
  return header_of_fieldblob(fname)['offset'][1]

def extent_of_fieldblob(fname):
  return header_of_fieldblob(fname)['extent']

def extentx_of_fieldblob(fname):
  return header_of_fieldblob(fname)['extent'][0]

def extenty_of_fieldblob(fname):
  return header_of_fieldblob(fname)['extent'][1]

def nt_of_fieldblob(fname):
  return header_of_fieldblob(fname)['nt'][0]

def t_of_fieldblob(fname):
  return header_of_fieldblob(fname)['t'][0]

def B0_of_fieldblob(fname):
  return header_of_fieldblob(fname)['B0'][0]

def vte_of_fieldblob(fname):
  return header_of_fieldblob(fname)['vte'][0]

def vti_of_fieldblob(fname):
  return header_of_fieldblob(fname)['vti'][0]

def qe_of_fieldblob(fname):
  return header_of_fieldblob(fname)['qe'][0]

def qi_of_fieldblob(fname):
  return header_of_fieldblob(fname)['qi'][0]

def me_of_fieldblob(fname):
  return header_of_fieldblob(fname)['me'][0]

def mi_of_fieldblob(fname):
  return header_of_fieldblob(fname)['mi'][0]

def min_of_fieldblob(fname):
  return header_of_fieldblob(fname)['min'][0]

def max_of_fieldblob(fname):
  return header_of_fieldblob(fname)['max'][0]

def field_of_fieldblob(fname):
  n = n_of_fieldblob(fname)
  field_dtype = '{e}f8'.format(e = endianness_of_fieldblob(fname))
  with open(fname, 'r') as f:
    f.seek(FIELDBLOB_HEADER_SKIP)
    y = np.fromfile(f, dtype = field_dtype, count = n[0] * n[1])
  
  y = np.reshape(y, (n[0], n[1]), 'F')
  return y
    
