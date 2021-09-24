#!/usr/bin/env python

# This script parses spglib_f.c, searches for all symbols, and writes a new
# wrapper spglib_f_meta.c with all combinations of name mangling.
# This script should only be called when spglib is updated.

# Felipe H. da Jornada, Nov 2013

import re

f_in = open('spglib_f.c')
str_in = f_in.read()
f_in.close()

f_out = open('spglib_f_meta.c', 'w')
f_out.write(
'''\
/******************************************************************************
*
*  spglib_f_meta           Originally by FHJ     Last Modified 11/08/2012 (FHJ)
*
*
*  SHORT
*    BerkeleyGW meta wrapper for spglib.
*
*  LONG
*    spglib`s Fortran wrappers assume the Fortran compiler will append an
*    underscore after each symbol`s name, which is not true for all compilers,
*    such as XL. This can lead to terrible results, especially since the symbols
*    without underscores are also available in the library, but take different
*    arguments. This file fixes this problem by exporting symbols with both
*    "_f" and "_f_" appended to the names, so that they should work regardless
*    of the name mangling.
*
*    This file is part of the BerkeleyGW package.
*
******************************************************************************/

''')

for match in re.findall(r'^void spg_[^)]*\);', str_in, re.DOTALL+re.MULTILINE):
    obj = re.search(r'(spg_[^\s(]*)(\([^)]*\))', match)
    try:
        func_name = obj.group(1)
        arg_list = obj.group(2)
    except:
        continue
    arg_list_stripped = arg_list.replace('*', '').replace('&', '')
    subs = ['short','long','unsigned','signed','char','int','float','double','void','const']
    for sub in subs:
        arg_list_stripped = re.sub(r'\W'+sub+'\W', '', arg_list_stripped)
    arg_list_stripped = re.sub('\[[^]]*\]', '', arg_list_stripped)
    print('Exporting symbol %s'%(func_name))
    f_out.write('\n/* Original symbol from spglib_f.c */\n')
    f_out.write('extern %s\n'%(match))
    f_out.write('/* Wrapped symbol without underscore */\n')
    f_out.write('''void %s%s {
    %s%s;
}
'''%(func_name+'f', arg_list, func_name, arg_list_stripped))
    f_out.write('/* Wrapped symbol with underscore */\n')
    f_out.write('''void %s%s {
    %s%s;
}
'''%(func_name+'f_', arg_list, func_name, arg_list_stripped))
