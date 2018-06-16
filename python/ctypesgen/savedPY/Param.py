'''Wrapper for Param.h

Generated with:
ctypesgen.py -lneon /home/xzhang11/Planet/Libs/ParaDiSDll/include/Param.h -o Param.py

Do not modify this file.
'''

__docformat__ =  'restructuredtext'

# Begin preamble

import ctypes, os, sys
from ctypes import *

_int_types = (c_int16, c_int32)
if hasattr(ctypes, 'c_int64'):
    # Some builds of ctypes apparently do not have c_int64
    # defined; it's a pretty good bet that these builds do not
    # have 64-bit pointers.
    _int_types += (c_int64,)
for t in _int_types:
    if sizeof(t) == sizeof(c_size_t):
        c_ptrdiff_t = t
del t
del _int_types

class c_void(Structure):
    # c_void_p is a buggy return type, converting to int, so
    # POINTER(None) == c_void_p is actually written as
    # POINTER(c_void), so it can be treated as a real pointer.
    _fields_ = [('dummy', c_int)]

def POINTER(obj):
    p = ctypes.POINTER(obj)

    # Convert None to a real NULL pointer to work around bugs
    # in how ctypes handles None on 64-bit platforms
    if not isinstance(p.from_param, classmethod):
        def from_param(cls, x):
            if x is None:
                return cls()
            else:
                return x
        p.from_param = classmethod(from_param)

    return p

class UserString:
    def __init__(self, seq):
        if isinstance(seq, basestring):
            self.data = seq
        elif isinstance(seq, UserString):
            self.data = seq.data[:]
        else:
            self.data = str(seq)
    def __str__(self): return str(self.data)
    def __repr__(self): return repr(self.data)
    def __int__(self): return int(self.data)
    def __long__(self): return long(self.data)
    def __float__(self): return float(self.data)
    def __complex__(self): return complex(self.data)
    def __hash__(self): return hash(self.data)

    def __cmp__(self, string):
        if isinstance(string, UserString):
            return cmp(self.data, string.data)
        else:
            return cmp(self.data, string)
    def __contains__(self, char):
        return char in self.data

    def __len__(self): return len(self.data)
    def __getitem__(self, index): return self.__class__(self.data[index])
    def __getslice__(self, start, end):
        start = max(start, 0); end = max(end, 0)
        return self.__class__(self.data[start:end])

    def __add__(self, other):
        if isinstance(other, UserString):
            return self.__class__(self.data + other.data)
        elif isinstance(other, basestring):
            return self.__class__(self.data + other)
        else:
            return self.__class__(self.data + str(other))
    def __radd__(self, other):
        if isinstance(other, basestring):
            return self.__class__(other + self.data)
        else:
            return self.__class__(str(other) + self.data)
    def __mul__(self, n):
        return self.__class__(self.data*n)
    __rmul__ = __mul__
    def __mod__(self, args):
        return self.__class__(self.data % args)

    # the following methods are defined in alphabetical order:
    def capitalize(self): return self.__class__(self.data.capitalize())
    def center(self, width, *args):
        return self.__class__(self.data.center(width, *args))
    def count(self, sub, start=0, end=sys.maxint):
        return self.data.count(sub, start, end)
    def decode(self, encoding=None, errors=None): # XXX improve this?
        if encoding:
            if errors:
                return self.__class__(self.data.decode(encoding, errors))
            else:
                return self.__class__(self.data.decode(encoding))
        else:
            return self.__class__(self.data.decode())
    def encode(self, encoding=None, errors=None): # XXX improve this?
        if encoding:
            if errors:
                return self.__class__(self.data.encode(encoding, errors))
            else:
                return self.__class__(self.data.encode(encoding))
        else:
            return self.__class__(self.data.encode())
    def endswith(self, suffix, start=0, end=sys.maxint):
        return self.data.endswith(suffix, start, end)
    def expandtabs(self, tabsize=8):
        return self.__class__(self.data.expandtabs(tabsize))
    def find(self, sub, start=0, end=sys.maxint):
        return self.data.find(sub, start, end)
    def index(self, sub, start=0, end=sys.maxint):
        return self.data.index(sub, start, end)
    def isalpha(self): return self.data.isalpha()
    def isalnum(self): return self.data.isalnum()
    def isdecimal(self): return self.data.isdecimal()
    def isdigit(self): return self.data.isdigit()
    def islower(self): return self.data.islower()
    def isnumeric(self): return self.data.isnumeric()
    def isspace(self): return self.data.isspace()
    def istitle(self): return self.data.istitle()
    def isupper(self): return self.data.isupper()
    def join(self, seq): return self.data.join(seq)
    def ljust(self, width, *args):
        return self.__class__(self.data.ljust(width, *args))
    def lower(self): return self.__class__(self.data.lower())
    def lstrip(self, chars=None): return self.__class__(self.data.lstrip(chars))
    def partition(self, sep):
        return self.data.partition(sep)
    def replace(self, old, new, maxsplit=-1):
        return self.__class__(self.data.replace(old, new, maxsplit))
    def rfind(self, sub, start=0, end=sys.maxint):
        return self.data.rfind(sub, start, end)
    def rindex(self, sub, start=0, end=sys.maxint):
        return self.data.rindex(sub, start, end)
    def rjust(self, width, *args):
        return self.__class__(self.data.rjust(width, *args))
    def rpartition(self, sep):
        return self.data.rpartition(sep)
    def rstrip(self, chars=None): return self.__class__(self.data.rstrip(chars))
    def split(self, sep=None, maxsplit=-1):
        return self.data.split(sep, maxsplit)
    def rsplit(self, sep=None, maxsplit=-1):
        return self.data.rsplit(sep, maxsplit)
    def splitlines(self, keepends=0): return self.data.splitlines(keepends)
    def startswith(self, prefix, start=0, end=sys.maxint):
        return self.data.startswith(prefix, start, end)
    def strip(self, chars=None): return self.__class__(self.data.strip(chars))
    def swapcase(self): return self.__class__(self.data.swapcase())
    def title(self): return self.__class__(self.data.title())
    def translate(self, *args):
        return self.__class__(self.data.translate(*args))
    def upper(self): return self.__class__(self.data.upper())
    def zfill(self, width): return self.__class__(self.data.zfill(width))

class MutableString(UserString):
    """mutable string objects

    Python strings are immutable objects.  This has the advantage, that
    strings may be used as dictionary keys.  If this property isn't needed
    and you insist on changing string values in place instead, you may cheat
    and use MutableString.

    But the purpose of this class is an educational one: to prevent
    people from inventing their own mutable string class derived
    from UserString and than forget thereby to remove (override) the
    __hash__ method inherited from UserString.  This would lead to
    errors that would be very hard to track down.

    A faster and better solution is to rewrite your program using lists."""
    def __init__(self, string=""):
        self.data = string
    def __hash__(self):
        raise TypeError("unhashable type (it is mutable)")
    def __setitem__(self, index, sub):
        if index < 0:
            index += len(self.data)
        if index < 0 or index >= len(self.data): raise IndexError
        self.data = self.data[:index] + sub + self.data[index+1:]
    def __delitem__(self, index):
        if index < 0:
            index += len(self.data)
        if index < 0 or index >= len(self.data): raise IndexError
        self.data = self.data[:index] + self.data[index+1:]
    def __setslice__(self, start, end, sub):
        start = max(start, 0); end = max(end, 0)
        if isinstance(sub, UserString):
            self.data = self.data[:start]+sub.data+self.data[end:]
        elif isinstance(sub, basestring):
            self.data = self.data[:start]+sub+self.data[end:]
        else:
            self.data =  self.data[:start]+str(sub)+self.data[end:]
    def __delslice__(self, start, end):
        start = max(start, 0); end = max(end, 0)
        self.data = self.data[:start] + self.data[end:]
    def immutable(self):
        return UserString(self.data)
    def __iadd__(self, other):
        if isinstance(other, UserString):
            self.data += other.data
        elif isinstance(other, basestring):
            self.data += other
        else:
            self.data += str(other)
        return self
    def __imul__(self, n):
        self.data *= n
        return self

class String(MutableString, Union):

    _fields_ = [('raw', POINTER(c_char)),
                ('data', c_char_p)]

    def __init__(self, obj=""):
        if isinstance(obj, (str, unicode, UserString)):
            self.data = str(obj)
        else:
            self.raw = obj

    def __len__(self):
        return self.data and len(self.data) or 0

    def from_param(cls, obj):
        # Convert None or 0
        if obj is None or obj == 0:
            return cls(POINTER(c_char)())

        # Convert from String
        elif isinstance(obj, String):
            return obj

        # Convert from str
        elif isinstance(obj, str):
            return cls(obj)

        # Convert from c_char_p
        elif isinstance(obj, c_char_p):
            return obj

        # Convert from POINTER(c_char)
        elif isinstance(obj, POINTER(c_char)):
            return obj

        # Convert from raw pointer
        elif isinstance(obj, int):
            return cls(cast(obj, POINTER(c_char)))

        # Convert from object
        else:
            return String.from_param(obj._as_parameter_)
    from_param = classmethod(from_param)

def ReturnString(obj, func=None, arguments=None):
    return String.from_param(obj)

# As of ctypes 1.0, ctypes does not support custom error-checking
# functions on callbacks, nor does it support custom datatypes on
# callbacks, so we must ensure that all callbacks return
# primitive datatypes.
#
# Non-primitive return values wrapped with UNCHECKED won't be
# typechecked, and will be converted to c_void_p.
def UNCHECKED(type):
    if (hasattr(type, "_type_") and isinstance(type._type_, str)
        and type._type_ != "P"):
        return type
    else:
        return c_void_p

# ctypes doesn't have direct support for variadic functions, so we have to write
# our own wrapper class
class _variadic_function(object):
    def __init__(self,func,restype,argtypes):
        self.func=func
        self.func.restype=restype
        self.argtypes=argtypes
    def _as_parameter_(self):
        # So we can pass this variadic function as a function pointer
        return self.func
    def __call__(self,*args):
        fixed_args=[]
        i=0
        for argtype in self.argtypes:
            # Typecheck what we can
            fixed_args.append(argtype.from_param(args[i]))
            i+=1
        return self.func(*fixed_args+list(args[i:]))

# End preamble

_libs = {}
_libdirs = []

# Begin loader

# ----------------------------------------------------------------------------
# Copyright (c) 2008 David James
# Copyright (c) 2006-2008 Alex Holkner
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
#  * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in
#    the documentation and/or other materials provided with the
#    distribution.
#  * Neither the name of pyglet nor the names of its
#    contributors may be used to endorse or promote products
#    derived from this software without specific prior written
#    permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
# ----------------------------------------------------------------------------

import os.path, re, sys, glob
import platform
import ctypes
import ctypes.util

def _environ_path(name):
    if name in os.environ:
        return os.environ[name].split(":")
    else:
        return []

class LibraryLoader(object):
    def __init__(self):
        self.other_dirs=[]

    def load_library(self,libname):
        """Given the name of a library, load it."""
        paths = self.getpaths(libname)

        for path in paths:
            if os.path.exists(path):
                return self.load(path)

        raise ImportError("%s not found." % libname)

    def load(self,path):
        """Given a path to a library, load it."""
        try:
            # Darwin requires dlopen to be called with mode RTLD_GLOBAL instead
            # of the default RTLD_LOCAL.  Without this, you end up with
            # libraries not being loadable, resulting in "Symbol not found"
            # errors
            if sys.platform == 'darwin':
                return ctypes.CDLL(path, ctypes.RTLD_GLOBAL)
            else:
                return ctypes.cdll.LoadLibrary(path)
        except OSError,e:
            raise ImportError(e)

    def getpaths(self,libname):
        """Return a list of paths where the library might be found."""
        if os.path.isabs(libname):
            yield libname
        else:
            # FIXME / TODO return '.' and os.path.dirname(__file__)
            for path in self.getplatformpaths(libname):
                yield path

            path = ctypes.util.find_library(libname)
            if path: yield path

    def getplatformpaths(self, libname):
        return []

# Darwin (Mac OS X)

class DarwinLibraryLoader(LibraryLoader):
    name_formats = ["lib%s.dylib", "lib%s.so", "lib%s.bundle", "%s.dylib",
                "%s.so", "%s.bundle", "%s"]

    def getplatformpaths(self,libname):
        if os.path.pathsep in libname:
            names = [libname]
        else:
            names = [format % libname for format in self.name_formats]

        for dir in self.getdirs(libname):
            for name in names:
                yield os.path.join(dir,name)

    def getdirs(self,libname):
        '''Implements the dylib search as specified in Apple documentation:

        http://developer.apple.com/documentation/DeveloperTools/Conceptual/
            DynamicLibraries/Articles/DynamicLibraryUsageGuidelines.html

        Before commencing the standard search, the method first checks
        the bundle's ``Frameworks`` directory if the application is running
        within a bundle (OS X .app).
        '''

        dyld_fallback_library_path = _environ_path("DYLD_FALLBACK_LIBRARY_PATH")
        if not dyld_fallback_library_path:
            dyld_fallback_library_path = [os.path.expanduser('~/lib'),
                                          '/usr/local/lib', '/usr/lib']

        dirs = []

        if '/' in libname:
            dirs.extend(_environ_path("DYLD_LIBRARY_PATH"))
        else:
            dirs.extend(_environ_path("LD_LIBRARY_PATH"))
            dirs.extend(_environ_path("DYLD_LIBRARY_PATH"))

        dirs.extend(self.other_dirs)
        dirs.append(".")
        dirs.append(os.path.dirname(__file__))

        if hasattr(sys, 'frozen') and sys.frozen == 'macosx_app':
            dirs.append(os.path.join(
                os.environ['RESOURCEPATH'],
                '..',
                'Frameworks'))

        dirs.extend(dyld_fallback_library_path)

        return dirs

# Posix

class PosixLibraryLoader(LibraryLoader):
    _ld_so_cache = None

    def _create_ld_so_cache(self):
        # Recreate search path followed by ld.so.  This is going to be
        # slow to build, and incorrect (ld.so uses ld.so.cache, which may
        # not be up-to-date).  Used only as fallback for distros without
        # /sbin/ldconfig.
        #
        # We assume the DT_RPATH and DT_RUNPATH binary sections are omitted.

        directories = []
        for name in ("LD_LIBRARY_PATH",
                     "SHLIB_PATH", # HPUX
                     "LIBPATH", # OS/2, AIX
                     "LIBRARY_PATH", # BE/OS
                    ):
            if name in os.environ:
                directories.extend(os.environ[name].split(os.pathsep))
        directories.extend(self.other_dirs)
        directories.append(".")
        directories.append(os.path.dirname(__file__))

        try: directories.extend([dir.strip() for dir in open('/etc/ld.so.conf')])
        except IOError: pass

        unix_lib_dirs_list = ['/lib', '/usr/lib', '/lib64', '/usr/lib64']
        if sys.platform.startswith('linux'):
            # Try and support multiarch work in Ubuntu
            # https://wiki.ubuntu.com/MultiarchSpec
            bitage = platform.architecture()[0]
            if bitage.startswith('32'):
                # Assume Intel/AMD x86 compat
                unix_lib_dirs_list += ['/lib/i386-linux-gnu', '/usr/lib/i386-linux-gnu']
            elif bitage.startswith('64'):
                # Assume Intel/AMD x86 compat
                unix_lib_dirs_list += ['/lib/x86_64-linux-gnu', '/usr/lib/x86_64-linux-gnu']
            else:
                # guess...
                unix_lib_dirs_list += glob.glob('/lib/*linux-gnu')
        directories.extend(unix_lib_dirs_list)

        cache = {}
        lib_re = re.compile(r'lib(.*)\.s[ol]')
        ext_re = re.compile(r'\.s[ol]$')
        for dir in directories:
            try:
                for path in glob.glob("%s/*.s[ol]*" % dir):
                    file = os.path.basename(path)

                    # Index by filename
                    if file not in cache:
                        cache[file] = path

                    # Index by library name
                    match = lib_re.match(file)
                    if match:
                        library = match.group(1)
                        if library not in cache:
                            cache[library] = path
            except OSError:
                pass

        self._ld_so_cache = cache

    def getplatformpaths(self, libname):
        if self._ld_so_cache is None:
            self._create_ld_so_cache()

        result = self._ld_so_cache.get(libname)
        if result: yield result

        path = ctypes.util.find_library(libname)
        if path: yield os.path.join("/lib",path)

# Windows

class _WindowsLibrary(object):
    def __init__(self, path):
        self.cdll = ctypes.cdll.LoadLibrary(path)
        self.windll = ctypes.windll.LoadLibrary(path)

    def __getattr__(self, name):
        try: return getattr(self.cdll,name)
        except AttributeError:
            try: return getattr(self.windll,name)
            except AttributeError:
                raise

class WindowsLibraryLoader(LibraryLoader):
    name_formats = ["%s.dll", "lib%s.dll", "%slib.dll"]

    def load_library(self, libname):
        try:
            result = LibraryLoader.load_library(self, libname)
        except ImportError:
            result = None
            if os.path.sep not in libname:
                for name in self.name_formats:
                    try:
                        result = getattr(ctypes.cdll, name % libname)
                        if result:
                            break
                    except WindowsError:
                        result = None
            if result is None:
                try:
                    result = getattr(ctypes.cdll, libname)
                except WindowsError:
                    result = None
            if result is None:
                raise ImportError("%s not found." % libname)
        return result

    def load(self, path):
        return _WindowsLibrary(path)

    def getplatformpaths(self, libname):
        if os.path.sep not in libname:
            for name in self.name_formats:
                dll_in_current_dir = os.path.abspath(name % libname)
                if os.path.exists(dll_in_current_dir):
                    yield dll_in_current_dir
                path = ctypes.util.find_library(name % libname)
                if path:
                    yield path

# Platform switching

# If your value of sys.platform does not appear in this dict, please contact
# the Ctypesgen maintainers.

loaderclass = {
    "darwin":   DarwinLibraryLoader,
    "cygwin":   WindowsLibraryLoader,
    "win32":    WindowsLibraryLoader
}

loader = loaderclass.get(sys.platform, PosixLibraryLoader)()

def add_library_search_dirs(other_dirs):
    loader.other_dirs = other_dirs

load_library = loader.load_library

del loaderclass

# End loader

add_library_search_dirs([])

# Begin libraries

_libs["neon"] = load_library("neon")

# 1 libraries
# End libraries

# No modules

real8 = c_double # /home/xzhang11/Planet/Libs/ParaDiSDll/include/Typedefs.h: 17

# /home/xzhang11/Planet/Libs/ParaDiSDll/include/Cell.h: 15
class struct__cell(Structure):
    pass

Cell_t = struct__cell # /home/xzhang11/Planet/Libs/ParaDiSDll/include/Typedefs.h: 21

# /home/xzhang11/Planet/Libs/ParaDiSDll/include/Home.h: 161
class struct__home(Structure):
    pass

Home_t = struct__home # /home/xzhang11/Planet/Libs/ParaDiSDll/include/Typedefs.h: 23

# /home/xzhang11/Planet/Libs/ParaDiSDll/include/MirrorDomain.h: 17
class struct__mirrordomain(Structure):
    pass

MirrorDomain_t = struct__mirrordomain # /home/xzhang11/Planet/Libs/ParaDiSDll/include/Typedefs.h: 26

# /home/xzhang11/Planet/Libs/ParaDiSDll/include/Node.h: 46
class struct__node(Structure):
    pass

Node_t = struct__node # /home/xzhang11/Planet/Libs/ParaDiSDll/include/Typedefs.h: 27

# /home/xzhang11/Planet/Libs/ParaDiSDll/include/Node.h: 127
class struct__nodeblock(Structure):
    pass

NodeBlock_t = struct__nodeblock # /home/xzhang11/Planet/Libs/ParaDiSDll/include/Typedefs.h: 28

# /home/xzhang11/Planet/Libs/ParaDiSDll/include/OpList.h: 18
class struct__operate(Structure):
    pass

Operate_t = struct__operate # /home/xzhang11/Planet/Libs/ParaDiSDll/include/Typedefs.h: 29

# /home/xzhang11/Planet/Libs/ParaDiSDll/include/Param.h: 25
class struct__param(Structure):
    pass

Param_t = struct__param # /home/xzhang11/Planet/Libs/ParaDiSDll/include/Typedefs.h: 30

# /home/xzhang11/Planet/Libs/ParaDiSDll/include/RemoteDomain.h: 14
class struct__remotedomain(Structure):
    pass

RemoteDomain_t = struct__remotedomain # /home/xzhang11/Planet/Libs/ParaDiSDll/include/Typedefs.h: 31

# /home/xzhang11/Planet/Libs/ParaDiSDll/include/Tag.h: 12
class struct__tag(Structure):
    pass

Tag_t = struct__tag # /home/xzhang11/Planet/Libs/ParaDiSDll/include/Typedefs.h: 33

# /home/xzhang11/Planet/Libs/ParaDiSDll/include/Timer.h: 12
class struct__timer(Structure):
    pass

Timer_t = struct__timer # /home/xzhang11/Planet/Libs/ParaDiSDll/include/Typedefs.h: 34

enum_anon_32 = c_int # /home/xzhang11/Planet/Libs/ParaDiSDll/include/Typedefs.h: 41

BoundType_t = enum_anon_32 # /home/xzhang11/Planet/Libs/ParaDiSDll/include/Typedefs.h: 41

enum_anon_33 = c_int # /home/xzhang11/Planet/Libs/ParaDiSDll/include/Typedefs.h: 58

OpType_t = enum_anon_33 # /home/xzhang11/Planet/Libs/ParaDiSDll/include/Typedefs.h: 58

# /home/xzhang11/Planet/Libs/ParaDiSDll/include/Typedefs.h: 69
class struct_anon_34(Structure):
    pass

struct_anon_34.__slots__ = [
    'node',
    'next',
]
struct_anon_34._fields_ = [
    ('node', POINTER(Node_t)),
    ('next', c_int),
]

C2Qent_t = struct_anon_34 # /home/xzhang11/Planet/Libs/ParaDiSDll/include/Typedefs.h: 69

# /home/xzhang11/Planet/Libs/ParaDiSDll/include/Typedefs.h: 101
class struct_anon_36(Structure):
    pass

struct_anon_36.__slots__ = [
    'varName',
    'valType',
    'valCnt',
    'flags',
    'valList',
]
struct_anon_36._fields_ = [
    ('varName', c_char * 256),
    ('valType', c_int),
    ('valCnt', c_int),
    ('flags', c_int),
    ('valList', POINTER(None)),
]

VarData_t = struct_anon_36 # /home/xzhang11/Planet/Libs/ParaDiSDll/include/Typedefs.h: 101

# /home/xzhang11/Planet/Libs/ParaDiSDll/include/Typedefs.h: 107
class struct_anon_37(Structure):
    pass

struct_anon_37.__slots__ = [
    'paramCnt',
    'varList',
]
struct_anon_37._fields_ = [
    ('paramCnt', c_int),
    ('varList', POINTER(VarData_t)),
]

ParamList_t = struct_anon_37 # /home/xzhang11/Planet/Libs/ParaDiSDll/include/Typedefs.h: 107

struct__tag.__slots__ = [
    'domainID',
    'index',
]
struct__tag._fields_ = [
    ('domainID', c_int),
    ('index', c_int),
]

# /home/xzhang11/Planet/Libs/ParaDiSDll/include/FM.h: 42
class struct__fmcell(Structure):
    pass

FMCell_t = struct__fmcell # /home/xzhang11/Planet/Libs/ParaDiSDll/include/FM.h: 39

# /home/xzhang11/Planet/Libs/ParaDiSDll/include/FM.h: 71
class struct__fmlayer(Structure):
    pass

FMLayer_t = struct__fmlayer # /home/xzhang11/Planet/Libs/ParaDiSDll/include/FM.h: 40

struct__fmcell.__slots__ = [
    'cellID',
    'owningDom',
    'domCnt',
    'domList',
    'cellCtr',
    'mpCoeff',
    'taylorCoeff',
    'next',
    'prev',
]
struct__fmcell._fields_ = [
    ('cellID', c_int),
    ('owningDom', c_int),
    ('domCnt', c_int),
    ('domList', POINTER(c_int)),
    ('cellCtr', real8 * 3),
    ('mpCoeff', POINTER(real8)),
    ('taylorCoeff', POINTER(real8)),
    ('next', POINTER(FMCell_t)),
    ('prev', POINTER(FMCell_t)),
]

struct__fmlayer.__slots__ = [
    'lDim',
    'cellSize',
    'ownedCnt',
    'ownedMin',
    'ownedMax',
    'intersectCnt',
    'intersectMin',
    'intersectMax',
    'domBuf',
    'fmUpPassSendDomCnt',
    'fmUpPassSendDomList',
    'fmUpPassRecvDomCnt',
    'fmUpPassRecvDomList',
    'fmDownPassSendDomCnt',
    'fmDownPassSendDomList',
    'fmDownPassRecvDomCnt',
    'fmDownPassRecvDomList',
    'cellList',
    'numCells',
    'cellTable',
]
struct__fmlayer._fields_ = [
    ('lDim', c_int * 3),
    ('cellSize', real8 * 3),
    ('ownedCnt', c_int),
    ('ownedMin', c_int * 3),
    ('ownedMax', c_int * 3),
    ('intersectCnt', c_int),
    ('intersectMin', c_int * 3),
    ('intersectMax', c_int * 3),
    ('domBuf', POINTER(c_int)),
    ('fmUpPassSendDomCnt', c_int),
    ('fmUpPassSendDomList', POINTER(c_int)),
    ('fmUpPassRecvDomCnt', c_int),
    ('fmUpPassRecvDomList', POINTER(c_int)),
    ('fmDownPassSendDomCnt', c_int),
    ('fmDownPassSendDomList', POINTER(c_int)),
    ('fmDownPassRecvDomCnt', c_int),
    ('fmDownPassRecvDomList', POINTER(c_int)),
    ('cellList', POINTER(c_int)),
    ('numCells', c_int),
    ('cellTable', POINTER(FMCell_t) * 97),
]

struct__node.__slots__ = [
    'flags',
    'x',
    'y',
    'z',
    'fX',
    'fY',
    'fZ',
    'vX',
    'vY',
    'vZ',
    'oldx',
    'oldy',
    'oldz',
    'oldvX',
    'oldvY',
    'oldvZ',
    'currvX',
    'currvY',
    'currvZ',
    'myTag',
    'numNbrs',
    'nbrTag',
    'armfx',
    'armfy',
    'armfz',
    'burgX',
    'burgY',
    'burgZ',
    'nx',
    'ny',
    'nz',
    'sigbLoc',
    'sigbRem',
    'armCoordIndex',
    'constraint',
    'cellIdx',
    'cell2Idx',
    'cell2QentIdx',
    'native',
    'next',
    'nextInCell',
    'sgnv',
    'fem_Surface',
    'fem_Surface_Norm',
]
struct__node._fields_ = [
    ('flags', c_int),
    ('x', real8),
    ('y', real8),
    ('z', real8),
    ('fX', real8),
    ('fY', real8),
    ('fZ', real8),
    ('vX', real8),
    ('vY', real8),
    ('vZ', real8),
    ('oldx', real8),
    ('oldy', real8),
    ('oldz', real8),
    ('oldvX', real8),
    ('oldvY', real8),
    ('oldvZ', real8),
    ('currvX', real8),
    ('currvY', real8),
    ('currvZ', real8),
    ('myTag', Tag_t),
    ('numNbrs', c_int),
    ('nbrTag', POINTER(Tag_t)),
    ('armfx', POINTER(real8)),
    ('armfy', POINTER(real8)),
    ('armfz', POINTER(real8)),
    ('burgX', POINTER(real8)),
    ('burgY', POINTER(real8)),
    ('burgZ', POINTER(real8)),
    ('nx', POINTER(real8)),
    ('ny', POINTER(real8)),
    ('nz', POINTER(real8)),
    ('sigbLoc', POINTER(real8)),
    ('sigbRem', POINTER(real8)),
    ('armCoordIndex', POINTER(c_int)),
    ('constraint', c_int),
    ('cellIdx', c_int),
    ('cell2Idx', c_int),
    ('cell2QentIdx', c_int),
    ('native', c_int),
    ('next', POINTER(Node_t)),
    ('nextInCell', POINTER(Node_t)),
    ('sgnv', c_int),
    ('fem_Surface', c_int * 2),
    ('fem_Surface_Norm', real8 * 3),
]

struct__nodeblock.__slots__ = [
    'next',
    'nodes',
]
struct__nodeblock._fields_ = [
    ('next', POINTER(NodeBlock_t)),
    ('nodes', POINTER(Node_t)),
]

struct__cell.__slots__ = [
    'nodeQ',
    'nodeCount',
    'nbrList',
    'nbrCount',
    'domains',
    'domCount',
    'baseIdx',
    'xShift',
    'yShift',
    'zShift',
    'xIndex',
    'yIndex',
    'zIndex',
]
struct__cell._fields_ = [
    ('nodeQ', POINTER(Node_t)),
    ('nodeCount', c_int),
    ('nbrList', POINTER(c_int)),
    ('nbrCount', c_int),
    ('domains', POINTER(c_int)),
    ('domCount', c_int),
    ('baseIdx', c_int),
    ('xShift', real8),
    ('yShift', real8),
    ('zShift', real8),
    ('xIndex', c_int),
    ('yIndex', c_int),
    ('zIndex', c_int),
]

struct__remotedomain.__slots__ = [
    'domainIdx',
    'numExpCells',
    'expCells',
    'maxTagIndex',
    'nodeKeys',
    'inBufLen',
    'inBuf',
    'outBufLen',
    'outBuf',
]
struct__remotedomain._fields_ = [
    ('domainIdx', c_int),
    ('numExpCells', c_int),
    ('expCells', POINTER(c_int)),
    ('maxTagIndex', c_int),
    ('nodeKeys', POINTER(POINTER(Node_t))),
    ('inBufLen', c_int),
    ('inBuf', String),
    ('outBufLen', c_int),
    ('outBuf', String),
]

struct__mirrordomain.__slots__ = [
    'nodeKeys',
    'newNodeKeyPtr',
    'armX',
    'armY',
    'armZ',
]
struct__mirrordomain._fields_ = [
    ('nodeKeys', POINTER(POINTER(Node_t))),
    ('newNodeKeyPtr', c_int),
    ('armX', POINTER(real8)),
    ('armY', POINTER(real8)),
    ('armZ', POINTER(real8)),
]

struct__operate.__slots__ = [
    'type',
    'dom1',
    'idx1',
    'dom2',
    'idx2',
    'dom3',
    'idx3',
    'bx',
    'by',
    'bz',
    'x',
    'y',
    'z',
    'nx',
    'ny',
    'nz',
]
struct__operate._fields_ = [
    ('type', OpType_t),
    ('dom1', c_int),
    ('idx1', c_int),
    ('dom2', c_int),
    ('idx2', c_int),
    ('dom3', c_int),
    ('idx3', c_int),
    ('bx', real8),
    ('by', real8),
    ('bz', real8),
    ('x', real8),
    ('y', real8),
    ('z', real8),
    ('nx', real8),
    ('ny', real8),
    ('nz', real8),
]

struct__timer.__slots__ = [
    'startTime',
    'incr',
    'accum',
    'save',
    'started',
    'name',
]
struct__timer._fields_ = [
    ('startTime', real8),
    ('incr', real8),
    ('accum', real8),
    ('save', real8),
    ('started', c_int),
    ('name', String),
]

MPI_Comm = c_int # /home/xzhang11/Planet/Libs/ParaDiSDll/include/mpi.h: 291

MPI_Request = c_int # /home/xzhang11/Planet/Libs/ParaDiSDll/include/mpi.h: 455

# /home/xzhang11/Planet/Libs/ParaDiSDll/include/mpi.h: 602
class struct_MPI_Status(Structure):
    pass

struct_MPI_Status.__slots__ = [
    'count_lo',
    'count_hi_and_cancelled',
    'MPI_SOURCE',
    'MPI_TAG',
    'MPI_ERROR',
]
struct_MPI_Status._fields_ = [
    ('count_lo', c_int),
    ('count_hi_and_cancelled', c_int),
    ('MPI_SOURCE', c_int),
    ('MPI_TAG', c_int),
    ('MPI_ERROR', c_int),
]

MPI_Status = struct_MPI_Status # /home/xzhang11/Planet/Libs/ParaDiSDll/include/mpi.h: 602

# /home/xzhang11/Planet/Libs/ParaDiSDll/include/Home.h: 119
class struct_anon_39(Structure):
    pass

struct_anon_39.__slots__ = [
    'oldTag',
    'newTag',
]
struct_anon_39._fields_ = [
    ('oldTag', Tag_t),
    ('newTag', Tag_t),
]

TagMap_t = struct_anon_39 # /home/xzhang11/Planet/Libs/ParaDiSDll/include/Home.h: 119

# /home/xzhang11/Planet/Libs/ParaDiSDll/include/Home.h: 156
class struct_anon_41(Structure):
    pass

struct_anon_41.__slots__ = [
    'numBurgVectors',
    'numPlanes',
    'numPlanesPerBurg',
    'burgFirstPlaneIndex',
    'burgList',
    'planeList',
]
struct_anon_41._fields_ = [
    ('numBurgVectors', c_int),
    ('numPlanes', c_int),
    ('numPlanesPerBurg', POINTER(c_int)),
    ('burgFirstPlaneIndex', POINTER(c_int)),
    ('burgList', POINTER(real8 * 3)),
    ('planeList', POINTER(real8 * 3)),
]

BurgInfo_t = struct_anon_41 # /home/xzhang11/Planet/Libs/ParaDiSDll/include/Home.h: 156

struct__home.__slots__ = [
    'testHome',
    'myDomain',
    'numDomains',
    'cycle',
    'lastCycle',
    'param',
    'ctrlParamList',
    'dataParamList',
    'nativeNodeQ',
    'ghostNodeQ',
    'freeNodeQ',
    'lastFreeNode',
    'lastGhostNode',
    'nodeBlockQ',
    'nodeKeys',
    'newNodeKeyPtr',
    'newNodeKeyMax',
    'recycledNodeHeap',
    'recycledNodeHeapSize',
    'recycledNodeHeapEnts',
    'cellList',
    'cellCount',
    'nativeCellCount',
    'cellKeys',
    'remoteDomainCount',
    'secondaryRemoteDomainCount',
    'remoteDomains',
    'remoteDomainKeys',
    'decomp',
    'domXmin',
    'domXmax',
    'domYmin',
    'domYmax',
    'domZmin',
    'domZmax',
    'xMaxLevel',
    'yMaxLevel',
    'zMaxLevel',
    'maxPackSiz',
    'inRequests',
    'outRequests',
    'inStatus',
    'outStatus',
    'burgX',
    'burgY',
    'burgZ',
    'nburg',
    'mirrorDomainKeys',
    'currentMirrors',
    'inBuf',
    'outBuf',
    'opList',
    'OpCount',
    'OpListLen',
    'rcvOpList',
    'rcvOpCount',
    'timers',
    'cell2',
    'cell2QentArray',
    'cell2nx',
    'cell2ny',
    'cell2nz',
    'cellCharge',
    'tagMap',
    'tagMapSize',
    'tagMapEnts',
    'glPositions',
    'glWeights',
    'fmLayer',
    'fmNumMPCoeff',
    'fmNumTaylorCoeff',
    'cycleForceCalcCount',
    'rotMatrix',
    'rotMatrixInverse',
    'burgData',
    'ioGroupNum',
    'firstInIOGroup',
    'lastInIOGroup',
    'prevInIOGroup',
    'nextInIOGroup',
    'isFirstInIOGroup',
    'isLastInIOGroup',
    'commLastInIOGroup',
    'surfaceSegCount',
    'surfaceSegList',
]
struct__home._fields_ = [
    ('testHome', c_int),
    ('myDomain', c_int),
    ('numDomains', c_int),
    ('cycle', c_int),
    ('lastCycle', c_int),
    ('param', POINTER(Param_t)),
    ('ctrlParamList', POINTER(ParamList_t)),
    ('dataParamList', POINTER(ParamList_t)),
    ('nativeNodeQ', POINTER(Node_t)),
    ('ghostNodeQ', POINTER(Node_t)),
    ('freeNodeQ', POINTER(Node_t)),
    ('lastFreeNode', POINTER(Node_t)),
    ('lastGhostNode', POINTER(Node_t)),
    ('nodeBlockQ', POINTER(NodeBlock_t)),
    ('nodeKeys', POINTER(POINTER(Node_t))),
    ('newNodeKeyPtr', c_int),
    ('newNodeKeyMax', c_int),
    ('recycledNodeHeap', POINTER(c_int)),
    ('recycledNodeHeapSize', c_int),
    ('recycledNodeHeapEnts', c_int),
    ('cellList', POINTER(c_int)),
    ('cellCount', c_int),
    ('nativeCellCount', c_int),
    ('cellKeys', POINTER(POINTER(Cell_t))),
    ('remoteDomainCount', c_int),
    ('secondaryRemoteDomainCount', c_int),
    ('remoteDomains', POINTER(c_int)),
    ('remoteDomainKeys', POINTER(POINTER(RemoteDomain_t))),
    ('decomp', POINTER(None)),
    ('domXmin', real8),
    ('domXmax', real8),
    ('domYmin', real8),
    ('domYmax', real8),
    ('domZmin', real8),
    ('domZmax', real8),
    ('xMaxLevel', c_int),
    ('yMaxLevel', c_int),
    ('zMaxLevel', c_int),
    ('maxPackSiz', c_int),
    ('inRequests', POINTER(MPI_Request)),
    ('outRequests', POINTER(MPI_Request)),
    ('inStatus', POINTER(MPI_Status)),
    ('outStatus', POINTER(MPI_Status)),
    ('burgX', POINTER(real8)),
    ('burgY', POINTER(real8)),
    ('burgZ', POINTER(real8)),
    ('nburg', c_int),
    ('mirrorDomainKeys', POINTER(POINTER(MirrorDomain_t))),
    ('currentMirrors', c_int),
    ('inBuf', String),
    ('outBuf', String),
    ('opList', POINTER(Operate_t)),
    ('OpCount', c_int),
    ('OpListLen', c_int),
    ('rcvOpList', POINTER(Operate_t)),
    ('rcvOpCount', c_int),
    ('timers', POINTER(Timer_t)),
    ('cell2', POINTER(c_int)),
    ('cell2QentArray', POINTER(C2Qent_t)),
    ('cell2nx', c_int),
    ('cell2ny', c_int),
    ('cell2nz', c_int),
    ('cellCharge', POINTER(real8)),
    ('tagMap', POINTER(TagMap_t)),
    ('tagMapSize', c_int),
    ('tagMapEnts', c_int),
    ('glPositions', POINTER(real8)),
    ('glWeights', POINTER(real8)),
    ('fmLayer', POINTER(FMLayer_t)),
    ('fmNumMPCoeff', c_int),
    ('fmNumTaylorCoeff', c_int),
    ('cycleForceCalcCount', c_int),
    ('rotMatrix', (real8 * 3) * 3),
    ('rotMatrixInverse', (real8 * 3) * 3),
    ('burgData', BurgInfo_t),
    ('ioGroupNum', c_int),
    ('firstInIOGroup', c_int),
    ('lastInIOGroup', c_int),
    ('prevInIOGroup', c_int),
    ('nextInIOGroup', c_int),
    ('isFirstInIOGroup', c_int),
    ('isLastInIOGroup', c_int),
    ('commLastInIOGroup', MPI_Comm),
    ('surfaceSegCount', c_int),
    ('surfaceSegList', POINTER(real8)),
]

struct__param.__slots__ = [
    'nXdoms',
    'nYdoms',
    'nZdoms',
    'nXcells',
    'nYcells',
    'nZcells',
    'iCellNatMin',
    'iCellNatMax',
    'jCellNatMin',
    'jCellNatMax',
    'kCellNatMin',
    'kCellNatMax',
    'xBoundType',
    'yBoundType',
    'zBoundType',
    'xBoundMin',
    'xBoundMax',
    'yBoundMin',
    'yBoundMax',
    'zBoundMin',
    'zBoundMax',
    'minSideX',
    'maxSideX',
    'minSideY',
    'maxSideY',
    'minSideZ',
    'maxSideZ',
    'decompType',
    'DLBfreq',
    'numDLBCycles',
    'cycleStart',
    'maxstep',
    'timeStart',
    'timeNow',
    'timestepIntegrator',
    'deltaTT',
    'realdt',
    'nextDT',
    'maxDT',
    'dtIncrementFact',
    'dtDecrementFact',
    'dtExponent',
    'dtVariableAdjustment',
    'rTol',
    'rmax',
    'minSeg',
    'maxSeg',
    'remeshRule',
    'collisionMethod',
    'remeshAreaMax',
    'remeshAreaMin',
    'splitMultiNodeFreq',
    'fmEnabled',
    'fmNumLayers',
    'fmMPOrder',
    'fmTaylorOrder',
    'fmNumPoints',
    'fmCorrectionTbl',
    'Rijmfile',
    'RijmPBCfile',
    'TempK',
    'loadType',
    'appliedStress',
    'eRate',
    'indxErate',
    'edotdir',
    'cTimeOld',
    'netCyclicStrain',
    'dCyclicStrain',
    'numLoadCycle',
    'eAmp',
    'useLabFrame',
    'labFrameXDir',
    'labFrameYDir',
    'labFrameZDir',
    'mobilityLaw',
    'mobilityType',
    'materialType',
    'vacancyConc',
    'vacancyConcEquilibrium',
    'shearModulus',
    'pois',
    'burgMag',
    'YoungsModulus',
    'rc',
    'Ecore',
    'enforceGlidePlanes',
    'allowFuzzyGlidePlanes',
    'enableCrossSlip',
    'mobilityFunc',
    'MobScrew',
    'MobEdge',
    'MobClimb',
    'MobGlide',
    'MobLine',
    'sessileburgspec',
    'sessilelinespec',
    'includeInertia',
    'massDensity',
    'vAverage',
    'vStDev',
    'dirname',
    'writeBinRestart',
    'doBinRead',
    'numIOGroups',
    'skipIO',
    'armfile',
    'armfilefreq',
    'armfilecounter',
    'armfiledt',
    'armfiletime',
    'fluxfile',
    'fluxfreq',
    'fluxcounter',
    'fluxdt',
    'fluxtime',
    'fragfile',
    'fragfreq',
    'fragcounter',
    'fragdt',
    'fragtime',
    'gnuplot',
    'gnuplotfreq',
    'gnuplotcounter',
    'gnuplotdt',
    'gnuplottime',
    'polefigfile',
    'polefigfreq',
    'polefigcounter',
    'polefigdt',
    'polefigtime',
    'povray',
    'povrayfreq',
    'povraycounter',
    'povraydt',
    'povraytime',
    'atomeye',
    'atomeyefreq',
    'atomeyecounter',
    'atomeyedt',
    'atomeyetime',
    'atomeyesegradius',
    'psfile',
    'psfilefreq',
    'psfiledt',
    'psfiletime',
    'savecn',
    'savecnfreq',
    'savecncounter',
    'savecndt',
    'savecntime',
    'saveprop',
    'savepropfreq',
    'savepropdt',
    'saveproptime',
    'savetimers',
    'savetimersfreq',
    'savetimerscounter',
    'savetimersdt',
    'savetimerstime',
    'savedensityspec',
    'tecplot',
    'tecplotfreq',
    'tecplotcounter',
    'tecplotdt',
    'tecplottime',
    'velfile',
    'velfilefreq',
    'velfilecounter',
    'velfiledt',
    'velfiletime',
    'writeForce',
    'writeForceFreq',
    'writeForceCounter',
    'writeForceDT',
    'writeForceTime',
    'writeVisit',
    'writeVisitFreq',
    'writeVisitCounter',
    'writeVisitSegments',
    'writeVisitSegmentsAsText',
    'writeVisitNodes',
    'writeVisitNodesAsText',
    'writeVisitDT',
    'writeVisitTime',
    'winDefaultsFile',
    'Lx',
    'Ly',
    'Lz',
    'invLx',
    'invLy',
    'invLz',
    'springConst',
    'rann',
    'numBurgGroups',
    'partialDisloDensity',
    'disloDensity',
    'delSegLength',
    'densityChange',
    'TensionFactor',
    'elasticinteraction',
    'delpStrain',
    'delSig',
    'totpStn',
    'delpSpin',
    'totpSpn',
    'totstraintensor',
    'totedgepStrain',
    'totscrewpStrain',
    'dedgepStrain',
    'dscrewpStrain',
    'Ltot',
    'fluxtot',
    'dLtot',
    'dfluxtot',
    'FCC_Ltot',
    'FCC_fluxtot',
    'FCC_dLtot',
    'FCC_dfluxtot',
    'imgstrgrid',
    'node_data_file',
    'dataFileVersion',
    'numFileSegments',
    'nodeCount',
    'dataDecompType',
    'dataDecompGeometry',
    'minCoordinates',
    'maxCoordinates',
    'simVol',
    'burgVolFactor',
    'maxNumThreads',
    'mesh_type',
    'BC_type',
    'dirmax',
    'fem_nx',
    'fem_ny',
    'fem_nz',
    'fem_radius',
    'fem_height',
    'fem_nr',
    'fem_nh',
    'fem_void_radius',
    'fem_cube_edge_length',
    'fem_numelm_arc',
    'fem_base_diagonal_half',
    'fem_grid_base',
    'fem_grid_vertical',
    'fem_ageom_x',
    'fem_ageom_z',
    'fem_delSegLength',
]
struct__param._fields_ = [
    ('nXdoms', c_int),
    ('nYdoms', c_int),
    ('nZdoms', c_int),
    ('nXcells', c_int),
    ('nYcells', c_int),
    ('nZcells', c_int),
    ('iCellNatMin', c_int),
    ('iCellNatMax', c_int),
    ('jCellNatMin', c_int),
    ('jCellNatMax', c_int),
    ('kCellNatMin', c_int),
    ('kCellNatMax', c_int),
    ('xBoundType', BoundType_t),
    ('yBoundType', BoundType_t),
    ('zBoundType', BoundType_t),
    ('xBoundMin', real8),
    ('xBoundMax', real8),
    ('yBoundMin', real8),
    ('yBoundMax', real8),
    ('zBoundMin', real8),
    ('zBoundMax', real8),
    ('minSideX', real8),
    ('maxSideX', real8),
    ('minSideY', real8),
    ('maxSideY', real8),
    ('minSideZ', real8),
    ('maxSideZ', real8),
    ('decompType', c_int),
    ('DLBfreq', c_int),
    ('numDLBCycles', c_int),
    ('cycleStart', c_int),
    ('maxstep', c_int),
    ('timeStart', real8),
    ('timeNow', real8),
    ('timestepIntegrator', c_char * 256),
    ('deltaTT', real8),
    ('realdt', real8),
    ('nextDT', real8),
    ('maxDT', real8),
    ('dtIncrementFact', real8),
    ('dtDecrementFact', real8),
    ('dtExponent', real8),
    ('dtVariableAdjustment', c_int),
    ('rTol', real8),
    ('rmax', real8),
    ('minSeg', real8),
    ('maxSeg', real8),
    ('remeshRule', c_int),
    ('collisionMethod', c_int),
    ('remeshAreaMax', real8),
    ('remeshAreaMin', real8),
    ('splitMultiNodeFreq', c_int),
    ('fmEnabled', c_int),
    ('fmNumLayers', c_int),
    ('fmMPOrder', c_int),
    ('fmTaylorOrder', c_int),
    ('fmNumPoints', c_int),
    ('fmCorrectionTbl', c_char * 256),
    ('Rijmfile', c_char * 256),
    ('RijmPBCfile', c_char * 256),
    ('TempK', real8),
    ('loadType', c_int),
    ('appliedStress', real8 * 6),
    ('eRate', real8),
    ('indxErate', c_int),
    ('edotdir', real8 * 3),
    ('cTimeOld', real8),
    ('netCyclicStrain', real8),
    ('dCyclicStrain', real8),
    ('numLoadCycle', c_int),
    ('eAmp', real8),
    ('useLabFrame', c_int),
    ('labFrameXDir', real8 * 3),
    ('labFrameYDir', real8 * 3),
    ('labFrameZDir', real8 * 3),
    ('mobilityLaw', c_char * 256),
    ('mobilityType', c_int),
    ('materialType', c_int),
    ('vacancyConc', real8),
    ('vacancyConcEquilibrium', real8),
    ('shearModulus', real8),
    ('pois', real8),
    ('burgMag', real8),
    ('YoungsModulus', real8),
    ('rc', real8),
    ('Ecore', real8),
    ('enforceGlidePlanes', c_int),
    ('allowFuzzyGlidePlanes', c_int),
    ('enableCrossSlip', c_int),
    ('mobilityFunc', CFUNCTYPE(UNCHECKED(c_int), POINTER(Home_t), POINTER(Node_t))),
    ('MobScrew', real8),
    ('MobEdge', real8),
    ('MobClimb', real8),
    ('MobGlide', real8),
    ('MobLine', real8),
    ('sessileburgspec', real8 * 30),
    ('sessilelinespec', real8 * 30),
    ('includeInertia', c_int),
    ('massDensity', real8),
    ('vAverage', real8),
    ('vStDev', real8),
    ('dirname', c_char * 256),
    ('writeBinRestart', c_int),
    ('doBinRead', c_int),
    ('numIOGroups', c_int),
    ('skipIO', c_int),
    ('armfile', c_int),
    ('armfilefreq', c_int),
    ('armfilecounter', c_int),
    ('armfiledt', real8),
    ('armfiletime', real8),
    ('fluxfile', c_int),
    ('fluxfreq', c_int),
    ('fluxcounter', c_int),
    ('fluxdt', real8),
    ('fluxtime', real8),
    ('fragfile', c_int),
    ('fragfreq', c_int),
    ('fragcounter', c_int),
    ('fragdt', real8),
    ('fragtime', real8),
    ('gnuplot', c_int),
    ('gnuplotfreq', c_int),
    ('gnuplotcounter', c_int),
    ('gnuplotdt', real8),
    ('gnuplottime', real8),
    ('polefigfile', c_int),
    ('polefigfreq', c_int),
    ('polefigcounter', c_int),
    ('polefigdt', real8),
    ('polefigtime', real8),
    ('povray', c_int),
    ('povrayfreq', c_int),
    ('povraycounter', c_int),
    ('povraydt', real8),
    ('povraytime', real8),
    ('atomeye', c_int),
    ('atomeyefreq', c_int),
    ('atomeyecounter', c_int),
    ('atomeyedt', real8),
    ('atomeyetime', real8),
    ('atomeyesegradius', real8),
    ('psfile', c_int),
    ('psfilefreq', c_int),
    ('psfiledt', real8),
    ('psfiletime', real8),
    ('savecn', c_int),
    ('savecnfreq', c_int),
    ('savecncounter', c_int),
    ('savecndt', real8),
    ('savecntime', real8),
    ('saveprop', c_int),
    ('savepropfreq', c_int),
    ('savepropdt', real8),
    ('saveproptime', real8),
    ('savetimers', c_int),
    ('savetimersfreq', c_int),
    ('savetimerscounter', c_int),
    ('savetimersdt', real8),
    ('savetimerstime', real8),
    ('savedensityspec', c_int * 3),
    ('tecplot', c_int),
    ('tecplotfreq', c_int),
    ('tecplotcounter', c_int),
    ('tecplotdt', real8),
    ('tecplottime', real8),
    ('velfile', c_int),
    ('velfilefreq', c_int),
    ('velfilecounter', c_int),
    ('velfiledt', real8),
    ('velfiletime', real8),
    ('writeForce', c_int),
    ('writeForceFreq', c_int),
    ('writeForceCounter', c_int),
    ('writeForceDT', real8),
    ('writeForceTime', real8),
    ('writeVisit', c_int),
    ('writeVisitFreq', c_int),
    ('writeVisitCounter', c_int),
    ('writeVisitSegments', c_int),
    ('writeVisitSegmentsAsText', c_int),
    ('writeVisitNodes', c_int),
    ('writeVisitNodesAsText', c_int),
    ('writeVisitDT', real8),
    ('writeVisitTime', real8),
    ('winDefaultsFile', c_char * 256),
    ('Lx', real8),
    ('Ly', real8),
    ('Lz', real8),
    ('invLx', real8),
    ('invLy', real8),
    ('invLz', real8),
    ('springConst', real8),
    ('rann', real8),
    ('numBurgGroups', c_int),
    ('partialDisloDensity', POINTER(real8)),
    ('disloDensity', real8),
    ('delSegLength', real8),
    ('densityChange', real8 * 14),
    ('TensionFactor', real8),
    ('elasticinteraction', c_int),
    ('delpStrain', real8 * 6),
    ('delSig', real8 * 6),
    ('totpStn', real8 * 6),
    ('delpSpin', real8 * 6),
    ('totpSpn', real8 * 6),
    ('totstraintensor', real8 * 6),
    ('totedgepStrain', real8 * 6),
    ('totscrewpStrain', real8 * 6),
    ('dedgepStrain', real8 * 6),
    ('dscrewpStrain', real8 * 6),
    ('Ltot', (real8 * 4) * 4),
    ('fluxtot', (real8 * 7) * 4),
    ('dLtot', (real8 * 4) * 4),
    ('dfluxtot', (real8 * 7) * 4),
    ('FCC_Ltot', (real8 * 4) * 6),
    ('FCC_fluxtot', (real8 * 7) * 6),
    ('FCC_dLtot', (real8 * 4) * 6),
    ('FCC_dfluxtot', (real8 * 7) * 6),
    ('imgstrgrid', c_int * 6),
    ('node_data_file', c_char * 256),
    ('dataFileVersion', c_int),
    ('numFileSegments', c_int),
    ('nodeCount', c_int),
    ('dataDecompType', c_int),
    ('dataDecompGeometry', c_int * 3),
    ('minCoordinates', real8 * 3),
    ('maxCoordinates', real8 * 3),
    ('simVol', real8),
    ('burgVolFactor', real8),
    ('maxNumThreads', c_int),
    ('mesh_type', c_int),
    ('BC_type', c_int),
    ('dirmax', c_int),
    ('fem_nx', c_int),
    ('fem_ny', c_int),
    ('fem_nz', c_int),
    ('fem_radius', real8),
    ('fem_height', real8),
    ('fem_nr', c_int),
    ('fem_nh', c_int),
    ('fem_void_radius', real8),
    ('fem_cube_edge_length', real8),
    ('fem_numelm_arc', c_int),
    ('fem_base_diagonal_half', real8),
    ('fem_grid_base', c_int),
    ('fem_grid_vertical', c_int),
    ('fem_ageom_x', real8 * 3),
    ('fem_ageom_z', real8 * 3),
    ('fem_delSegLength', real8),
]

# /home/xzhang11/Planet/Libs/ParaDiSDll/include/Param.h: 538
for _lib in _libs.itervalues():
    if not hasattr(_lib, 'CtrlParamInit'):
        continue
    CtrlParamInit = _lib.CtrlParamInit
    CtrlParamInit.argtypes = [POINTER(Param_t), POINTER(ParamList_t)]
    CtrlParamInit.restype = None
    break

# /home/xzhang11/Planet/Libs/ParaDiSDll/include/Param.h: 539
for _lib in _libs.itervalues():
    if not hasattr(_lib, 'DataParamInit'):
        continue
    DataParamInit = _lib.DataParamInit
    DataParamInit.argtypes = [POINTER(Param_t), POINTER(ParamList_t)]
    DataParamInit.restype = None
    break

# /home/xzhang11/Planet/Libs/ParaDiSDll/include/Param.h: 21
try:
    HDF_DATA_FILE_SUFFIX = '.hdf'
except:
    pass

# /home/xzhang11/Planet/Libs/ParaDiSDll/include/Param.h: 22
try:
    NODEDATA_FILE_SUFFIX = '.data'
except:
    pass

# /home/xzhang11/Planet/Libs/ParaDiSDll/include/Param.h: 23
try:
    NODEDATA_FILE_VERSION = 4
except:
    pass

_param = struct__param # /home/xzhang11/Planet/Libs/ParaDiSDll/include/Param.h: 25

# No inserted files

