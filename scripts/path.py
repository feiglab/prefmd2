import os
import sys
import glob
import shutil
import subprocess as sp

class Path:
    def __init__(self, path):
        if isinstance(path, str):
            self._path = os.path.realpath(path)
        else:
            self._path = path._path
        self._isdir = os.path.isdir(self._path)
    def __repr__(self):
        return self._path
    def __len__(self):
        return len(self._path)
    def __getitem__(self, *args, **kwargs):
        return self._path.__getitem__(*args, **kwargs)
    def __cmp__(self, othr, *args):
        return cmp(self._path, othr._path, *args)
    def __lt__(self, othr):
        return self._path < othr._path
    def __contains__(self, key):
        return key in self._path
    def abspath(self):
        return self._path
    def path(self):
        return self._path
    def home(self):
        return Dir(os.path.dirname(self._path))
    def dirname(self):
        return self.home()
    def isdir(self):
        return self._isdir
    def short(self):
        return self.shortpath()
    def shortpath(self):
        if len(self.relpath()) < len(self._path):
            return self.relpath()
        else:
            return self._path
    def relpath(self):
        return os.path.relpath(self._path)
    def realpath(self):
        return Path(os.path.realpath(self._path))
    def endswith(self, key):
        return self._path.endswith(key)
    def status(self, size=0):
        if size < 0:
            return os.path.exists(self._path)
        else:
            if os.path.exists(self._path):
                return os.path.getsize(self._path) > size
            else:
                return False
    def exists(self, size=-1):
        return self.status(size=size)
    def getmtime(self):
        if self.status():
            return os.path.getmtime(self._path)
        else:
            return -1
    def split(self, *args, **kwarg):
        return self._path.split(*args, **kwarg)
    def prefix(self):
        return prefix(self._path)
    def suffix(self, k=-1):
        return suffix(self._path, k=k)
    def fname(self):
        return fname(self._path)
    def name(self):
        return name(self._path)
    def open(self, *args):
        if len(args) == 0:
            mode = 'r'
        else:
            mode = args[0]
        return open(self.path(), mode)
    def remove(self):
        if os.path.exists(self._path):
            os.remove(self._path)
    def copy(self, dest):
        shutil.copyfile(self.path(), "%s"%dest)
    def __eq__(self, other):
        if isinstance(other, Path):
            return self._path == other._path
        else:
            return self._path == other
    def __ne__(self, other):
        return not self.__eq__(other)
    @classmethod
    def glob(cls, pattern):
        fn_s = []
        for fn in glob.glob(pattern):
            if os.path.isdir(fn):
                fn_s.append(Dir(fn))
            else:
                fn_s.append(cls(fn))
        return fn_s

class Dir(Path):
    def __init__(self, path, build=False):
        self._isdir = True
        if not isinstance(path, Path):
            self._path = os.path.realpath(path)
        else:
            self._path = path._path
        if build:
            self.build()
    def build(self):
        if not os.path.exists(self._path):
            os.makedirs(self._path)
    def fn(self, fn):
        return Path("%s/%s"%(self.path(), fn))
    def subdir(self, name, **kwarg):
        return Dir("%s/%s"%(self.path(), name), **kwarg)
    def chdir(self):
        os.chdir(self.path())
    def glob(self, pattern):
        fn_s = []
        for fn in glob.glob('%s/%s'%(self.path(), pattern)):
            if os.path.isdir(fn):
                fn_s.append(Dir(fn))
            else:
                fn_s.append(Path(fn))
        return fn_s
    def remove(self):
        if os.path.exists(self._path):
            sp.call(["rm", "-rf", self._path])
    def realpath(self):
        return Dir(os.path.realpath(self._path))

def exists(path):
    if isinstance(path, Path):
        return path.exists()
    else:
        return os.path.exists(path)

def prefix(in_fn):
    pref = '.'.join(('%s'%in_fn).split(".")[:-1])
    if len(pref) != 0:
        return pref
    else:
        return in_fn

def fname(in_fn):
    return ('%s'%in_fn).split("/")[-1]

def name(in_fn):
    return prefix(fname(in_fn))

def suffix(in_fn, k=-1):
    if k == -1:
        return ('%s'%in_fn).split(".")[-1]
    else:
        return '.'.join(fname(in_fn).split(".")[k:])
