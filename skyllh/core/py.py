# -*- coding: utf-8 -*-

from __future__ import division

import copy
import inspect
import numpy as np
import sys

def typename(t):
    """Returns the name of the given type ``t``.
    """
    return t.__name__

def classname(obj):
    """Returns the name of the class of the class instance ``obj``.
    """
    return typename(type(obj))

def get_byte_size_prefix(size):
    """Determines the biggest size prefix for the given size in bytes such that
    the new size is still greater one.

    Parameters
    ----------
    size : int
        The size in bytes.

    Returns
    -------
    newsize : float
        The new byte size accounting for the byte prefix.
    prefix : str
        The biggest byte size prefix.
    """
    prefix_factor_list = [
        ('', 1), ('K', 1024), ('M', 1024**2), ('G', 1024**3), ('T', 1024**4)]

    prefix_idx = 0
    for (prefix, factor) in prefix_factor_list[1:]:
        if(size / factor < 1):
            break
        prefix_idx += 1

    (prefix, factor) = prefix_factor_list[prefix_idx]
    newsize = size / factor

    return (newsize, prefix)

def getsizeof(objects):
    """Determines the size in bytes the given objects have in memory.
    If an object is a sequence, the size of the elements of the sequence will
    be estimated as well and added to the result. This does not account for the
    multiple occurence of the same object.

    Parameters
    ----------
    objects : sequence of instances of object | instance of object.

    Returns
    -------
    memsize : int
        The memory size in bytes of the given objects.
    """
    if(not issequence(objects)):
        objects = [objects]

    memsize = 0
    for obj in objects:
        if(issequence(obj)):
            memsize += getsizeof(obj)
        else:
            memsize += sys.getsizeof(obj)

    return memsize

def issequence(obj):
    """Checks if the given object ``obj`` is a sequence or not. The definition of
    a sequence in this case is, that the function ``len`` is defined for the
    object.

    .. note::

        A str object is NOT considered as a sequence!

    :return True: If the given object is a sequence.
    :return False: If the given object is a str object or not a sequence.

    """
    if(isinstance(obj, str)):
        return False

    try:
        len(obj)
    except TypeError:
        return False

    return True

def issequenceof(obj, T):
    """Checks if the given object ``obj`` is a sequence with items being
    instances of type ``T``.
    """
    if(not issequence(obj)):
        return False
    for item in obj:
        if(not isinstance(item, T)):
            return False
    return True

def issequenceofsubclass(obj, T):
    """Checks if the given object ``obj`` is a sequence with items being
    sub-classes of class T.
    """
    if(not issequence(obj)):
        return False
    for item in obj:
        if(not issubclass(item, T)):
            return False
    return True

def isproperty(obj, name):
    """Checks if the given attribute is of type property. The attribute must
    exist in ``obj``.

    Parameters
    ----------
    obj : object
        The Python object whose attribute to check for being a property.
    name : str
        The name of the attribute.

    Returns
    -------
    check : bool
        True if the given attribute is of type property, False otherwise.
    """
    return isinstance(type(obj).__dict__[name], property)

def func_has_n_args(func, n):
    """Checks if the given function `func` has `n` arguments.

    Parameters
    ----------
    func : callable
        The function to check.
    n : int
        The number of arguments the function must have.

    Returns
    -------
    check : bool
        True if the given function has `n` arguments. False otherwise.
    """
    check = (len(inspect.getargspec(func)[0]) == n)
    return check

def int_cast(v, errmsg):
    """Casts the given value to an integer value. If the cast is impossible, a
    TypeError is raised with the given error message.
    """
    try:
        v = int(v)
    except:
        raise TypeError(errmsg)
    return v

def float_cast(v, errmsg):
    """Casts the given value to a float. If the cast is impossible, a TypeError
    is raised with the given error message.
    """
    try:
        v = float(v)
    except:
        raise TypeError(errmsg)
    return v

def str_cast(v, errmsg):
    """Casts the given value to a str object.
    If the cast is impossible, a TypeError is raised with the given error
    message.
    """
    try:
        v = str(v)
    except:
        raise TypeError(errmsg)
    return v

def list_of_cast(t, v, errmsg):
    """Casts the given value `v` to a list of items of type `t`.
    If the cast is impossible, a TypeError is raised with the given error
    message.
    """
    if(isinstance(v, t)):
        v = [v]
    if(not issequenceof(v, t)):
        raise TypeError(errmsg)
    v = list(v)
    return v

def get_smallest_numpy_int_type(values):
    """Returns the smallest numpy integer type that can represent the given
    integer values.

    Parameters
    ----------
    values : int | sequence of int
        The integer value(s) that need to be representable by the returned
        integer type.

    Returns
    -------
    inttype : numpy integer type
        The smallest numpy integer type that can represent the given values.
    """
    values = np.atleast_1d(values)

    vmin = np.min(values)
    vmax = np.max(values)

    if(vmin < 0):
        types = [np.int8, np.int16, np.int32, np.int64]
    else:
        types = [np.uint8, np.uint16, np.uint32, np.uint64]

    for inttype in types:
        ii = np.iinfo(inttype)
        if(vmin >= ii.min and vmax <= ii.max):
            return inttype

    raise ValueError("No integer type spans [%d, %d]!"%(vmin, vmax))

def _get_func_range():
    """Returns a lazy iterable `range` function to be consistent with Python 3.

    Returns
    -------
    func : function
        Lazy iterable `range` function.
    """
    try:
        func = xrange
    except NameError:
        func = range

    return func

# Overwrite of built-in `range` function to be consistent with Python 3.
range = _get_func_range()


class ObjectCollection(object):
    """This class provides a collection of objects of a specific type. Objects
    can be added to the collection via the ``add`` method or can be removed
    from the collection via the ``pop`` method. The objects of another object
    collection can be added to this object collection via the ``add`` method as
    well.
    """
    def __init__(self, obj_type, obj_list=None):
        """Constructor of the ObjectCollection class. Must be called by the
        derived class.

        Parameters
        ----------
        obj_type : type
            The type of the objects, which can be added to the collection.
        obj_list : list of obj_t | None
            The list of objects of type ``obj_t`` with which this collection
            should get initialized with.
        """
        if(not issubclass(obj_type, object)):
            raise TypeError('The obj_t argument must be a subclass of object!')
        self._obj_type = obj_type
        self._objects = []

        # Add given list of objects.
        if(obj_list is not None):
            if(not issequence(obj_list)):
                raise TypeError('The obj_list argument must be a sequence!')
            for obj in obj_list:
                self.add(obj)

    @property
    def obj_type(self):
        """(read-only) The object type.
        """
        return self._obj_type

    @property
    def objects(self):
        """(read-only) The list of objects of this object collection.
        All objects are of the same type as specified through the ``obj_type``
        property.
        """
        return self._objects

    def __len__(self):
        """Returns the number of objects being in this object collection.
        """
        return len(self._objects)

    def __getitem__(self, key):
        return self._objects[key]

    def __iter__(self):
        return iter(self._objects)

    def __add__(self, other):
        """Implementation to support the operation ``oc = self + other``, where
        ``self`` is this ObjectCollection object and ``other`` something useful
        else. This creates a copy ``oc`` of ``self`` and adds ``other``
        to ``oc``.

        Parameters
        ----------
        other : obj_type | ObjectCollection of obj_type

        Returns
        -------
        oc : ObjectCollection
            The new ObjectCollection object with object from self and other.
        """
        oc = self.copy()
        oc.add(other)
        return oc

    def __str__(self):
        """Pretty string representation of this object collection.
        """
        return classname(self)+ ': ' + str(self._objects)

    def copy(self):
        """Creates a copy of this ObjectCollection. The objects of the
        collection are not copied!
        """
        oc = ObjectCollection(self._obj_type)
        oc._objects = copy.copy(self._objects)
        return oc

    def add(self, obj):
        """Adds the given object to the collection.

        Parameters
        ----------
        obj : obj_type | ObjectCollection of obj_type
            An instance of ``obj_type`` that should be added to the collection.
            If given an ObjectCollection for objects of type obj_type, it will
            add all objects of the given collection to this collection.

        Returns
        -------
        self : ObjectCollection
            The instance of this ObjectCollection, in order to be able to chain
            several ``add`` calls.
        """
        if(isinstance(obj, ObjectCollection)):
            if(typename(obj.obj_type) != typename(self._obj_type)):
                raise TypeError('Cannot add objects from ObjectCollection for objects of type "%s" to this ObjectCollection for objects of type "%s"!'%(typename(obj.obj_type), typename(self._obj_type)))
            self._objects.extend(obj.objects)
            return self

        if(not isinstance(obj, self._obj_type)):
            raise TypeError('The object of type "%s" cannot be added to the object collection for objects of type "%s"!'%(classname(obj), typename(self._obj_type)))

        self._objects.append(obj)
        return self
    __iadd__ = add

    def index(self, obj):
        return self._objects.index(obj)

    def pop(self, index=None):
        """Removes and returns object at given index (default last).
        Raises IndexError if the collection is empty or index is out of range.

        Parameters
        ----------
        index : int | None
            The index of the object to remove. If set to None, the index of the
            last object is used.

        Returns
        -------
        obj : obj_type
            The removed object.
        """
        if(index is None):
            index = len(self._objects)-1
        obj = self._objects.pop(index)
        return obj
