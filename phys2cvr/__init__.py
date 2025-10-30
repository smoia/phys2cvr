"""Hopefully importing everything."""

import pkgutil

from . import _version

__version__ = _version.get_versions()['version']

__all__ = []
for loader, module_name, _is_pkg in pkgutil.walk_packages(__path__):
    __all__.append(module_name)

    try:
        _module = loader.find_module(module_name).load_module(module_name)
    except AttributeError:
        _module = loader.find_spec(module_name).loader.load_module(module_name)
    globals()[module_name] = _module
