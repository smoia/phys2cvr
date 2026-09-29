.. _api:

.. currentmodule:: phys2cvr


phys2cvr package
=============

:mod:`phys2cvr.workflows` - Primary workflows
------------------------------------------

.. currentmodule:: phys2cvr.workflows

.. autosummary::
   :template: function.rst
   :toctree: generated/

   phys2cvr

:mod:`phys2cvr.regressors` - Create regressors
-------------------------------------------

.. currentmodule:: phys2cvr.regressors

.. autosummary::
   :template: function.rst
   :toctree: generated/

   create_legendre
   compute_bulk_shift
   create_fine_shift_regressors
   create_physio_regressor

:mod:`phys2cvr.stats` - Statistics
-------------------------------

.. currentmodule:: phys2cvr.stats

.. autosummary::
   :template: function.rst
   :toctree: generated/

   R2MODEL
   x_corr
   ols
   regression

:mod:`phys2cvr.signal` - Signal processing related
-----------------------------------------------

.. currentmodule:: phys2cvr.signal

.. autosummary::
   :template: function.rst
   :toctree: generated/

   spc
   create_hrf
   filter_signal
   compute_petco2hrf
   resample_signal_samples
   resample_signal_freqs

:mod:`phys2cvr.io` - I/O functions: load
-------------------------------------

.. currentmodule:: phys2cvr.io

.. autosummary::
   :template: function.rst
   :toctree: generated/

   load_nifti_get_mask
   load_txt
   load_mat
   load_xls
   load_array
   load_physio

:mod:`phys2cvr.io` - I/O functions: export
---------------------------------------

.. currentmodule:: phys2cvr.io

.. autosummary::
   :template: function.rst
   :toctree: generated/

   export_regressor
   export_nifti

:mod:`phys2cvr.io` - I/O functions: supported extensions
-----------------------------------------------------

.. currentmodule:: phys2cvr.io

.. autosummary::
   :template: attribute.rst
   :toctree: generated/

   EXT_1D
   EXT_MAT
   EXT_ARRAY
   EXT_PHYS
   EXT_NIFTI
   EXT_GIFTI
   EXT_NIMG
   EXT_ALL

:mod:`phys2cvr.viz` - Visualisations
---------------------------------

.. currentmodule:: phys2cvr.viz

.. autosummary::
   :template: function.rst
   :toctree: generated/

   plot_timeseries
   plot_xcorr

:mod:`phys2cvr.utils` - Utility functions
--------------------------------------

.. currentmodule:: phys2cvr.utils

.. autosummary::
   :template: function.rst
   :toctree: generated/

   if_declared_force_type
   check_ext
   check_nifti_dim
   check_array_dim
   save_bash_call
