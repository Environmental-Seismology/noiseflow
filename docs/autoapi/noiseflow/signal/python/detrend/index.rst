:py:mod:`noiseflow.signal.python.detrend`
=========================================

.. py:module:: noiseflow.signal.python.detrend


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   noiseflow.signal.python.detrend.detrend
   noiseflow.signal.python.detrend.detrend_py
   noiseflow.signal.python.detrend.spline



.. py:function:: detrend(data, type='linear', order=2, dspline=1000)

   Detrend data.
   :param data: The data to detrend.
   :type data: :class:`numpy.ndarray`
   :param type: The type of detrending to perform. Can be one of ``'demean'``,
       ``'linear'``, or ``'spline'``.
   :type type: str
   :param overwrite_data: If True the data will be detrended in place.
   :type overwrite_data: bool
   :param dtype: The data type to use for the detrended data. Defaults to
       ``np.float32``.
   :type dtype: :class:`numpy.dtype`
   :param order: The order of the spline to fit. Defaults to 2.
   :type order: int
   :param dspline: The distance between spline knots. Defaults to 1000.
   :type dspline: int

   :returns: The detrended data.
   :rtype: :class:`numpy.ndarray`

   .. rubric:: Example

   >>> import numpy as np
   >>> from noiseflow.signal import detrend
   >>> data = np.array([1.1,2.7,3,3,3,6,7,8,9,10], overwrite_data=True)
   >>> detrend(data, type='spline')



.. py:function:: detrend_py(data, type='linear', flag=False, flag_gap=None)


.. py:function:: spline(data, order=2, dspline=1000)


