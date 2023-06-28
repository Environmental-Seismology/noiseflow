:py:mod:`noiseflow.cc.wrapper`
==============================

.. py:module:: noiseflow.cc.wrapper


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   noiseflow.cc.wrapper.corr
   noiseflow.cc.wrapper.rfft
   noiseflow.cc.wrapper.stack



Attributes
~~~~~~~~~~

.. autoapisummary::

   noiseflow.cc.wrapper.NOISEFLOW_USE_CPP
   noiseflow.cc.wrapper.config
   noiseflow.cc.wrapper.config_path
   noiseflow.cc.wrapper.env


.. py:function:: corr(rfft_data, dt, corr_method, corr_pair, maxlag, smoothspect_N=10, flag=False, flag_gap=None, threads=1, jobs=1, py=False)


.. py:function:: rfft(raw_data, dt, cc_len, cc_step, time_norm, clip_std, smooth_N, freq_norm, freqmin, freqmax, whiten_npad, smoothspect_N, flag=False, flag_gap=None, threads=1, jobs=1, py=False)


.. py:function:: stack(corr_data, dt, stack_method, par=None, stack_all=True, stack_len=0, stack_step=0, pick=False, median_high=10, median_low=0.1, flag=False, flag_gap=None, threads=1, jobs=1, py=False)


.. py:data:: NOISEFLOW_USE_CPP
   :value: True

   

.. py:data:: config

   

.. py:data:: config_path

   

.. py:data:: env

   

