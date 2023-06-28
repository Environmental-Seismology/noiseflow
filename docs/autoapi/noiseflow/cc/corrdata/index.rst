:py:mod:`noiseflow.cc.corrdata`
===============================

.. py:module:: noiseflow.cc.corrdata


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   noiseflow.cc.corrdata.CorrData_Class




.. py:class:: CorrData_Class(corr_data, dt, corr_method, corr_pair, maxlag, smoothspect_N, flag, flag_gap, threads, jobs, py)

   Bases: :py:obj:`object`

   .. py:method:: plot(pair_indx=0, t_min=UTCDateTime('1970-01-01T00:00:00.0'), cc_len=None, cc_step=None, win_start=None, win_end=None, lag_start=None, lag_end=None, amp_normalize=True, amp_scale=1, filter=False, f1=None, f2=None, corners=4, zerophase=True, win_interval=None, mode='waveform', cmap='seismic', linewidth=0.8, yticklabel_num=5, figsize=(10, 6), save=False, save_path=None, dpi=100)


   .. py:method:: save(save_path, format='npz', compression=False, h5_compression_format='gzip', h5_compression_opts=3)



