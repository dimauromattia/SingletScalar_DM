
.. DO NOT EDIT.
.. THIS FILE WAS AUTOMATICALLY GENERATED BY SPHINX-GALLERY.
.. TO MAKE CHANGES, EDIT THE SOURCE PYTHON FILE:
.. "examples_gallery/plot_indirect_detection.py"
.. LINE NUMBERS ARE GIVEN BELOW.

.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        :ref:`Go to the end <sphx_glr_download_examples_gallery_plot_indirect_detection.py>`
        to download the full example code

.. rst-class:: sphx-glr-example-title

.. _sphx_glr_examples_gallery_plot_indirect_detection.py:


Indirect detection
==================

This example studies the likelihoods for different indirect detection searches.

.. GENERATED FROM PYTHON SOURCE LINES 7-15

.. code-block:: default


    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.interpolate import interp2d
    from singletscalar_dm import *

    table_int = np.loadtxt(import_data_file('SHP_sigmav_table.dat'))








.. GENERATED FROM PYTHON SOURCE LINES 16-17

Here we define the vectors for the annihilation cross section, mass and lambda_hs

.. GENERATED FROM PYTHON SOURCE LINES 17-22

.. code-block:: default


    sigmav_vec = np.power(10.,np.arange(-29.,-20.96,0.04))
    mass_vec = massz_vec[19:224]
    lambdap_vec = np.power(10.,np.arange(-6.,0.,0.04))








.. GENERATED FROM PYTHON SOURCE LINES 23-24

Here we create the likelihood profile for dwarf.

.. GENERATED FROM PYTHON SOURCE LINES 24-39

.. code-block:: default


    table_dwarf = np.loadtxt(import_data_file('LogLike_stacked_paper_SHP.dat'))
    LogLike_table_dwarf = np.zeros(shape=(len(sigmav_vec),len(mass_vec)))
    for t in range(len(mass_vec)):
        for u in range(len(sigmav_vec)):
            LogLike_table_dwarf[u,t] = -2.*table_dwarf[t*len(sigmav_vec)+u,3]
    funcint_dwarf = interp2d(mass_vec,sigmav_vec,LogLike_table_dwarf)

    LogLikel_table_dwarf = np.zeros(shape=(len(lambdap_vec),len(mass_vec)))
    for t in range(len(mass_vec)):
        for u in range(len(lambdap_vec)):
            sigmav_val = lambda2sigmav(mass_vec[t],lambdap_vec[u],table_int)
            LogLikel_table_dwarf[u,t] = funcint_dwarf(mass_vec[t],sigmav_val)









.. GENERATED FROM PYTHON SOURCE LINES 40-41

Here we create the likelihood profile for dwarf.

.. GENERATED FROM PYTHON SOURCE LINES 41-56

.. code-block:: default


    table_pbar = np.loadtxt(import_data_file('LogLike_Manconi2021_pbar_paper.dat'))
    LogLike_table_pbar = np.zeros(shape=(len(sigmav_vec),len(mass_vec)))
    for t in range(len(mass_vec)):
        for u in range(len(sigmav_vec)):
            LogLike_table_pbar[u,t] = table_pbar[t*len(sigmav_vec)+u,2]
    funcint_pbar = interp2d(mass_vec,sigmav_vec,LogLike_table_pbar)

    LogLikel_table_pbar = np.zeros(shape=(len(lambdap_vec),len(mass_vec)))
    for t in range(len(mass_vec)):
        for u in range(len(lambdap_vec)):
            sigmav_val = lambda2sigmav(mass_vec[t],lambdap_vec[u],table_int)
            LogLikel_table_pbar[u,t] = funcint_pbar(mass_vec[t],sigmav_val)









.. GENERATED FROM PYTHON SOURCE LINES 57-58

Here we consider the :math:`\chi^2` for the GCE, relative to the MED model.

.. GENERATED FROM PYTHON SOURCE LINES 58-73

.. code-block:: default


    table_gce = np.loadtxt(import_data_file('Chi_table_Cholis_GCE_MED_paper.dat'))
    LogLike_table_gce = np.zeros(shape=(len(sigmav_vec),len(mass_vec)))
    for t in range(len(mass_vec)):
        for u in range(len(sigmav_vec)):
            LogLike_table_gce[u,t] = -table_gce[t*len(sigmav_vec)+u,2]
    funcint_gce = interp2d(mass_vec,sigmav_vec,LogLike_table_gce)

    LogLikel_table_gce = np.zeros(shape=(len(lambdap_vec),len(mass_vec)))
    for t in range(len(mass_vec)):
        for u in range(len(lambdap_vec)):
            sigmav_val = lambda2sigmav(mass_vec[t],lambdap_vec[u],table_int)
            LogLikel_table_gce[u,t] = funcint_gce(mass_vec[t],sigmav_val)









.. GENERATED FROM PYTHON SOURCE LINES 74-75

Further, we combine the likelihoods.

.. GENERATED FROM PYTHON SOURCE LINES 75-79

.. code-block:: default


    LogLike_table_combined = LogLike_table_gce+LogLike_table_pbar+LogLike_table_dwarf
    LogLikel_table_combined = LogLikel_table_gce+LogLikel_table_pbar+LogLikel_table_dwarf








.. GENERATED FROM PYTHON SOURCE LINES 80-81

And we plot the results.

.. GENERATED FROM PYTHON SOURCE LINES 81-103

.. code-block:: default


    fig, ax = plt.subplots(figsize=(8,6))

    dlin = ( 80+2 )/30.
    scale_vec = np.arange( -80,2, dlin )
    scale_cb = np.arange( -80,2, dlin*10.)

    cf = ax.contourf(mass_vec, lambdap_vec, LogLikel_table_combined-LogLikel_table_combined.max(), 30, levels=list(scale_vec), cmap='seismic')
    cbar = fig.colorbar(cf, ax = ax)
    cbar.ax.set_ylabel(r'$-\Delta\chi^2$ Combined AstroP', fontsize="large")

    plt.ylabel(r'$\lambda_{HS}$', fontsize=18)
    plt.xlabel(r'$m_{S}$ [GeV]', fontsize=18)
    plt.axis([35.,90,1e-4,1.0])
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.grid(True)
    plt.yscale('log')
    plt.xscale('linear') 
    fig.tight_layout(pad=0.5)
    plt.show()




.. image-sg:: /examples_gallery/images/sphx_glr_plot_indirect_detection_001.png
   :alt: plot indirect detection
   :srcset: /examples_gallery/images/sphx_glr_plot_indirect_detection_001.png
   :class: sphx-glr-single-img





.. GENERATED FROM PYTHON SOURCE LINES 105-129

.. code-block:: default


    fig, ax = plt.subplots(figsize=(8,6))

    dlin = ( 80+2 )/30.
    scale_vec = np.arange( -80,2, dlin )
    scale_cb = np.arange( -80,2, dlin*10.)

    cf = ax.contourf(mass_vec-62.50, lambdap_vec, LogLikel_table_combined-LogLikel_table_combined.max(), 30, levels=list(scale_vec), cmap='seismic')
    cbar = fig.colorbar(cf, ax = ax)
    cbar.ax.set_ylabel(r'$-\Delta\chi^2$ Combined AstroP', fontsize="large")

    plt.text(0.11,2e-5,'MED', fontsize=18)
    plt.ylabel(r'$\lambda_{HS}$', fontsize=18)
    plt.xlabel(r'$m_{S}-m_h/2$ [GeV]', fontsize=18)
    plt.axis([-0.2,0.2,1e-5,0.005])
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.grid(True)
    plt.yscale('log')
    plt.xscale('linear') 
    fig.tight_layout(pad=0.5)
    plt.show()





.. image-sg:: /examples_gallery/images/sphx_glr_plot_indirect_detection_002.png
   :alt: plot indirect detection
   :srcset: /examples_gallery/images/sphx_glr_plot_indirect_detection_002.png
   :class: sphx-glr-single-img





.. GENERATED FROM PYTHON SOURCE LINES 130-131

Here we create the likelihood profile for antiprotons from Balan et al. 2023.

.. GENERATED FROM PYTHON SOURCE LINES 131-145

.. code-block:: default


    table_pbar = np.loadtxt(import_data_file('LogLike_sigmav_Manconi2023_model1_corr_pbar_paper.dat'))
    LogLike_table_pbar = np.zeros(shape=(len(sigmav_vec),len(mass_vec)))
    for t in range(len(mass_vec)):
        for u in range(len(sigmav_vec)):
            LogLike_table_pbar[u,t] = table_pbar[t*len(sigmav_vec)+u,2]
    funcint_pbar = interp2d(mass_vec,sigmav_vec,LogLike_table_pbar)

    LogLikel_table_pbar = np.zeros(shape=(len(lambdap_vec),len(mass_vec)))
    for t in range(len(mass_vec)):
        for u in range(len(lambdap_vec)):
            sigmav_val = lambda2sigmav(mass_vec[t],lambdap_vec[u],table_int)
            LogLikel_table_pbar[u,t] = funcint_pbar(mass_vec[t],sigmav_val)








.. GENERATED FROM PYTHON SOURCE LINES 146-147

We then combine the likelihoods.

.. GENERATED FROM PYTHON SOURCE LINES 147-151

.. code-block:: default


    LogLike_table_combined = LogLike_table_gce+LogLike_table_pbar+LogLike_table_dwarf
    LogLikel_table_combined = LogLikel_table_gce+LogLikel_table_pbar+LogLikel_table_dwarf








.. GENERATED FROM PYTHON SOURCE LINES 152-153

And we plot the results.

.. GENERATED FROM PYTHON SOURCE LINES 153-174

.. code-block:: default


    fig, ax = plt.subplots(figsize=(8,6))

    dlin = ( 80+2 )/30.
    scale_vec = np.arange( -80,2, dlin )
    scale_cb = np.arange( -80,2, dlin*10.)

    cf = ax.contourf(mass_vec, lambdap_vec, LogLikel_table_combined-LogLikel_table_combined.max(), 30, levels=list(scale_vec), cmap='seismic')
    cbar = fig.colorbar(cf, ax = ax)
    cbar.ax.set_ylabel(r'$-\Delta\chi^2$ Combined AstroP', fontsize="large")

    plt.ylabel(r'$\lambda_{HS}$', fontsize=18)
    plt.xlabel(r'$m_{S}$ [GeV]', fontsize=18)
    plt.axis([35.,90,1e-4,1.0])
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.grid(True)
    plt.yscale('log')
    plt.xscale('linear') 
    fig.tight_layout(pad=0.5)
    plt.show()



.. image-sg:: /examples_gallery/images/sphx_glr_plot_indirect_detection_003.png
   :alt: plot indirect detection
   :srcset: /examples_gallery/images/sphx_glr_plot_indirect_detection_003.png
   :class: sphx-glr-single-img






.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 0 minutes  37.051 seconds)


.. _sphx_glr_download_examples_gallery_plot_indirect_detection.py:

.. only:: html

  .. container:: sphx-glr-footer sphx-glr-footer-example




    .. container:: sphx-glr-download sphx-glr-download-python

      :download:`Download Python source code: plot_indirect_detection.py <plot_indirect_detection.py>`

    .. container:: sphx-glr-download sphx-glr-download-jupyter

      :download:`Download Jupyter notebook: plot_indirect_detection.ipynb <plot_indirect_detection.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
