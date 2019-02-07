MULTIPLY prior-engine
=======================

|buildstatus| |docstatus|

.. htmlonly::

    `View on GitHub <https://github.com/multiply-org/prior-engine>`_


Scope of MULTIPLY
------------------

The MULTIPLY project will “develop a new platform for joint and consistent retrieval of Copernicus
SENTINEL data and beyond”.

This documentation covers the prior engine for the MULTIPLY main platform.
This module provides *a priori* information to the `Inference Engine <https://github.com/multiply-org/KaFKA-InferenceEngine>`_ to support land surface parameter retrieval.

The `prior engine specific documentation <https://multiply-prior-engine.readthedocs.io/en/latest/>`_ is hosted on ReadTheDocs. It is part of the `MULTIPLY core documentation <http://multiply.readthedocs.io/en/latest/>`_.
Please find the latest pdf version of this documentation `here <https://readthedocs.org/projects/multiply-prior-engine/downloads/pdf/latest/>`_.

First Steps
-----------

Getting Started
^^^^^^^^^^^^^^^
Please find instructions on how to download and install the prior engine in the :ref:`Installation` section.

.. note::
   TBD: Getting started with python, bayes theorem, ..


Testing and Contribution
^^^^^^^^^^^^^^^^^^^^^^^^^

You are welcome to test and contribute to the MULTIPLY Prior Engine.

Please find corresponding guidelines and further information on how to do so in the :ref:`Contribution` section and on the `project GitHub page <https://github.com/multiply-org/prior-engine>`_.



Content
---------

.. toctree::
   :maxdepth: 2

   Introduction
   Installation
   Usage
   Processing



Developer Documentation
------------------------

.. toctree::
   :maxdepth: 2

   Changelog
   Contribution
   Testing
   multiply_prior_engine
   License


Indices and tables
-------------------

* :ref:`modindex`
* :ref:`genindex`
* :ref:`search`




.. |logo| image:: https://raw.githubusercontent.com/multiply-org/multiply-core/master/doc/source/_static/logo/Multiply_multicolour.png
   :width: 10%

.. |buildstatus| image:: https://travis-ci.org/multiply-org/prior-engine.svg?branch=master
    :target: https://travis-ci.org/multiply-org/prior-engine

.. |docstatus| image:: https://readthedocs.org/projects/multiply-prior-engine/badge/?version=latest
    :target: https://multiply-prior-engine.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status
