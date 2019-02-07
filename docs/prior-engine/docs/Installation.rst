.. _installation:

Installation
============

Download
--------
If not already done so, the first step is to clone the latest code and change directory:

.. code-block:: bash
    :linenos:

    git clone https://github.com/multiply-org/prior-engine.git
    cd prior-engine

.. note::
    The MULTIPLY platform has been developed against Python 3.6.
    It cannot be guaranteed to work with previous Python versions.

Installation procedure
--------------------------
The MULTIPLY prior engine can be run from sources directly.
To install the MULTIPLY prior engine into an existing Python environment just for the current user, use

.. code-block:: bash

    python setup.py install --user

To install the MULTIPLY Core for development and for the current user, use

.. code-block:: bash

    python setup.py develop --user

Using Conda
^^^^^^^^^^^
.. note::
   TBD


Module requirements
-------------------
from `requirements.txt`:

.. literalinclude:: ../requirements.txt
