Usage
=====

Python Package
------------------

MULTIPLY prior engine is available as Python Package.
To import it into your python application, use

.. code-block:: python

   import multiply_prior_engine


Command Line Interface
------------------------

There is a Command Line Interface (CLI) integrated to allow for the following actions:

- add user defined prior data,
- import user defined prior data,
- remove/un-select prior data from configuration,
- show configuration.

The CLI's help can be accessed via `-h` flag:

.. code-block:: bash

   user_prior -h


The help and description of the above mentioned sub-commands can be accessed via, e.g.:

.. code-block:: bash

   user_prior add -h


Current limitations in the user defined priors
----------------------------------------------------

So far, priors can only be added as point data for specific variables. User defined prior data has to be passed to the engine in the form of comma separated values (csv) with dates in the first column and the parameter values in the second column.
There is the requirement for a header line specifying the variable name (lai, sm, ...) and geolocation (latitude, longitude) of the data.
E.g.:

.. code-block:: csv

    lai, 10.5564, 48.3124
    2017-06-01, 1.01
    2017-06-02, 1.01
    2017-06-03, 1.2
    2017-06-04, 1.25
    2017-06-05, 1.4
    ....
