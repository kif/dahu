Getting started with Dahu
=========================

Register Dahu with an existing Tango database
---------------------------------------------

Assuming you have a tango database running on port 10000,
you can register a server instance (called "dahu" in this example)

.. code-block:: bash

    TANGO_HOST=localhost:10000 dahu-register --instance dahu

Note: you can register more than one Dahu tango device with the same server.
Type `dahu-register --help` for more information.

Dahu server with database
-------------------------

To start a Dahu server that has been registered with a
tango database running on port 10000:

.. code-block:: bash

    TANGO_HOST=localhost:10000 dahu-server dahu

To use a tango device managed by this server:

.. code-block:: python

    from PyTango import DeviceProxy
    p = DeviceProxy("tango://localhost:10000/id00/dahu/1")


Dahu server without database
----------------------------

To start a Dahu server that has not been registered with a
tango database:

.. code-block:: bash

    dahu-server dahu -ORBendPoint giop:tcp::10001 -nodb -dlist id00/dahu/1 -v4

To use a tango device managed by this server:

.. code-block:: python

    from PyTango import DeviceProxy
    proxy = DeviceProxy("tango://localhost:10001/id00/dahu/1#dbase=no")
