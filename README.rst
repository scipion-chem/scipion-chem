================================
Atom_struct_utils scipion plugin
================================

Atomic model utilties on the web

===================
Install this plugin
===================

You will need to use `3.0.0 <https://github.com/I2PC/scipion/releases/tag/v3.0>`_ version of Scipion to run these protocols. To install the plugin, you have two options:

- **Stable version**  

.. code-block:: 

      scipion installp -p scipion-em-atom_struct_utils_web
      
OR

  - through the plugin manager GUI by launching Scipion and following **Configuration** >> **Plugins**
      
- **Developer's version** 

1. Download repository: 

.. code-block::

            git clone https://github.com/scipion-em/scipion-em-atomstructutilsWeb.git

2. Install:

.. code-block::

            scipion installp -p path_to_scipion-em-atomstructutilsWeb --devel

- **Binary files** 

Atom_struct_utils plugin is a pure Python module, no binary files are required. 

- **Tests**

To check the installation, simply run the following Scipion test:

===============
Buildbot status
===============

Status devel version: 

.. image:: http://scipion-test.cnb.csic.es:9980/badges/atomstructutilsWeb_devel.svg

Status production version: 

.. image:: http://scipion-test.cnb.csic.es:9980/badges/atomstructutilsWeb_prod.svg
