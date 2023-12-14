================================
CHEM scipion plugin
================================

**Documentation under development, sorry for the inconvenience**

Base Scipion plugin defining objects and protocols for CHEMoinformatics and virtual drug screening

Full documentation can be found in the `Scipion Chem official documentation page <https://scipion-chem.github.io/docs/index.html>`_.

===================
Install this plugin
===================

You will need to use `3.0.0 <https://github.com/I2PC/scipion/releases/tag/v3.0>`_ version of Scipion
to run these protocols. To install the plugin, you have two options:

- **Stable version (Currently not available)**

.. code-block:: 

      scipion3 installp -p scipion-chem
      
OR

  - through the plugin manager GUI by launching Scipion and following **Configuration** >> **Plugins**
      
- **Developer's version** 

1. **Download repository**:

.. code-block::

            git clone https://github.com/scipion-chem/scipion-chem.git

2. **Switch to the desired branch** (master or devel):

Scipion-chem is constantly under development.
If you want a relatively older an more stable version, use master branch (default).
If you want the latest changes and developments, user devel branch.

.. code-block::

            cd scipion-chem
            git checkout devel

3. **Install**:
The following comand will launch the installation of the plugin, together with some third-party programs

.. code-block::

            scipion3 installp -p path_to_scipion-chem --devel

- **External software**

External software currently installed by scipion-chem:

- **OpenBabel** and **RDKit**: the main small molecule handlers and converters
- **MGLTools**: additional utils for small molecules, docking, ... (includes AutoDockTools)
- **JChemPaint**: Java program to manually draw small molecules.
- **PyMol**: main viewer of Scipion-Chem for small molecules and structures
- **VMD**: secondary viewer of Scipion-Chem for structures and Molecular Dynamics
- **AliView**: main viewer for sequences
- **PLIP**: specialized viewer for docking interactions in PyMol

===============
Buildbot status
===============

Status devel version: 

.. image:: http://scipion-test.cnb.csic.es:9980/badges/chem_devel.svg

Status production version: 

.. image:: http://scipion-test.cnb.csic.es:9980/badges/chem_prod.svg
