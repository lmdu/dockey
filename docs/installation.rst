Installation
============

Prior to installing Dockey, you should first install a docking engine, any of `AutoDock4 <https://autodock.scripps.edu/>`_, `AutoDock Vina <https://vina.scripps.edu/>`_ or `QuickVina-W <https://qvina.github.io/>`_. We have tested Dockey on these systems:

.. list-table::
   :widths: 15 25 60
   :header-rows: 1

   * - 
     - System
     - Versions
   * - |windows|
     - Windows
     - 10, 11
   * - |linux|
     - Ubuntu
     - 20.04
   * - |macos|
     - MacOS
     - Big Sur 11.7

Windows
-------

**AutoDock4:**

Go to `https://autodock.scripps.edu/download-autodock4/ <https://autodock.scripps.edu/download-autodock4/>`_ page, click on ``AutoDock4 Windows`` to download ``autodocksuite-4.2.6.i86Windows.exe`` installer, and then double click installer to install AutoDock4.

**AutoDock Vina:**

Go to `https://github.com/ccsb-scripps/AutoDock-Vina/releases <https://github.com/ccsb-scripps/AutoDock-Vina/releases>`_ page, click on ``vina_1.2.3_windows_x86_64.exe`` to download it, and put it any where.

**QuickVina-W:**

The QuickVina homepage only provide linux executable. For Windows, you can install QuickVina-W by `anaconda <https://anaconda.org/conda-forge/qvina>`_.

**Dockey:**

Go to `https://github.com/lmdu/dockey/releases <https://github.com/lmdu/dockey/releases>`_ page, click on ``Dockey-version-win64.exe`` to download it. Then double click the downloaded installer to install the program following the on-screen instructions.

Linux
-----

**AutoDock4:**

Go to `https://autodock.scripps.edu/download-autodock4/ <https://autodock.scripps.edu/download-autodock4/>`_ page, click on ``AutoDock v4.2.6 (Linux 64bit)`` to download ``autodocksuite-4.2.6-x86_64Linux2.tar`` package, Then uncompress the package:

.. code:: shell

	tar xvf autodocksuite-4.2.6-x86_64Linux2.tar

After uncompress, you will get a folder named ``x86_64Linux2`` which contains autogrid4 and autodock4 executables. You can move these two executables to any where, for example, move to /usr/local/bin:

.. code:: shell

	cd x86_64Linux2
	sudo mv autogrid4 autodock4 /usr/local/bin

**AutoDock Vina:**

Go to `https://github.com/ccsb-scripps/AutoDock-Vina/releases <https://github.com/ccsb-scripps/AutoDock-Vina/releases>`_ page, click on ``vina_1.2.3_linux_x86_64`` to download it, and put it any where.

**QuickVina-W:**

Go to `https://qvina.github.io/ <https://qvina.github.io/>`_ page, and then click ``Download binary QuickVina-W`` button on the right side of page to download ``qvina-w`` executable, and put it any where.

**Dockey:**

.. note::

	Currently, we only tested Dockey on Ubuntu systems. It may also support other Linux systems based on Debian, for example, Deepin.

#. Go to `https://github.com/lmdu/dockey/releases <https://github.com/lmdu/dockey/releases>`_ page, click on ``Dockey-version-amd64.deb`` to download it. Then double click the downloaded installer to install the program following the on-screen instructions.

#. Or, you can click on ``Dockey-version-linux.tar.gz`` to download the compressed package, and then uncompress it to run Dockey:

.. code:: shell

	tar xzvf Dockey-version-linux.tar.gz
	cd Dockey
	./Dockey

MacOS
-----

Because the docking engines and Dockey are unsigned applications which can not be installed from apple store, the gatekeeper of MacOS will prevent the installation and running of Dockey and AutoDock. So, before installation, you should disable gatekeeper from command line in MacOS.

To disable gatekeeper, follow these steps:

#. Launch **Terminal** from **Applications** > **Utilities**.
#. Enter the following command:

	.. code:: shell

		sudo spctl --master-disable

#. Press **Enter** and type your admin password.
#. Press **Enter** again.

Now, the Anywhere option should be available under the **Allow apps downloaded from** section of **System Preferences** > **Security & Privacy** > **General**.

.. note::
	If you want to re-enable gatekeeper, you can do with a simple command:

	.. code:: shell

		sudo spclt --master-enable


**AutoDock4:**

Go to `https://autodock.scripps.edu/download-autodock4/ <https://autodock.scripps.edu/download-autodock4/>`_ page, click on ``AutoDock4 (Mac OS X)`` to download ``autodocksuite-4.2.6-MacOSX.tar`` package, then uncompress the package:

.. code:: shell

	tar xvf autodocksuite-4.2.6-MacOSX.tar

After uncompress, you will get a folder named ``MacOSX`` which contains autogrid4 and autodock4 executables. You can move these two executables to any where, for example, move to /usr/local/bin:

.. code:: shell

	cd MacOSX
	sudo mv autogrid4 autodock4 /usr/local/bin

After installation, you should set the permissions so they can work following these commands:

.. code:: shell

	sudo xattr -r -d com.apple.quarantine /usr/local/bin/autogrid4
	sudo xattr -r -d com.apple.quarantine /usr/local/bin/autodock4

**AutoDock Vina:**

Go to `https://github.com/ccsb-scripps/AutoDock-Vina/releases <https://github.com/ccsb-scripps/AutoDock-Vina/releases>`_ page, click on ``vina_1.2.3_mac_x86_64`` to download it, and put it any where. And then set the permissions following:

.. code:: shell

	sudo xattr -r -d com.apple.quarantine vina_1.2.3_mac_x86_64

**QuickVina-W:**

The QuickVina homepage only provide linux executable. For MacOS, you can install QuickVina-W by `anaconda <https://anaconda.org/conda-forge/qvina>`_.

**Dockey:**

Go to `https://github.com/lmdu/dockey/releases <https://github.com/lmdu/dockey/releases>`_ page, click on ``Dockey-version-macos.dmg`` to download it. Then double click the downloaded installer to install the program following the on-screen instructions.

After installation, you should set the permissions following:

.. code:: shell

	sudo xattr -r -d com.apple.quarantine /Applications/Dockey.app


.. |windows| image:: _static/windows.png
	:width: 24
.. |linux| image:: _static/ubuntu.png
	:width: 24
.. |macos| image:: _static/apple.png
	:width: 24