Installation
============

Prior to installing Dockey, you should first install a docking engine, any of `AutoDock4 <https://autodock.scripps.edu/>`_, `AutoDock Vina <https://vina.scripps.edu/>`_ and `QuickVina-W <https://qvina.github.io/>`_.

Dockey download: `https://github.com/lmdu/dockey/releases <https://github.com/lmdu/dockey/releases>`_

中国镜像下载地址: `https://big.cdu.edu.cn/software/dockey <https://big.cdu.edu.cn/software/dockey>`_

We have tested Dockey on these systems:

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
     - 20.04, 22.04
   * - |macos|
     - MacOS
     - Big Sur 11.7

On Windows
----------

**Install AutoDock4:**

Go to `https://autodock.scripps.edu/download-autodock4/ <https://autodock.scripps.edu/download-autodock4/>`_ page, click on ``AutoDock4 Windows`` to download ``autodocksuite-4.2.6.i86Windows.exe`` installer, and then double click installer to install AutoDock4.

**Or install AutoDock Vina:**

Go to `https://github.com/ccsb-scripps/AutoDock-Vina/releases <https://github.com/ccsb-scripps/AutoDock-Vina/releases>`_ page, click on ``vina_1.2.5_win.exe`` to download it, and put it any where.

**Or install QuickVina-W:**

The QuickVina homepage only provide linux executable. For Windows, you can install QuickVina-W by `anaconda <https://anaconda.org/conda-forge/qvina>`_.

**Install Dockey:**

#. Go to `https://github.com/lmdu/dockey/releases <https://github.com/lmdu/dockey/releases>`_ page, click on ``Dockey-version-win64.exe`` to download it. Then double click the downloaded installer to install the program following the on-screen instructions.

#. Or, you can click on ``Dockey-version-win64.7z`` to download it. Then uncompress the 7z file using `7-zip <https://www.7-zip.org/>`_ and double click the Dockey.exe to run it.

On Linux
--------

**Install AutoDock4:**

Go to `https://autodock.scripps.edu/download-autodock4/ <https://autodock.scripps.edu/download-autodock4/>`_ page, click on ``AutoDock v4.2.6 (Linux 64bit)`` to download ``autodocksuite-4.2.6-x86_64Linux2.tar`` package, Then uncompress the package:

.. code:: shell

	tar xvf autodocksuite-4.2.6-x86_64Linux2.tar

After uncompress, you will get a folder named ``x86_64Linux2`` which contains autogrid4 and autodock4 executables. You can move these two executables to any where, for example, move to /usr/local/bin:

.. code:: shell

	cd x86_64Linux2
	sudo mv autogrid4 autodock4 /usr/local/bin

**Or install AutoDock Vina:**

Go to `https://github.com/ccsb-scripps/AutoDock-Vina/releases <https://github.com/ccsb-scripps/AutoDock-Vina/releases>`_ page, click on ``vina_1.2.5_linux_x86_64`` to download it, and put it any where.

**Or install QuickVina-W:**

Go to `https://qvina.github.io/ <https://qvina.github.io/>`_ page, and then click ``Download binary QuickVina-W`` button on the right side of page to download ``qvina-w`` executable, and put it any where.

**Install Dockey:**

.. note::

	Currently, we only tested Dockey on Ubuntu systems. It may also support other Linux systems based on Debian, for example, Deepin.

#. Go to `https://github.com/lmdu/dockey/releases <https://github.com/lmdu/dockey/releases>`_ page, click on ``Dockey-version-ubuntu.deb`` to download it. Then double click the downloaded installer to install the program following the on-screen instructions.

#. You can also install Dockey using command line tool like this:

	.. code:: shell

		sudo dpkg -i Dockey-version-ubuntu.deb

#. Or, you can click on ``Dockey-version-ubuntu.AppImage`` to download it and run it like this:

	.. code:: shell

		chmod +x Dockey-version-ubuntu.AppImage
		./Dockey-version-ubuntu.AppImage

	You can also double click the .AppImage file to run it.

#. Or, you can click on ``Dockey-version-ubuntu.tar.gz`` to download the compressed package, and then uncompress it to run Dockey:

.. code:: shell

	tar xzvf Dockey-version-ubuntu.tar.gz
	cd Dockey
	./Dockey

.. note::
	The installation file with ``-ubuntu20.04`` can only run on Ubuntu 20.04, while file with ``-ubuntu22.04`` can run on Ubuntu >= 22.04

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


**Install AutoDock4:**

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

**Or install AutoDock Vina:**

Go to `https://github.com/ccsb-scripps/AutoDock-Vina/releases <https://github.com/ccsb-scripps/AutoDock-Vina/releases>`_ page, click on ``vina_1.2.3_mac_x86_64`` to download it, and put it any where. And then set the permissions following:

.. code:: shell

	sudo xattr -r -d com.apple.quarantine vina_1.2.5_mac_x86_64

**Or install QuickVina-W:**

The QuickVina homepage only provide linux executable. For MacOS, you can install QuickVina-W by `anaconda <https://anaconda.org/conda-forge/qvina>`_.

**Install Dockey:**

Go to `https://github.com/lmdu/dockey/releases <https://github.com/lmdu/dockey/releases>`_ page, click on ``Dockey-version-macos.dmg`` to download it. Then double click the downloaded installer to install the program following the on-screen instructions.

After installation, you should set the permissions following:

.. code:: shell

	sudo xattr -r -d com.apple.quarantine /Applications/Dockey.app


.. |windows| image:: _static/windows.svg
	:width: 24
.. |linux| image:: _static/ubuntu.svg
	:width: 24
.. |macos| image:: _static/apple.svg
	:width: 24
