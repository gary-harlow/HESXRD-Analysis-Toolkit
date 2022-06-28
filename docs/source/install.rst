Quick install guide
===================

Before you can use HAT, you'll need to get it installed. 

Install Python
--------------

Being a Python application HAT requires python Python. 

Get the latest version of Python at https://www.python.org/downloads/ or with
your operating system's package manager. 

You should install at least python 3.8.

You can verify that Python is installed by typing ``python`` from your shell;
you should see something like::

    Python 3.8.y
    [GCC 4.x] on linux
    Type "help", "copyright", "credits" or "license" for more information.
    >>>

In windows we recommend you use something like Visual Studio Code or Anaconda


Install HAT
--------------

Install from PyPi The Python Package Index
````````````````

Typically in Linux or Windows if you use pip you can just type:

``python -m pip install –-user xray-hat``

And everything should be taken care of. Although it is recommonded you instead do this
inside of a dedicated virtual python enviroment to keep control of package versions. 

A nice guide can be found 'here <https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/#creating-a-virtual-environment>`_.

Install from source
````````````````

Download the source files from GitHub yourself 'here <https://github.com/gary-harlow/HESXRD-Analysis-Toolkit/archive/refs/heads/main.zip>'`_.

Then make sure you have all the required packages and versions that are listed in the pyproject.toml file. 

That's it!
==========

That's it -- you can now :doc:`move onto the tutorial </tutorial>`.