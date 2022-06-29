#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
from PyQt6.QtWidgets import QApplication
from xrayhat import interface


def main():  
    app = QApplication(sys.argv)
    window = interface.MainWindow()
    window.show()

    app.exec()    

if __name__ == '__main__':
    main()
