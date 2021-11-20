import sys

from PySide6.QtGui import *
from PySide6.QtCore import *
from PySide6.QtWidgets import *

def createIntroPage():
    page = QWizardPage();
    #page.setTitle("Introduction")

    label = QLabel("This wizard will help you register your copy of Super Product Two.")
    label.setWordWrap(True)

    layout = QVBoxLayout()
    layout.addWidget(label)
    page.setLayout(layout)

    return page

if __name__ == '__main__':
    app = QApplication(sys.argv)

    wizard = QWizard()
    wizard.setWizardStyle(QWizard.ClassicStyle)
    wizard.addPage(createIntroPage())

    wizard.setWindowTitle("Run AutoDock")
    wizard.show()

    sys.exit(app.exec())
