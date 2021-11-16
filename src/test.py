from PySide6.QtWidgets import QTextBrowser, QApplication


if __name__ == '__main__':
    import sys

    app = QApplication(sys.argv)

    text_browser = QTextBrowser()
    str_html = """
        <!DOCTYPE html>
        <html>
        <body>

        <h1 style="color:blue;">Hello World!</h1>
        <p style="color:red;">Lorem ipsum dolor sit amet.</p>
        <input type="text" name="dd">

        </body>
        </html>
        """
    text_browser.setText(str_html)
    text_browser.show()
    text_browser.raise_()

    sys.exit(app.exec())