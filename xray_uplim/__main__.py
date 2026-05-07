"""
xray_uplim.__main__
--------------------
Default entry point — launched by the 'xray_uplim' console script or
by running:

    python -m xray_uplim

Tries to open the PySide6 GUI.  If PySide6 is not installed, falls back
gracefully to the command-line interface (same as 'xray_uplim-cli').
"""

import os
import sys


def _fix_qt_plugins():
    """
    On macOS (and sometimes Linux) with pip-installed PySide6, Qt cannot
    locate its own platform plugins (cocoa / xcb) because the search path
    is not set automatically.  We resolve this by pointing
    QT_QPA_PLATFORM_PLUGIN_PATH at the plugins/ directory that ships
    inside the PySide6 wheel before any Qt code is imported.
    """
    if 'QT_QPA_PLATFORM_PLUGIN_PATH' in os.environ:
        return   # already set by the user, leave it alone
    try:
        import PySide6
        pyside6_dir  = os.path.dirname(PySide6.__file__)
        # pip wheels bundle plugins at  PySide6/Qt/plugins/
        plugins_dir  = os.path.join(pyside6_dir, 'Qt', 'plugins')
        if os.path.isdir(plugins_dir):
            os.environ['QT_QPA_PLATFORM_PLUGIN_PATH'] = plugins_dir
    except Exception:
        pass   # nothing can be done; let Qt report the error normally


def main():
    _fix_qt_plugins()

    try:
        from PySide6.QtWidgets import QApplication
        from xray_uplim.gui.app import MainWindow

        app = QApplication(sys.argv)
        app.setApplicationName('xray_uplim')
        win = MainWindow()
        win.show()
        sys.exit(app.exec())

    except ImportError:
        print(
            "─────────────────────────────────────────────────────\n"
            "  PySide6 is not installed — falling back to CLI mode.\n"
            "  To enable the GUI:  pip install PySide6\n"
            "─────────────────────────────────────────────────────\n"
        )
        from xray_uplim.cli import main as cli_main
        cli_main()


if __name__ == '__main__':
    main()
