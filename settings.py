import os

_OUT_PATH = "."


def TMP_PATH():
    os.system("mkdir -p .mana")
    return ".mana"


def OUT_PATH():
    os.system("mkdir -p " + _OUT_PATH)
    return _OUT_PATH


def LOG_FILE():
    return OUT_PATH() + os.sep + "log.txt"


def REPORT_FILE():
    return OUT_PATH() + os.sep + "report.txt"
